################################################################################
# STUDY 2: Effect of Discretisation Mismatch
# Author: Clara Panchaud
# 
# Description:
# This script simulates data on a fine grid from the UMove-SCR model and fits the model
# using multiple coarser grids to evaluate the impact of grid resolution
# on parameter recovery and computational efficiency.
################################################################################

# Load required libraries ------------------------------------------------------
library(ggplot2)
library(secr)
library(TMB)
library(expm)
library(tidyverse)
library(viridis)

# Source custom functions ------------------------------------------------------
source("Move-SCR/Functions/sim_MMPP.r")
source("Move-SCR/Functions/Extra_Func.r")

# Compile and load TMB model ---------------------------------------------------
compile("Move-SCR/Functions/like_MMPP.cpp")
dyn.load(dynlib("Move-SCR/Functions/like_MMPP"))

# Set random seed for reproducibility
set.seed(1)

# Load and prepare trap data ---------------------------------------------------

# Read trap locations
traps <- read.csv("Data/Study2/marten_traps2.csv")

# Create trap polygon
trap <- make.poly(x = traps$x, y = traps$y)

# Create fine simulation grid
mask_sim <- make.mask(trap, buffer = 2, spacing = 0.4, type = "trapbuffer")
S_sim <- nrow(mask_sim)

# Print grid size
print(paste("Simulation grid states (S):", S_sim))

# Convert traps to matrix format
traps_mat <- as.matrix(traps)

# Map traps to the state space
distance_traps_sim <- proxy::dist(traps_mat, mask_sim)
traps_loc_sim <- apply(distance_traps_sim, 1, which.min)

# Handle multiple cameras in same cell (does not happen with the large S chosen here, but code still includes the option)
camera_count_sim <- rep(1, S_sim)
trap_counts_sim <- table(traps_loc_sim)
camera_count_sim[as.numeric(names(trap_counts_sim))] <- as.numeric(trap_counts_sim)

# Set simulation parameters ----------------------------------------------------

# Study duration
Time <- 264  # hours (11 days)

# True population size
N <- 30

# Movement parameter (log scale)
q <- log(0.2)

# Detection rate parameter
l <- 0.14

# Create generator matrix Q for simulation -------------------------------------

# Get neighbour matrix
neighbour_sim <- neighbour_matrix(mask_sim)

# Build Q matrix with movement rates
Q_sim <- neighbour_sim * exp(q)

# Set diagonal elements (rates of leaving each state)
diag(Q_sim) <- -rowSums(Q_sim)

# Set up detection rates at each location
lambda_sim <- rep(0, S_sim)
lambda_sim[traps_loc_sim] <- l * camera_count_sim[traps_loc_sim]

# Define grid spacings and simulation settings ---------------------------------

# Different grid spacings to test the misspecification of state space
grid_spacings <- c(0.4, 0.5, 0.6, 1, 2)

# Number of simulations per grid
n_sim <- 50

# Initialize results data frame
total_rows <- n_sim * length(grid_spacings)
results <- data.frame(matrix(ncol = 12, nrow = total_rows))
colnames(results) <- c("sim_id", "grid_spacing", "S_fit", "n.observed", 
                       "q", "SE.q", "l", "SE.l", "N", "SE.N", 
                       "time", "detections")

# Run simulations --------------------------------------------------------------

row_idx <- 1

for (i in 1:n_sim) {
  cat("Simulation", i, "of", n_sim, "\n")
  
  # Sample initial locations for N individuals
  s0 <- sample(1:S_sim, N, rep = TRUE)
  
  # Simulate continuous-time Markov chain (movement paths)
  obj <- sim_CTMC(Q_sim, Time, s0, N)
  
  # Simulate detections using Markov-modulated Poisson process
  mmpp_original <- sim_MMPP(obj, lambda_sim)
  
  # Filter to observed individuals only
  observed_idx <- which(!vapply(mmpp_original$tt, is.null, logical(1)))
  mmpp_clean <- list(
    tt = mmpp_original$tt[observed_idx],
    ss = mmpp_original$ss[observed_idx]
  )
  observed_ind <- length(mmpp_clean$tt) - 1
  
  # Fit model with multiple grid resolutions ----------------------------------
  
  for (j in seq_along(grid_spacings)) {
    cat("  Fitting with spacing", grid_spacings[j], "\n")
    
    start_time <- Sys.time()
    
    # Create fitting grid
    mask_fit <- make.mask(trap, buffer = 2, spacing = grid_spacings[j], 
                          type = "trapbuffer")
    S_fit <- nrow(mask_fit)
    
    # Map traps to fitting grid
    distance_traps_fit <- proxy::dist(traps_mat, mask_fit)
    traps_loc_fit <- apply(distance_traps_fit, 1, which.min)
    
    # Get neighbour matrix for fitting grid
    neighbour_fit <- neighbour_matrix(mask_fit)
    
    # Map detections from simulation grid to fitting grid
    mapping <- setNames(traps_loc_fit, traps_loc_sim)
    mmpp_mapped <- mmpp_clean
    mmpp_mapped$ss <- lapply(mmpp_clean$ss, function(vec) {
      if (length(vec) == 0) return(NULL)
      mapped_values <- mapping[as.character(vec)]
      return(unname(mapped_values[!is.na(mapped_values)]))
    })
    
    # Handle multiple cameras at same location on fitting grid
    camera_count_fit <- rep(1, S_fit)
    trap_counts_fit <- table(traps_loc_fit)
    camera_count_fit[as.numeric(names(trap_counts_fit))] <- as.numeric(trap_counts_fit)
    
    # Construct Q and lambda templates
    lambda_fixed <- rep(1, S_fit)
    lambda_fixed[traps_loc_fit] <- 0  # Free parameters at active trap locations
    lambda_template <- rep(0, S_fit)
    lambda_template[traps_loc_fit] <- 2
    
    # Build TMB data list
    tmbdata <- list(
      Us = lengths(mmpp_mapped$tt),
      t = unlist(mmpp_mapped$tt),
      s = as.integer(unlist(mmpp_mapped$ss) - 1L),  # 0-based indexing for C++
      f = rep(1/S_fit, S_fit),  # Initial distribution (uniform)
      n_obs = observed_ind,
      t0 = 0,
      Time = Time,
      n_states = S_fit,
      n_indiv = length(mmpp_mapped$tt),
      Q_template = neighbour_fit,
      Q_fixed = 1 - neighbour_fit,
      lambda_template = lambda_template,
      lambda_fixed = lambda_fixed,
      camera_count = camera_count_fit
    )
    
    # Handle individuals with no detection history
    tmbdata$Us[unlist(lapply(mmpp_mapped$tt, \(x) all(is.na(x))))] <- 0
    
    # Fit MMPP model ---------------------------------------------------------
    
    # Initial parameter values
    theta_init <- rep(-2.5, 2)
    
    # Create autodiff function
    obj_tmb <- MakeADFun(data = tmbdata, parameters = list(theta = theta_init), 
                         DLL = "like_MMPP")
    
    tryCatch({
      # Perform optimization (suppress output)
      invisible(capture.output({
        fit <- nlminb(start = theta_init, objective = obj_tmb$fn, 
                      gradient = obj_tmb$gr)
      }))
      
      # Get parameter estimates and standard errors
      report <- sdreport(obj_tmb, par.fixed = fit$par)
      param <- summary(report)
      
      
      # Calculate population size estimate with confidence interval
      Pop_size <- confint_pop_mmmpp(fit$par, Time, observed_ind, S_fit, 
                                    neighbour_fit, traps_loc_fit, 
                                    rep(1/S_fit, S_fit), 
                                    report$cov, camera_count_fit,)
      
      

      
      end_time <- Sys.time()
      
      # Store results
      results[row_idx, "sim_id"] <- i
      results[row_idx, "grid_spacing"] <- grid_spacings[j]
      results[row_idx, "S_fit"] <- S_fit
      results[row_idx, "n.observed"] <- observed_ind
      results[row_idx, c("q", "SE.q")] <- param[3, 1:2]
      results[row_idx, c("l", "SE.l")] <- param[4, 1:2]
      results[row_idx, c("N", "SE.N")] <- Pop_size[1:2]
      results[row_idx, "time"] <- as.numeric(difftime(end_time, start_time, 
                                                      units = "secs"))
      results[row_idx, "detections"] <- sum(lengths(mmpp_mapped$tt))
      
    }, error = function(e) {
      cat("    Error in fitting:", e$message, "\n")
      # Fill with NAs for failed fits
      results[row_idx, "sim_id"] <- i
      results[row_idx, "grid_spacing"] <- grid_spacings[j]
      results[row_idx, "S_fit"] <- S_fit
      results[row_idx, "n.observed"] <- observed_ind
    })
    
    row_idx <- row_idx + 1
  }
}



# Calculate coverage and performance metrics -----------------------------------

# True population size

# Clean results (remove failed fits)
results_clean <- results[!is.na(results$N), ]


# Calculate comprehensive summary including AIC
coverage_results <- results_clean %>%
  group_by(grid_spacing, S_fit) %>%
  summarise(
    n_sims = n(),
    mean_N = mean(N, na.rm = TRUE),
    mean_SE_N = mean(SE.N, na.rm = TRUE),
    
    # Coverage calculation
    coverage_N = mean({
      CI_lower <- N - 1.96 * SE.N
      CI_upper <- N + 1.96 * SE.N
      (CI_lower <= true_N) & (CI_upper >= true_N)
    }, na.rm = TRUE) * 100,
    
    # Bias and RMSE
    bias_N = mean(N - true_N, na.rm = TRUE),
    rel_bias_N = (bias_N / true_N) * 100,
    rmse_N = sqrt(mean((N - true_N)^2, na.rm = TRUE)),
    
    # Other metrics
    mean_q = mean(q, na.rm = TRUE),
    mean_l = mean(l, na.rm = TRUE),
    mean_time = mean(time, na.rm = TRUE),
    .groups = "drop"
  )


# Print results
print("Coverage and Performance by Grid Spacing :")
print(coverage_results)

# Save results -----------------------------------------------------------------

# Save detailed results
#write.csv(results, "results_grid_comparison.csv", row.names = FALSE)

# Save summary results
#write.csv(coverage_results, "summary_grid_comparison.csv", row.names = FALSE)