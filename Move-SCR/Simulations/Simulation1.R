################################################################################
# STUDY 1: PARAMETER RECOVERY AND COMPUTATIONAL SCALABILITY
# Author: Clara Panchaud
# 
# Description:
# This script simulates data from the UMove-SCR model and fits the model
# to the simulated data to ensure good parameter recovery and scalability when the state space increases.
################################################################################

# Load required libraries ------------------------------------------------------
library(ggplot2)
library(secr)
library(TMB)
library(expm)
library(tidyverse)

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

# Create mask (state space) with buffer and spacing - Change spacing to 1, 0.8 and 0.6 for the other values of S
mask <- make.mask(trap, buffer = 2, spacing = 1.5, type = "trapbuffer")

# Convert traps to matrix format
traps <- as.matrix(traps)

# Calculate distances between traps and mask points
distance_traps <- proxy::dist(traps, mask)

# Find nearest mask point for each trap
traps_loc <- apply(distance_traps, 1, which.min)

# Get number of states (mask points)
S <- dim(mask)[1]

# Count cameras at each location (handle multiple traps at same location)
camera_count <- rep(1, S)
trap_counts <- table(traps_loc)
camera_count[as.numeric(names(trap_counts))] <- as.numeric(trap_counts)

# Print number of states
print(paste("Number of states (S):", S))


# Set simulation parameters ----------------------------------------------------

set.seed(1)

# Movement parameter (log scale) - Modify as S changes according to Table 4.3 of thesis
q <- log(0.03) 

# Detection rate parameter - Modify as well with S
l <- 0.008

# Number of simulations
n_sim <- 50

# Initialize results data frame
results <- data.frame(matrix(ncol = 9, nrow = n_sim))
colnames(results) <- c("n.observed", "q", "SE.q", "l", "SE.l", "N", "SE.N", "time", 
                       "detections" )


# Set up detection rates at each location
lambda <- rep(0, S)
lambda[traps_loc] <- l * camera_count[traps_loc]

# Study duration
Time <- 264 # hours, represents 11 days

# True population size for the simulation study
N <- 30


# Create generator matrix Q ---------------------------------------------------

# Get neighbour matrix 
neighbour <- neighbour_matrix(mask)

# Build Q matrix with movement rates
Q <- neighbour * exp(q)

# Set diagonal elements (rates of leaving each state)
diag(Q) <- -rowSums(Q) 


# Run simulations --------------------------------------------------------------

sim <- 1

for (j in 1:n_sim) {
  
  # Sample initial locations for N individuals
  s0 <- sample(1:S, N, rep = TRUE)
  
  # Simulate continuous-time Markov chain (movement paths)
  obj <- sim_CTMC(Q, Time, s0, N)
  
  # Simulate detections using Markov-modulated Poisson process
  mmpp <- sim_MMPP(obj, lambda)
  
  # Filter to observed individuals only (those with at least one detection)
  observed_idx <- which(!vapply(mmpp$tt, is.null, logical(1)))
  
  # Subset detection times and locations
  mmpp$tt <- mmpp$tt[observed_idx]
  mmpp$ss <- mmpp$ss[observed_idx]
  
  # Count observed individuals
  observed_ind <- length(mmpp$tt) - 1
  
  # Total number of detections
  total_detections <- sum(lengths(mmpp$tt))
  
  cat("  Observed individuals:", observed_ind, "\n")
  cat("  Total detections:", total_detections, "\n")

  # Prepare data for MMPP model fitting ----------------------------------------
  
  # Construct Q and lambda templates
  Q_fixed <- 1 - neighbour         # Positions of structural zeros
  Q_template <- neighbour          # Positions to fill with transition rates
  
  lambda_fixed <- rep(1, S)
  lambda_fixed[traps_loc] <- 0     # Free parameters at active trap locations only
  lambda_template <- rep(0, S)
  lambda_template[traps_loc] <- 2
  
  # Initial distribution (uniform across all states)
  f <- rep(1/S, S)
  
  # Build TMB data list
  tmbdata <- list(
    Us = lengths(mmpp$tt),
    t = unlist(mmpp$tt),
    s = as.integer(unlist(mmpp$ss) - 1L),  # 0-based indexing for C++
    f = f,
    n_obs = observed_ind,
    t0 = 0,
    Time = Time,
    n_states = S,
    n_indiv = length(mmpp$tt),
    Q_template = Q_template,
    Q_fixed = Q_fixed,
    lambda_template = lambda_template,
    lambda_fixed = lambda_fixed,
    camera_count = camera_count
  )
  
  # Handle individuals with no detection history (set length to zero)
  tmbdata$Us[unlist(lapply(mmpp$tt, \(x) all(is.na(x))))] <- 0
  
  
  # Fit MMPP model -------------------------------------------------------------
  
  # Initial parameter values
  theta_init <- rep(-2.5, 2)
  
  # Record start time
  start <- Sys.time()
  
  # Create autodiff function
  obj <- MakeADFun(data = tmbdata, parameters = list(theta = theta_init), 
                   DLL = "like_MMPP")
  
  # Perform optimization (suppress output)
  invisible(capture.output({
    fit <- nlminb(start = theta_init, objective = obj$fn, gradient = obj$gr)
  }))
  
  # Get parameter estimates and standard errors
  report <- sdreport(obj, par.fixed = fit$par)
  param <- summary(report)
  
  # Record end time
  end <- Sys.time()
  
  # Calculate population size estimate with confidence interval
  Pop_size <- confint_pop_mmmpp(fit$par, Time, observed_ind, S, neighbour, 
                                traps_loc, f, report$cov)
  
  
  # Store MMPP results
  results[sim, 1] <- observed_ind
  results[sim, 2:3] <- param[3, 1:2]      # q estimate and SE
  results[sim, 4:5] <- param[4, 1:2]      # l estimate and SE
  results[sim, 6] <- Pop_size[1]          # N estimate
  results[sim, 7] <- Pop_size[2]          # SE of N
  results[sim, 8] <- as.numeric(difftime(end, start, units = "secs"))
  results[sim, 9] <- total_detections

  # Increment simulation counter
  sim <- sim + 1
}

# Save results -----------------------------------------------------------------

# Optionally save results to CSV
# write.csv(results, "simulation_results.csv", row.names = FALSE)


# Calculate bias metrics and store results -------------------------------------------------------

# True parameter values
true_q <- exp(q) 
true_l <- l  
true_N <- 30

# Calculate averages
avg_observed <- mean(results$n.observed, na.rm = TRUE)
avg_detections <- mean(results$detections, na.rm = TRUE)
avg_q <- mean(results$q, na.rm = TRUE)
avg_SE_q <- mean(results$SE.q, na.rm = TRUE)
avg_l <- mean(results$l, na.rm = TRUE)
avg_SE_l <- mean(results$SE.l, na.rm = TRUE)
avg_N <- mean(results$N, na.rm = TRUE)
avg_SE_N <- mean(results$SE.N, na.rm = TRUE)
avg_runtime <- mean(results$time, na.rm = TRUE)

# MMPP model bias
bias_N_MMPP <- mean(results$N, na.rm = TRUE) - 30
rel_bias_N_MMPP <- (bias_N_MMPP / 30) * 100

# Calculate coverage and confidence interval width -----------------------------

true_N <- 30

# MMPP model coverage
CI_lower_N_MMPP <- results$N - 1.96 * results$SE.N
CI_upper_N_MMPP <- results$N + 1.96 * results$SE.N
coverage_N_MMPP <- mean(CI_lower_N_MMPP <= true_N & CI_upper_N_MMPP >= true_N, 
                        na.rm = TRUE) * 100
CI_width_N_MMPP <- mean(CI_upper_N_MMPP - CI_lower_N_MMPP, na.rm = TRUE)


# Calculate RMSE for N
RMSE_N <- sqrt(mean((results$N - true_N)^2, na.rm = TRUE))


# Create summary table ---------------------------------------------------------

summary_table <- data.frame(
  Metric = c(
    "True value",
    "Avg. Estimate ± SE"
  ),
  log_alpha = c(
    sprintf("%.4f", true_q),
    sprintf("%.4f ± %.4f", avg_q, avg_SE_q)
  ),
  lambda = c(
    sprintf("%.4f", true_l),
    sprintf("%.4f ± %.4f", avg_l, avg_SE_l)
  ),
  N = c(
    sprintf("%.1f", true_N),
    sprintf("%.2f ± %.2f", avg_N, avg_SE_N)
  )
)

# Print summary results (for one column of Table 4.3 of thesis) -----------------------------------

cat("\n")
cat(strrep("=", 80), "\n")
cat("SIMULATION STUDY RESULTS (n =", n_sim, "simulations)\n")
cat(strrep("=", 80), "\n\n")

cat("Study Design:\n")
cat("  Number of states (S):", S, "\n")
cat("  Study duration:", Time, "hours\n")
cat("  True population size:", true_N, "\n\n")

cat("Detection Summary:\n")
cat("  Avg. observed individuals:", round(avg_observed, 2), "/", true_N, "\n")
cat("  Avg. total detections:", round(avg_detections, 2), "\n\n")

cat("Parameter Estimates:\n\n")
print(summary_table, row.names = FALSE, right = FALSE)

cat("\n")
cat("Population Size Performance Metrics:\n")
cat("  Relative Bias (%):", round(rel_bias_N_MMPP, 2), "%\n")
cat("  Coverage (%):", round(coverage_N_MMPP, 2), "%\n")
cat("  RMSE:", round(RMSE_N, 2), "\n")

cat("\n")
cat("Computational Performance:\n")
cat("  Average runtime:", round(avg_runtime, 2), "seconds\n")

cat("\n")
cat(strrep("=", 80), "\n\n")






