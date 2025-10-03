# ==============================================================================
# OU SIMULATION STUDY
# By Clara Panchaud
# ==============================================================================
# This code simulates SCR data from an OU process and compares the results of fitting two 
# continuous time spatial capture-recapture models:
# 1. MSCR (Memory Spatial Capture-Recapture) - incorporates memory effects
# 2. CT SCR (Continuous Time SCR) - no memory
#
# The study evaluates how different values of the memory parameter tau in the OU movement model
# affect population estimation accuracy across these models.
# ==============================================================================

# Load required libraries
# ------------------------------------------------------------------------------
require(secr)
library(ggplot2)
library(dplyr)

# Source custom functions for likelihood calculations and simulation
# ------------------------------------------------------------------------------
Rcpp::sourceCpp("MSCR/Functions/LikelihoodC.cpp")
source("MSCR/Functions/Sim_Func.R")
source("MSCR/Functions/Fit_Func.R")

# set seed for reproducible results
set.seed(1)

# ==============================================================================
# STUDY SETUP AND PARAMETERS
# ==============================================================================

# Load the trap locations (from a real survey) and setup the spatial configuration
traps <- read.csv("Data/Study2/marten_traps2.csv", col.names = c("x","y"))
trap = make.poly(x=traps$x, y=traps$y)
trap <- trap[-31,]
cams <- read.traps(data = traps, detector = "count")


noccasions<-1584 # Number of capture occasions, one represents a 10 minutes period
T_days <- # Total survey times of 11 days
N<-20 # True total number of individuals 

# Create spatial meshes for activity center simulation and model fitting
# ------------------------------------------------------------------------------
# Fine mesh for sampling activity centers (high resolution)
samplemesh = make.mask(trap,buffer=2,spacing=0.05,type="trapbuffer" )
samplemask<-as.matrix(samplemesh)

# Coarser mesh for model fitting (computational efficiency)
mask = make.mask(trap,buffer=2,spacing=0.5,type="trapbuffer")
meshmat<-as.matrix(mask)

# Calculate study area
# ------------------------------------------------------------------------------
dim(mask)
a = attr(mask,"a") # Area of each cell of the mesh
D<-dim(mask)[1] # Number of mesh points
A = a*D # Total study area

# ==============================================================================
# SIMULATION LOOP SET UP
# ==============================================================================

# Define number of simulations per beta value
n_sim <- 20


tau <- c(0, 1, 5, 10, 50, 100, 200, 500) # memory parameter 1/beta in an OU process
                       # so the larger the less memory

sigma <- 0.5 # corresponds to movement per iteration if no memory 
epsilon <- 0.1 # km, how far from a trap an individual can be detected (100m )
p_capture <- 0.8 # How likely an individual is to be seen given it is within epsilon-km from a trap

theta_init<- c(1.00,  0.01 , 0.00) # Initial values for the optimisation

# Initialize results dataframe
# ------------------------------------------------------------------------------
results <- data.frame(matrix(ncol = 18, nrow = n_sim * length(tau)))
colnames(results) <- c(
  # MSCR model results (with memory)
  "h0", "sigma", "beta",                    # Parameter estimates
  "SE.h0", "SE.sigma", "SE.beta",          # Standard errors
  "N", "SE.N",                             # Population estimates
  
  # CT SCR model results (no memory)
  "h0.nomem", "sigma.nomem",               # Parameter estimates
  "SE.h0.nomem", "SE.sigma.nomem",         # Standard errors
  "N.nomem", "SE.N.nomem",                 # Population estimates
  
  # Additional information
  "n",                                     # Number of observed individuals
  "true.tau",                             # True beta value used
  "time.memory", "time.no.memory"  # Computation times
)

j<-1 # Results row counter 

# ==============================================================================
# MAIN SIMULATION LOOP
# ==============================================================================

# Loop over different tau values (memory parameter)
for (k in 1:length(tau)){
  
    tau_k <- tau[k] # set tau to the current value
  
  # Run multiple simulations for each value of beta
  for (i in 1:n_sim){

      # ----------------------------------------------------------------------
      # DATA SIMULATION
      # ----------------------------------------------------------------------
      
      # Sample activity centers for the N individuals
      sampling<-sample(nrow(samplemask),N,replace=TRUE)
      ac_coords<-as.matrix(samplemask[sampling,])
  
      # Get the right format to simulate the capture history, define the activity centres of the population as "pop".
      pop <- data.frame(x = ac_coords[,1], y = ac_coords[,2])
      class(pop) <- c("popn", "data.frame")
      attr(pop, "boundingbox") <- apply(pop[,1:2], 2, range)
      
      # Simulate paths for each individual of the population. 
      # This function also simulates capture histories but not in the way we want so we will only use the paths
      ch <- simOU.capthist(cams, pop, 
                           detectpar = list(tau = tau_k, sigma = sigma, epsilon = epsilon), 
                           noccasions = noccasions, 
                           savepath = TRUE)
      
      # Extract the paths from the previously simulated object
      paths <- attr(ch, "path")
    
      # Prepare the simulated data frame where each row is an observation, with variables:
      # id: unique individual ID
      # Time: time of capture (in days)
      # y: unique ID of the trap where the capture occurred
      data<-data.frame(matrix(ncol = 3, nrow = 0))
      names<-c("id","Time","y")
      colnames(data)<-names
      
      # From the simulated paths above, simulate the observations
      rows <- list()
      l <- 1
      for (i in 1:N) {
        path_i <- paths[[i]]
        for (t in 1:noccasions) {
          for (j_r in 1:nrow(cams)) {
            d <- sqrt((path_i[t,1] - cams$x[j_r])^2 + (path_i[t,2] - cams$y[j_r])^2)
            if (d <= epsilon && runif(1) < p_capture) {
              rows[[l]] <- data.frame(id = i, Time = t, y = j_r)
              l <- l + 1
            }
          }
        }
      }
      
      data <- do.call(rbind, rows)
      df_sim <- re_id(data) # Re-assign individual IDs since unobserved individuals were removed 
      
      # Add trap coordinates to the capture histories
      df_sim$trap_x<-traps[df_sim$y,1]
      df_sim$trap_y<-traps[df_sim$y,2]
      df_sim<-arrange(df_sim,id,Time)
      
      df_sim$Time <- df_sim$Time / 144 # We want the time in days, not minutes
      
      # ----------------------------------------------------------------------
      # DATA PREPARATION FOR MODEL FITTING
      # ----------------------------------------------------------------------

      # Define time discretisation for model fitting 
      r2<-20 
      n<-length(unique(df_sim$id)) # Number of observed individuals 
      
      # Discretise the capture histories for model fitting 
      n_rows <- sum(table(df_sim$id))  
      ddf_sim <- data.frame(
        t = numeric(n_rows),
        y = integer(n_rows),
        id = integer(n_rows)
      )
      
      row_idx <- 1
      for(i in unique(df_sim$id)){
        data <- discretize(df_sim[df_sim$id == i, ],T_days, r2)
        data$id <- i
        n_new <- nrow(data)
        ddf_sim[row_idx:(row_idx + n_new - 1), ] <- data  
        row_idx <- row_idx + n_new
      }
      # Convert to the format necessary for C++ functions
      ddfmat_sim = as.matrix(ddf_sim)
      dfrows_sim = as.numeric(table(ddf_sim$id))
      
      # ----------------------------------------------------------------------
      # MODEL 1: MSCR
      # ----------------------------------------------------------------------
      
      # Fit MSCR model with auxiliarly C++ function
      start_time <- Sys.time()
      fit <- optim(theta_init, LikelihoodC, trap = as.matrix(cams),df = ddfmat_sim, dfrows = dfrows_sim, mesh = meshmat, endt = T_days, hessian=TRUE)
      end_time <- Sys.time()
      time_memory<- end_time - start_time
      
      
      # Calculate population estimate and confidence interval 
      N_est<-confint_pop(fit,T_days,traps,mask,n)
      
      
      # Store MSCR results
      results[j,1:3] <- fit$par                         # Parameter estimates
      results[j,7:8] <- N_est[1:2]                      # Population estimate and SE
      results[j,4:6] <- sqrt(diag(solve(fit$hessian)))  # Standard errors
      
      
      # ----------------------------------------------------------------------
      # DATA FILTERING FOR NON-MEMORY MODELS
      # ----------------------------------------------------------------------
      
      # Filter data to have at most one detection per hour per individual
      # This removes the fine temporal structure that memory models use, and is common practice
      df_sim <- df_sim %>%
        arrange(id, Time) %>% 
        filter(!(id == lag(id) & (Time - lag(Time) < 0.04166667))) # 0.04166667 = 1/24 (1 hour)
      
      # Re-discretise filtered data
      n_rows <- sum(table(df_sim$id))
      ddf_sim <- data.frame(
        t = numeric(n_rows),
        y = integer(n_rows),
        id = integer(n_rows)
      )
      
      row_idx <- 1
      
      for(i in unique(df_sim$id)){
        data <- discretize(df_sim[df_sim$id == i, ],T_days, r2)
        data$id <- i
        n_new <- nrow(data)
        ddf_sim[row_idx:(row_idx + n_new - 1), ] <- data
        row_idx <- row_idx + n_new
        ddfmat_sim = as.matrix(ddf_sim)
        dfrows_sim = as.numeric(table(ddf_sim$id))
      }
      
      # ----------------------------------------------------------------------
      # MODEL 2: CT SCR
      # ----------------------------------------------------------------------
      
      start_time <- Sys.time()
      # Fit CT SCR model (only h0 and sigma, no memory parameter)
      fit_nomem <- optim(theta_init[1:2], LikelihoodCnoMem, trap = as.matrix(cams),df = ddfmat_sim, dfrows =dfrows_sim, mesh = meshmat, endt = T_days, hessian=TRUE)
      end_time <- Sys.time()
      time_nomemory<- end_time - start_time
      
      # Calculate population estimate
      N_est_nomem<-confint_pop(fit_nomem,T_days,trap,mask,n)
      
      # Store CT SCR results
      results[j, 9:10] <- fit_nomem$par                         # Parameter estimates
      results[j, 11:12] <- sqrt(diag(solve(fit_nomem$hessian))) # Standard errors
      results[j, 13:14] <- N_est_nomem[1:2]                     # Population estimate and SE
      
      # ----------------------------------------------------------------------
      # STORE ADDITIONAL INFORMATION
      # ----------------------------------------------------------------------
      
      
      results[j,15]<-n # Number of observed individuals
      results[j,16]<-tau[k] # True beta value
      results[j,17:18]<-c(time_memory,time_nomemory) # Computation times
      
      
      j <- j+1 # Move to next results row
  }
}

write.csv(results,"results_OU.csv")
results


library(dplyr)
library(ggplot2)
library(tidyr)

# Group by tau and calculate summary statistics
summary_by_tau <- results %>%
  group_by(true.tau) %>%
  summarise(
    # MSCR model summaries
    mean_N = mean(N, na.rm = TRUE),
    sd_N = sd(N, na.rm = TRUE),
    median_N = median(N, na.rm = TRUE),
    mean_h0 = mean(h0, na.rm = TRUE),
    mean_sigma = mean(sigma, na.rm = TRUE),
    mean_beta = mean(beta, na.rm = TRUE),
    
    # CT SCR model summaries
    mean_N_nomem = mean(N.nomem, na.rm = TRUE),
    sd_N_nomem = sd(N.nomem, na.rm = TRUE),
    median_N_nomem = median(N.nomem, na.rm = TRUE),
    
    # Other summaries
    mean_n_observed = mean(n, na.rm = TRUE),
    n_sims = n()
  )

print(summary_by_tau)

summary_by_tau

# Boxplot for N estimates by tau - MSCR model
ggplot(results, aes(x = factor(true.tau), y = N))  +
  geom_boxplot(fill = "lightblue") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Population Estimates by Tau (MSCR Model)",
       x = "True Tau (Memory Parameter)",
       y = "Estimated N") +
  theme_minimal()

# Boxplot for N estimates by tau - CT SCR model
ggplot(results, aes(x = factor(true.tau), y = N.nomem)) +
  geom_boxplot(fill = "lightgreen") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Population Estimates by Tau (CT SCR Model)",
       x = "True Tau (Memory Parameter)",
       y = "Estimated N") +
  theme_minimal()

# Side-by-side comparison
results_long <- results %>%
  select(true.tau, N, N.nomem) %>%
  pivot_longer(cols = c(N, N.nomem), 
               names_to = "Model", 
               values_to = "N_estimate") %>%
  mutate(Model = ifelse(Model == "N", "MSCR (with memory)", "CT SCR (no memory)"))

ggplot(results_long, aes(x = factor(true.tau), y = N_estimate, fill = Model)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Population Estimates by Tau: Model Comparison",
       x = "True Tau (Memory Parameter)",
       y = "Estimated N",
       fill = "Model Type") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen"))

# More detailed summary table
detailed_summary <- results %>%
  group_by(true.tau) %>%
  summarise(
    # MSCR
    MSCR_mean = mean(N, na.rm = TRUE),
    MSCR_sd = sd(N, na.rm = TRUE),
    MSCR_bias = mean(N - 20, na.rm = TRUE),
    MSCR_rmse = sqrt(mean((N - 20)^2, na.rm = TRUE)),
    
    # CT SCR
    CTSCR_mean = mean(N.nomem, na.rm = TRUE),
    CTSCR_sd = sd(N.nomem, na.rm = TRUE),
    CTSCR_bias = mean(N.nomem - 20, na.rm = TRUE),
    CTSCR_rmse = sqrt(mean((N.nomem - 20)^2, na.rm = TRUE)),
    
    n_sims = n()
  )

print(detailed_summary)






noccasions<-1584

samplemesh = make.mask(trap,buffer=2,spacing=0.05,type="trapbuffer" )
samplemask<-as.matrix(samplemesh)


sigma <- 0.5 # corresponds to movement per iteration if no memory 
epsilon <- 0.1 # km, how far from a trap an individual can be detected (100m )
p_capture <- 0.8 # How likely an individual is to be seen given it is within epsilon-km from a trap


tau_k <-0 # memory

sampling<-sample(nrow(samplemask),N,replace=TRUE)
ac_coords<-as.matrix(samplemask[sampling,])

# Get the right format to simulate the capture history, define the activity centres of the population as "pop".
pop <- data.frame(x = ac_coords[,1], y = ac_coords[,2])
class(pop) <- c("popn", "data.frame")
attr(pop, "boundingbox") <- apply(pop[,1:2], 2, range)

# Simulate paths for each individual of the population. 
# This function also simulates capture histories but not in the way we want so we will only use the paths
ch <- simOU.capthist(cams, pop, 
                     detectpar = list(tau = tau_k, sigma = sigma, epsilon = epsilon), 
                     noccasions = noccasions, 
                     savepath = TRUE)

# Extract the paths from the previously simulated object
paths <- attr(ch, "path")

# Prepare the simulated data frame where each row is an observation, with variables:
# id: unique individual ID
# Time: time of capture (in days)
# y: unique ID of the trap where the capture occurred
data<-data.frame(matrix(ncol = 3, nrow = 0))
names<-c("id","Time","y")
colnames(data)<-names

# From the simulated paths above, simulate the observations
rows <- list()
l <- 1
for (i in 1:N) {
  path_i <- paths[[i]]
  for (t in 1:noccasions) {
    for (j_r in 1:nrow(cams)) {
      d <- sqrt((path_i[t,1] - cams$x[j_r])^2 + (path_i[t,2] - cams$y[j_r])^2)
      if (d <= epsilon && runif(1) < p_capture) {
        rows[[l]] <- data.frame(id = i, Time = t, y = j_r)
        l <- l + 1
      }
    }
  }
}

data <- do.call(rbind, rows)
df_sim <- re_id(data)
df_sim





# Set up a 4x5 grid (20 panels)
par(mfrow = c(4, 5), mar = c(2,2,2,1), oma = c(4,4,4,2), pty = 's')

for (individual_id in 1:N) {
  path_i <- paths[[individual_id]]
  ind_data <- data[data$id == individual_id, ]
  trap_totals <- table(factor(ind_data$y, levels = 1:nrow(cams)))
  captured_traps <- which(trap_totals > 0)

  # Plot limits (slightly larger than min/max of path + AC + traps)
  all_x <- c(cams$x, pop$x[individual_id], path_i[,1])
  all_y <- c(cams$y, pop$y[individual_id], path_i[,2])
  xlim <- range(all_x) + c(-0.5, 0.5)
  ylim <- range(all_y) + c(-0.5, 0.5)

  plot(path_i[,1], path_i[,2], type = 'n',
       xlab = "", ylab = "",
       xlim = xlim, ylim = ylim, asp = 1,
       main = paste("Ind", individual_id), cex.main = 0.9)

  # Path
  lines(path_i[,1], path_i[,2], col = 'blue', lwd = 1)

  # Start
  points(path_i[1,1], path_i[1,2], pch = 16, col = 'green', cex = 1.2)

  # Activity center
  points(pop$x[individual_id], pop$y[individual_id], pch = 17, col = 'red', cex = 1.5)

  # Captured traps
  if(length(captured_traps) > 0){
    points(cams$x[captured_traps], cams$y[captured_traps],
           pch = 16, col = 'orange',
           cex = 1 + trap_totals[captured_traps] / max(trap_totals[captured_traps]) * 1.5)
  }

  # Optional: all traps as background
  points(cams$x, cams$y, pch = 18, col = 'maroon2', cex = 1)
}

# Add shared x/y labels
mtext("X", side = 1, outer = TRUE, line = 2.5)
mtext("Y", side = 2, outer = TRUE, line = 2.5)
