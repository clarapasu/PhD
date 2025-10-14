# ==============================================================================
# MARTEN POPULATION DENSITY ESTIMATION (SURVEY 2) USING ACMOVE-SCR MODEL
# ==============================================================================
# Author: Clara Panchaud
# Description: Analysis of the second survey of marten camera trap data using Markov-Modulated 
#              Poisson Process (MMPP) to estimate population density. AcMove-SCR model used here to model
#              the activity centre reverting movement. 
# ==============================================================================

# ------------------------------------------------------------------------------
#  SETUP AND DEPENDENCIES
# ------------------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(secr)
library(TMB)
library(expm)
library(tidyverse)
library(fields)
library(dplyr)
library(purrr)

# Source custom functions
source("MSCR/Functions/Fit_Func.r")
source("MSCR/Functions/Sim_Func.r")
source("Move-SCR/Functions/Extra_Func.r")

# Compile and load TMB model
compile("Move-SCR/Functions/like_MMPP_ac.cpp")
dyn.load(dynlib("Move-SCR/Functions/like_MMPP_ac"))

# Set random seed for reproducibility
set.seed(1)

# ------------------------------------------------------------------------------
# DEFINE HELPER FUNCTIONS
# ------------------------------------------------------------------------------

#' Create Neighbour Matrix for Mesh
#' Identifies neighbouring points in a spatial mesh based on minimum distances
#' 
#' @param mesh Matrix of spatial coordinates
#' @param tol Tolerance for considering points as neighbours (default: 1e-6)
#' @return Binary matrix where 1 indicates neighbouring relationship
neighbour_matrix <- function(mesh, tol = 1e-6) {
  
  # Calculate the pairwise distances between all mesh points
  dist_mat <- as.matrix(dist(mesh))
  
  # Set the diagonal to infinity to exclude self-distances
  diag(dist_mat) <- Inf
  
  # Find minimum distance for each point
  min_dist <- apply(dist_mat, 1, min)
  
  # Initialise neighbour matrix with zeros
  neighbour <- matrix(0, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
  
  # Mark neighbours: points within tolerance of minimum distance
  for (i in 1:nrow(dist_mat)) {
    neighbour[i, dist_mat[i, ] - min_dist[i] < tol] <- 1
  }
  
  return(neighbour)
}

# ------------------------------------------------------------------------------
#  LOAD AND PREPARE DATA
# ------------------------------------------------------------------------------

# Read the camera trap data and the trap locations
data <- read.csv("Data/Study2/marten_data2.csv")
traps <- read.csv("Data/Study2/marten_traps2.csv")

# Create a trap polygon
trap <- make.poly(x = traps$x, y = traps$y)

# ------------------------------------------------------------------------------
# CREATE DUAL MESH SYSTEM
# ------------------------------------------------------------------------------

# === DUAL MESH APPROACH ===
# Create TWO different meshes with different purposes:

# 1. FINE MESH: For state space (animal movement between locations)
#    This represents the spatial resolution of the biological process
state_spacing <- 0.5  # Fine spacing for detailed spatial process 
mask_states <- make.mask(trap, buffer = 2, spacing = state_spacing, type = "trapbuffer")

# Cell spacing (km)
h <- attr(mask_states, "spacing")

# Area per cell (km²)
a <- h * h

# Total study area (km²)
A <- nrow(mask_states) * a

# 2. COARSE MESH: For activity center integration (numerical integration only)
#    This reduces computational load while maintaining accuracy
ac_spacing <- 1.5  # Coarser spacing to reduce computational burden 
mask_ac <- make.mask(trap, buffer = 2, spacing = ac_spacing, type = "trapbuffer")


# ------------------------------------------------------------------------------
# CALCULATE SPATIAL DISTANCES
# ------------------------------------------------------------------------------

# === DISTANCE CALCULATIONS ===
# CRITICAL: Use the same distance scaling for both meshes to maintain 
#           parameter interpretability across the dual mesh system

# Convert traps to matrix format
traps <- as.matrix(traps)

# Calculate raw distances for both meshes
mesh_dist_raw_states <- rdist(mask_states, mask_states)
mesh_dist_raw_ac <- rdist(mask_ac, mask_ac)

# Use the state space mesh to define the distance scale (the "biological" scale)
mean_dist <- mean(mesh_dist_raw_states[mesh_dist_raw_states > 0])
cat("Distance scale (from state mesh):", mean_dist, "\n")

# Scale both distance matrices using the same scaling factor
mesh_dist_states <- mesh_dist_raw_states / mean_dist  # For spatial transitions
mesh_dist_ac <- mesh_dist_raw_ac / mean_dist         # For AC integration

# Distance from AC mesh to state mesh (for mapping between meshes)
ac_to_states_dist <- rdist(mask_ac, mask_states) / mean_dist

# ------------------------------------------------------------------------------
# MAP TRAPS TO STATE SPACE
# ------------------------------------------------------------------------------

# === TRAP LOCATION MAPPING ===
# Map traps to the state space mesh (where animals are detected)
distance_traps <- proxy::dist(traps, mask_states)
traps_loc <- apply(distance_traps, 1, which.min)

# Get mesh dimensions
S <- dim(mask_states)[1]  # Number of states = state space mesh size
S_AC <- dim(mask_ac)[1]   # Number of AC integration points

# Remap trap identifiers in data to mask locations
data$y <- traps_loc[data$y]

# Set observation time period
Time <- 11

# Build neighbour matrix for STATE SPACE mesh
neighbour <- neighbour_matrix(mask_states)

# ------------------------------------------------------------------------------
# PREPARE CAPTURE HISTORY DATA
# ------------------------------------------------------------------------------

# Group the data by individual ID and sort by time
df_grouped <- data %>% 
  group_by(id) %>% 
  arrange(Time, .by_group = TRUE)

# Extract the capture histories (spatial locations and times) for each individual
mmpp <- list(
  ss = df_grouped %>% group_split() %>% map(~ .x$y),
  tt = df_grouped %>% group_split() %>% map(~ .x$Time)
)

# Add an NA at the end as a convention
mmpp$ss[[length(mmpp$ss) + 1]] <- NA
mmpp$tt[[length(mmpp$tt) + 1]] <- NA

# Number of observed individuals
observed_ind <- length(mmpp$ss) - 1

# ------------------------------------------------------------------------------
# CONSTRUCT MODEL MATRICES
# ------------------------------------------------------------------------------

# Lambda (observation rate) structure:
# - Fixed at 1 for non-trap locations (no detections possible)
# - Free parameters at trap locations
lambda_fixed <- rep(1, S)
lambda_fixed[traps_loc] <- 0

lambda_template <- rep(0, S)
lambda_template[traps_loc] <- 1

# Count cameras at each location (handle multiple traps in same cell)
camera_count <- rep(1, S)
trap_counts <- table(traps_loc)
camera_count[as.numeric(names(trap_counts))] <- as.numeric(trap_counts)

# Q matrix structure:
# - Fixed at 1 for non-neighbours (structural zeros)
# - Template has 2 for neighbours (to be filled with transition rates)
Q_fixed <- 1 - neighbour
Q_template <- 2 * neighbour

# ------------------------------------------------------------------------------
# BUILD TMB DATA OBJECT WITH DUAL MESH INFORMATION
# ------------------------------------------------------------------------------

tmbdata <- list()

# Individual capture history data
tmbdata$Us <- as.integer(unlist(lapply(mmpp$tt, length)))
tmbdata$t <- unlist(mmpp$tt)
tmbdata$s <- as.integer(unlist(mmpp$ss) - 1)  # 0-indexed for C++
tmbdata$Us[unlist(lapply(mmpp$tt, \(x) all(is.na(x))))] <- 0

# Time and basic parameters
tmbdata$t0 <- 0
tmbdata$Time <- Time
tmbdata$n_states <- S        # State space mesh size
tmbdata$n_ac <- S_AC         # Activity center mesh size
tmbdata$n_indiv <- length(mmpp$tt)
tmbdata$n_obs <- observed_ind

# Area calculations
a_states <- attr(mask_states, "a")
a_ac <- attr(mask_ac, "a")
A_total <- a_states * S  # Total area

tmbdata$area_cell_ac <- a_ac
tmbdata$area_region <- A_total

# Mesh distance matrices
tmbdata$mesh_dist_states <- mesh_dist_states      # S x S: for spatial transitions
tmbdata$mesh_dist_ac <- mesh_dist_ac             # S_AC x S_AC: for AC-specific Q matrices
tmbdata$ac_to_states_dist <- ac_to_states_dist   # S_AC x S: mapping AC to states

# Q and lambda templates (based on state space mesh)
tmbdata$Q_template <- Q_template
tmbdata$Q_fixed <- Q_fixed
tmbdata$lambda_template <- lambda_template
tmbdata$lambda_fixed <- lambda_fixed
tmbdata$camera_count <- camera_count

# ------------------------------------------------------------------------------
#  MODEL FITTING
# ------------------------------------------------------------------------------

# Initial parameter values: (log_lambda, theta1, theta2)
theta_init <- c(0.1, -1, 0.01)

# Record start time
start <- Sys.time()

# Create the TMB autodiff function
obj <- MakeADFun(
  data = tmbdata, 
  parameters = list(theta = theta_init), 
  DLL = "like_MMPP_ac", 
  silent = FALSE
)

# Optimise the model
out <- nlminb(theta_init, obj$fn, obj$gr)

# Record end time and calculate duration
end <- Sys.time()
optimization_time <- difftime(end, start, units = "hours")

# ------------------------------------------------------------------------------
# COMPUTE STANDARD ERRORS
# ------------------------------------------------------------------------------

# Get standard errors using sdreport
tryCatch({
  rep <- sdreport(obj, par.fixed = out$par)
}, error = function(e) {
  cat("Standard error computation failed:", e$message, "\n")
})

# Calculate AIC
nll <- obj$fn(out$par)
AIC <- 2 * 3 + nll * 2  
AIC
---------------------------------------------------------------------------
# COMPUTE CONFIDENCE INTERVALS
# ------------------------------------------------------------------------------

# Calculate population size and confidence intervals using custom function
Pop <- confint_pop_mmmpp_ac(
  out$par, 
  Time, 
  observed_ind, 
  S, 
  neighbour, 
  traps_loc, 
  camera_count, 
  rep$cov, 
  S_AC, 
  ac_to_states_dist)


# Calculate density estimates (convert to per 100 ha)
D <- Pop[1] / A          # Density (animals/100 ha)
SE.D <- Pop[2] / A       # Standard error of density
CI_lower_D <- D - 1.96 * SE.D       # Lower 95% CI
CI_upper_D <- D + 1.96 * SE.D       # Upper 95% CI

# ------------------------------------------------------------------------------
# DISPLAY RESULTS
# ------------------------------------------------------------------------------

# Extract parameter estimates with confidence intervals
param_estimates <- summary(rep)

param_1 <- exp(param_estimates[1, "Estimate"])
param_1_SE <- param_estimates[1, "Std. Error"] * param_1
param_1_CI_lower <- param_1 - 1.96 * param_1_SE
param_1_CI_upper <- param_1 + 1.96 * param_1_SE

param_2 <- param_estimates[2, "Estimate"]
param_2_SE <- param_estimates[2, "Std. Error"]
param_2_CI_lower <- param_2 - 1.96 * param_2_SE
param_2_CI_upper <- param_2 + 1.96 * param_2_SE

param_3 <- param_estimates[3, "Estimate"]
param_3_SE <- param_estimates[3, "Std. Error"]
param_3_CI_lower <- param_3 - 1.96 * param_3_SE
param_3_CI_upper <- param_3 + 1.96 * param_3_SE

# Create summary table
results_table <- data.frame(
  Parameter = c("Density (animals/km²)", "lambda", "theta1", "theta2"),
  Estimate = c(round(D, 4), round(param_1, 4), round(param_2, 4), round(param_3, 4)),
  SE = c(round(SE.D, 4), round(param_1_SE, 4), round(param_2_SE, 4), round(param_3_SE, 4)),
  CI_Lower = c(round(CI_lower_D, 4), round(param_1_CI_lower, 4), 
               round(param_2_CI_lower, 4), round(param_3_CI_lower, 4)),
  CI_Upper = c(round(CI_upper_D, 4), round(param_1_CI_upper, 4), 
               round(param_2_CI_upper, 4), round(param_3_CI_upper, 4)),
  stringsAsFactors = FALSE
)



# Print the results
cat("=== FINAL RESULTS ===\n")
print(results_table, row.names = FALSE)

# Display density-specific results
cat("\n=== DENSITY ESTIMATES ===\n")
cat("Density:", round(D, 4), "animals/100 ha\n")
cat("Standard Error:", round(SE.D, 4), "\n")
cat("95% CI: [", round(CI_lower_D, 4), ",", round(CI_upper_D, 4), "]\n")

# Display optimisation results
cat("\n=== OPTIMISATION TIME ===\n")
cat("Time:", round(optimization_time, 2), "hours\n")

# Display AIC
cat("AIC:", AIC, "\n\n")