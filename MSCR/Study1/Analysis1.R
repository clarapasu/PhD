# ==============================================================================
# ANALYSIS OF PINE MARTEN DATA SET (SURVEY 1)
# By Clara Panchaud
# ==============================================================================
# This analysis fits the data of the Survey 1 of American martens with two different models:
# 1. MSCR (Memory Spatial Capture-Recapture) - incorporates memory effects
# 2. CT SCR (Continuous Time SCR) - no memory, continuous time
#
# The results make up Table 3.3 of the thesis. 
# ==============================================================================

# Load required libraries
# ------------------------------------------------------------------------------
library(lubridate)
library(ggplot2)
library(dplyr)
library(tibble)
library(secr)
library(sf)

# Source custom functions for likelihood calculations
# ------------------------------------------------------------------------------
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Fit_Func.R")
source("Functions/Sim_Func.R") # Needed for the half normal function

# ==============================================================================
# LOADING THE DATA
# ==============================================================================

# load the American marten data
# - traps consists of the UTM coordinates (in km) of the camera traps. The ID of a trap is its row number.
# - in data we have the capture histories, where each row represents an observation and 
#   consists of three variables (id, y, Time)
#   representing an individual's unique ID, the index of the camera trap
#   where the capture occurred, and the time of the capture.
data<-read.csv("Study1/Data/marten_data1.csv")
traps<-read.csv("Study1/Data/marten_traps1.csv")


# Create spatial meshes for model fitting
# ------------------------------------------------------------------------------
trap = make.poly(x=traps$x, y=traps$y) # Units are in km 
trap <- trap[-31,] # Remove the last element as it was repeated 
mask = make.mask(trap,buffer=2,spacing=0.2,type="trapbuffer") # Resolution determined in the file time_tests.R, spacing of 0.2 km
sp <- attr(mask, "spacing") 
a <- sp^2 # Area of one cell of the mask (0.04 km^2)
A = a*dim(mask)[1] # Total area covered by the mask (100.12 km^2)
meshmat<-as.matrix(mask)
traps<-as.matrix(traps)


# ==============================================================================
# PLOTTING THE LANDSCAPE
# Code to uncomment to produce Figure 1.1(a) of the thesis 
# ==============================================================================


# mask_plot = make.mask(trap, buffer = 2, spacing = 0.05, type = "trapbuffer")
# 
# mask_edge <- mask_plot %>%
#   filter(!duplicated(paste(x, y)))  
# 
# mask_hull <- mask_edge[chull(mask_edge$x, mask_edge$y), ]
# 
# trap_plot <- trap %>%
#   mutate(Type = "Trap")
# ggplot() +
#   geom_polygon(data = mask_hull, aes(x = x, y = y),
#                fill = "grey80", color = "grey30", size = 1) +
#   geom_point(data = trap_plot, aes(x = x, y = y, shape = Type, color = Type),
#              size = 4) +
#   scale_shape_manual(values = 18) +  
#   scale_color_manual(values = "grey10") +
#   coord_fixed(
#     xlim = c(310, 320.5),
#     ylim = c(4950, 4970)
#   ) +
#   labs(
#     x = "Easting (km)",
#     y = "Northing (km)",
#     shape = "",    
#     color = ""    
#   ) +
#   theme( panel.background = element_rect(fill = "white"),
#          panel.grid.major = element_line(color = "grey80", size = 0.3), 
#          panel.grid.minor = element_line(color = "grey90", size = 0.2),
#          legend.position = c(-0.02, 0.92),     
#          legend.justification = c("left","top"),
#          legend.background = element_rect(fill = alpha("white", 0.6), color = "black"),
#          legend.key = element_blank(),
#          legend.spacing = unit(0.1, "cm"),     
#          legend.margin = margin(-10, 3, 1, 2),    
#          legend.box.margin = margin(0, 0, 0, 0)
#   )
# 
# ggsave(
#   filename = "Landscape1.pdf",
#   plot = last_plot(),     
#   width = 6,
#   height = 5
# )

# ==============================================================================
# SET UP FOR MODEL FITTING 
# ==============================================================================

T_days<-12                 # Known survey duration (days)
n<-length(unique(data$id)) # Number of observed individuals
L<-100                     # Number of bins in the time discretisation 

# Discretise the capture histories for model fitting
ddf <- data.frame(t = as.numeric(), y = as.integer(), id = as.integer())
for(i in unique(data$id)){
  df <- discretize(data[data$id == i, ],T_days, L)
  df$id <- i
  ddf <- rbind(ddf, df)
}
# Convert to the format necessary for C++ functions
ddfmat = as.matrix(ddf)
dfrows = as.numeric(table(ddf$id))

# ----------------------------------------------------------------------
# MODEL 1: MSCR
# ----------------------------------------------------------------------

# Inital parameters values for the optimisation
theta_init<-c(1.5, -2 ,1.5) 

# Fit MSCR model with auxiliarly C++ function
start_time <- Sys.time()
fit <- optim(theta_init, LikelihoodC, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T_days,hessian=TRUE) 
end_time <- Sys.time()
fit_time<-end_time - start_time

# Extract the estimated parameters and confidence intervals
theta_est<-fit$par
param<-confint_param(fit)
N_est<-confint_pop(fit,T_days,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
D_est <- N_est/A
rownames(D_est) <- "D"

# ----------------------------------------------------------------------
# MODEL 2: CT SCT
# ----------------------------------------------------------------------

# Fit CT SCR model with auxiliarly C++ function
start_time <- Sys.time()
fit_nomem <- optim(theta_init[1:2], LikelihoodCnoMem, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T_days,hessian=TRUE) 
fit_time_nomem<-end_time <- Sys.time()
time_nomem<- end_time - start_time

# Extract the estimated parameters and confidence intervals
param_nm<-confint_param(fit_nomem)
N_nm<-confint_pop(fit_nomem,T_days,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
D_nm <- N_nm/A
rownames(D_nm) <- "D"


# ==============================================================================
# RESULTS EXTRACTION AND PRESENTATION
# ==============================================================================

# ----------------------------------------------------------------------
# PARAMETER ESTIMATES
# ----------------------------------------------------------------------

# Function to extract parameter estimates and CIs from parameters on the log scale
extract_real_scale_params <- function(fit_object, param_names) {
  
  # Get parameter estimates (log scale)
  log_estimates <- fit_object$par
  
  # Transform to real scale
  real_estimates <- exp(log_estimates)
  
  # Calculate standard errors using delta method
  # SE on real scale = exp(log_estimate) * SE_log_scale
    vcov_matrix <- solve(fit_object$hessian)  # Get variance-covariance matrix
    se_log <- sqrt(diag(vcov_matrix))
    se_real <- exp(log_estimates) * se_log
    
    # Calculate 95% CI on real scale
    # Transform log-scale CI bounds
    log_lower <- log_estimates - 1.96 * se_log
    log_upper <- log_estimates + 1.96 * se_log
    
    ci_lower <- exp(log_lower)
    ci_upper <- exp(log_upper)
    
    # Create results dataframe
    results_df <- data.frame(
      Parameter = param_names,
      Estimate = real_estimates,
      SE = se_real,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper
    )
    return(results_df)
}

# Extract parameters for MSCR model
mscr_params <- extract_real_scale_params(fit, c("h0", "sigma", "beta"))

# Extract parameters for CT-SCR model
ct_scr_params <- extract_real_scale_params(fit_nomem, c("h0", "sigma"))

# ----------------------------------------------------------------------
# AIC comparison
# ----------------------------------------------------------------------

# Calculate AIC values for both models (fit$value contains the negative log-likelihood)
AIC_MSCR <- 2 * fit$value + 2 * 3      # 3 parameters: h0, sigma, beta
AIC_CT_SCR <- 2 * fit_nomem$value + 2 * 2  # 2 parameters: h0, sigma

# Calculate delta AIC (difference from best model)
delta_AIC <- AIC_MSCR - AIC_CT_SCR


# ----------------------------------------------------------------------
# Printing results
# ----------------------------------------------------------------------

cat("Density estimate of MSCR: \n")
print(D_est)

cat("Density estimate of CT SCR: \n")
print(D_nm)

cat("Parameter estimates of MSCR: \n")
print(mscr_params, digits = 4)

cat("Parameter estimates of CT SCR: \n")
print(ct_scr_params, digits = 4)
    
cat("Delta AIC (MSCR - CT-SCR):", round(delta_AIC, 2), "\n\n")

