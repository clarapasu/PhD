# ==============================================================================
# ANALYSIS OF PINE MARTEN DATA SET (SURVEY 1)
# By Clara Panchaud
# ==============================================================================
# This analysis fits the data of the Survey 1 of American martens with two different models:
# 1. MSCR (Memory Spatial Capture-Recapture) - incorporates memory effects
# 2. CT SCR (Continuous Time SCR) - no memory, continuous time
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


# ==============================================================================
# LOADING THE DATA
# ==============================================================================

# load the American marten data
# - data collection started on the 27th of February 2017
# - traps consists of UTM coordinates, the ID of a trap is its row number
# - in data we have the capture histories consisting of three variables 
#   representing an individual's unique ID, the time of the capture
#   and the ID of the camera trap where the capture occurred.
data<-read.csv("Study1/marten_data1.csv")
traps<-read.csv("Study1/marten_traps1.csv")


# Create spatial meshes for model fitting
# ------------------------------------------------------------------------------
trap = make.poly(x=traps$x, y=traps$y)
trap <- trap[-31,] # Remove the last element as it was repeated 
mask = make.mask(trap,buffer=2,spacing=0.2,type="trapbuffer") # Resolution determined in the file time_tests.R
meshmat<-as.matrix(mask)
traps<-as.matrix(traps)

# ==============================================================================
# PLOTTING THE LANDSCAPE
# Code to uncomment to produce Figure 1.1(a) of the thesis 
# ==============================================================================

# mask_plot = make.mask(trap, buffer = 2, spacing = 0.01, type = "trapbuffer")
# plot(mask_plot, dots = FALSE, border = 1, ppoly = FALSE, asp = 1, 
#     xlab = "Longitude", ylab = "Latitude", main = "Study Area")
# plotMaskEdge(mask_plot, add = TRUE, col = "black", lwd = 1.5) 
# points(trap, pch = 4, col = "black", cex = 1.5)
# title(main = "Study Area with Trap Locations", col.main = "black", font.main = 2, cex.main = 1.5)


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
param<-confint_param(fit,T_days,traps,mask)
N_est<-confint_pop(fit,T_days,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)


# ----------------------------------------------------------------------
# MODEL 2: CT SCT
# ----------------------------------------------------------------------

# Fit CT SCR model with auxiliarly C++ function
start_time <- Sys.time()
fit_nomem <- optim(theta_init[1:2], LikelihoodCnoMem, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T_days,hessian=TRUE) 
fit_time_nomem<-end_time <- Sys.time()
time_nomem<- end_time - start_time

# Extract the estimated parameters and confidence intervals
param_nm<-confint_param(fit_nomem,T_days,traps,mask)
N_nm<-confint_pop(fit_nomem,T_days,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)

# ----------------------------------------------------------------------
# AIC comparison
# ----------------------------------------------------------------------

# Calculate the difference in AIC
theta_est_nomem<-fit_nomem$par
(2*fit_nomem$value-2*2)-(2*fit$value - 2*3)
