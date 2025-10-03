# ==============================================================================
# DISCRETISATIONS TESTS FOR SURVEY 1
# By Clara Panchaud
# ==============================================================================
# This code tests different time and space discretisations until finding
# the stable values used in Analysis2.R. The results make up Table 3.4 of the thesis. 
#
# The purpose is to conduct a sensitivity analysis on the spatial (mask spacing)
# and temporal (parameter L) settings for MSCR and SCR models to ensure parameter stability. 
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
Rcpp::sourceCpp("MSCR/Functions/LikelihoodC.cpp")
source("MSCR/Functions/Sim_Func.R")
source("MSCR/Functions/Fit_Func.R")

# ==============================================================================
# LOADING THE DATA
# ==============================================================================

# Load trap locations and capture data
data<-read.csv("Data/Study2/marten_data2.csv")
traps_og<-read.csv("Data/Study2/marten_traps2.csv")

# ==============================================================================
# SENSITIVITY ANALYSIS SETUP
# ==============================================================================

# Define the parameter ranges for the sensitivity analysis
# Test different spatial mask spacing values (in km)
space_list <- c(0.1, 0.2,0.5,1)

# Test different temporal discretization parameters 
# L is number of time intervals for discretisation, so that
# the hazard function of MSCR is constant over intervals of length T_days/L
L_list <- c(10,50,100,200)

# Initialise the results data frame to store all combinations
# 4 spacing × 4 L values = 16 total combinations
results<-data.frame(matrix(ncol = 8, nrow = 16))
colnames(results)<-c( "spacing", "L", "time_MSCR", "time_CTSCR", "N_MSCR", "N_CTSCR", "Likelihood_MSCR", "Likelihood_CTSCR")

# Counter for filling the results data frame
j<-1

# ==============================================================================
# MAIN SENSITIVITY ANALYSIS LOOP
# ==============================================================================

# Nested loop to test all combinations of spatial and temporal parameters
for (space in space_list){
  for (L in L_list){
    
    # ----------------------------------------------------------------------
    # SPATIAL SETUP
    # ----------------------------------------------------------------------
    
    # Create the spatial mesh and the right format of traps
    trap = make.poly(x=traps_og$x, y=traps_og$y)
    trap <- trap[-31,]
    mask = make.mask(trap,buffer=2,spacing=space,type="trapbuffer") # Use the current mask spacing
    meshmat<-as.matrix(mask)
    traps<-as.matrix(trap)
    dim(meshmat)
    
    T_days<-11 # Survey length (days)
    n<-length(unique(data$id))  # Number of observed individuals
    
    # Discretise each individual capture history for model fitting, using the time discretisation by L
    ddf <- data.frame(t = as.numeric(), y = as.integer(), id = as.integer())
    for(i in unique(data$id)){
      df <- discretize(data[data$id == i, ],T_days, L)
      df$id <- i
      ddf <- rbind(ddf, df)
    }
    ddfmat = as.matrix(ddf)
    dfrows = as.numeric(table(ddf$id))
    
    # ----------------------------------------------------------------------
    # MODEL FITTING: MSCR
    # ----------------------------------------------------------------------
    
    # Inital parameters values for the optimisation
    theta_init<-c(1.5, -2 ,1.5) 
    
    # Fit MSCR model and record computation time
    start_time <- Sys.time()
    fit <- optim(theta_init, LikelihoodC, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T_days,hessian=TRUE) 
    end_time <- Sys.time()
    fit_time<-end_time - start_time
    
    # Extract the results (negative log-likelihood value at the optimum and population estimate)
    Likelihood.m <- LikelihoodC(fit$par,traps,ddfmat, dfrows, meshmat, T_days)
    N_est<-confint_pop(fit,T_days,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
    
    # ----------------------------------------------------------------------
    # MODEL FITTING: CT SCR 
    # ----------------------------------------------------------------------
    
    ## Fit CT SCR model and record computation time
    start_time <- Sys.time()
    fit_nomem <- optim(theta_init[1:2], LikelihoodCnoMem, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T_days,hessian=TRUE) 
    fit_time_nomem<-end_time <- Sys.time()
    time_nomem<- end_time - start_time
    
    # Extract the results (negative log-likelihood value at the optimum and population estimate)
    N_nm<-confint_pop(fit_nomem,T_days,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
    Likelihood.nm <- LikelihoodCnoMem(fit_nomem$par,traps, ddfmat, dfrows, meshmat, T_days)
    
    # ----------------------------------------------------------------------
    # STORE RESULTS
    # ----------------------------------------------------------------------
    
    # Convert computation times to minutes
    fit_time<-as.numeric(fit_time, units="mins")
    time_nomem<-as.numeric(time_nomem, units="mins")
    
    # Store the results for this parameter combination
    results[j,]<-c(space, L, fit_time, time_nomem, N_est[1], N_nm[1], Likelihood.m, Likelihood.nm)
    j <- j+1
  }
  }

# Display the results (Extract time_MSCR and Likelihood_MSCR to obtain Table 3.4)
print(results)



