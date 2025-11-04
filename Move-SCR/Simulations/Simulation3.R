################################################################################
# STUDY 3: Cross-Model Performance
# Author: Clara Panchaud
# 
# Description:
# This script simulates data from the ACMove-SCR model and fits three models: 
# ACMove-SCR, UMove-SCR and the standard CT SCT model. 
################################################################################

# Load required libraries ------------------------------------------------------
library(ggplot2)
library(secr)
library(TMB)
library(expm)
library(fields)
library(tidyverse)


# Source custom functions ------------------------------------------------------
source("Move-SCR/Functions/sim_MMPP.r")
source("Move-SCR/Functions/Extra_Func.r")

# Functions from MSCR to fit CT SCR
source("MSCR/Functions/Fit_Func.r")
source("MSCR/Functions/Sim_Func.r")
Rcpp::sourceCpp("MSCR/Functions/LikelihoodC.cpp")


# Compile and load TMB model for UMove and ACMove ------------------------------
compile("Move-SCR/Functions/like_MMPP.cpp")
dyn.load(dynlib("Move-SCR/Functions/like_MMPP"))

compile("Move-SCR/Functions/like_MMPP_ac.cpp")
dyn.load(dynlib("Move-SCR/Functions/like_MMPP_ac"))

# Set random seed for reproducibility
set.seed(1)

# Load and prepare trap data ---------------------------------------------------

# Read trap locations
traps <- read.csv("Data/Study2/marten_traps2.csv")

# Create trap polygon
trap <- make.poly(x = traps$x, y = traps$y)

# Make the mask/state space and some of its metrics
mask = make.mask(trap,buffer=2,spacing=1,type="trapbuffer")
mesh<-mask

S <-dim(mask)[1] #State space size
a_states <- attr(mask, "a") # area of one cell
A_total <- a_states * S #T otal area covered

# Print state space size
print(paste("Simulation grid states (S):", S))

# Convert traps to matrix format
traps<-as.matrix(traps)

# Map traps to the state space
distance_traps<-proxy::dist(traps,mask)
traps_loc<-apply(distance_traps,1,which.min)

# Handle multiple cameras in same cell 
camera_count <- rep(1, S)
trap_counts <- table(traps_loc)
camera_count[as.numeric(names(trap_counts))] <- as.numeric(trap_counts)


# Set simulation parameters ----------------------------------------------------------------------------------------

# Study duration
Time <- 264  # hours (11 days)

# True population size
N <- 30

# Movement parameter (log scale)
q <- log(0.05)

# Detection rate parameter
l <- 0.02

# Create detection matrix Lambda for simulation --------------------------------------------------------------------

lambda <- rep(0,S)
lambda[traps_loc]<-l *camera_count[traps_loc]

# set up the neighbour matrix for later building Q
neighbour<-neighbour_matrix(mask)

# Define the meshes for the activity centres  --------------------------------------------------------------

# Fine mesh to simulate the acs
mesh_ac = make.mask(trap,buffer=2,spacing=0.01,type="trapbuffer")

# Coarser mesh for the fitting
mesh_ac_fit = make.mask(trap,buffer=2,spacing=0.8,type="trapbuffer")
a_ac <- attr(mesh_ac_fit, "a") # area of one cell    # Note: Unsure of this should be mesh_ac or mesh_ac_fit.
n_ac<-dim(mesh_ac_fit)[1] # number of grid points

# Set up the distances between state space and activity centre mesh, divide by constant to scale -------------------

Dist0 <- 0.25
Dist_fit <- fields::rdist(mesh_ac_fit, mesh)
ac_fit_dist <- Dist_fit / Dist0 

# Define simulation settings ---------------------------------

# Different beta values we are testing
beta_grid <- c(0 , 0.05,  0.2) 

# Number of simulations per beta
n_sim<-2

# Initialise results data frame
col_names <- c(
  "alpha_true", "beta_true" ,"N_true", "Time", "S", "n_ac",
  "n_obs", "n_events",
  "lambda_hat", "alpha_hat", "beta_hat",
  "se_lambda", "se_alpha", "se_beta", "N_hat", "SE_N",
  "q", "lambda_UMove", "q_SE", "lambda_Umove_SE", "N_Umove", "SE_N_Umove",
  "h0", "sigma", "se_h0", "se_sigma", "N_SCR", "SE_N_SCR", 
  "nll", "AIC", "elapsed_sec", "elapsed_UMove", "elapsed_SCR"
)
results <- as.data.frame(matrix(NA, nrow = n_sim*9, ncol = length(col_names)))
colnames(results) <- col_names

# Run simulations --------------------------------------------------------------

i_sim <- 1



for (beta in beta_grid){
    
    theta<- c(q , beta)
    
    for (i in 1:n_sim){
      
      
      ac<-mesh_ac[sample(1:nrow(mesh_ac), N, replace = TRUE), ]
      D_ac <- fields::rdist(ac, mesh)
      ac_dist <- D_ac / Dist0 
      
      pre_data<-make_Q(theta, ac_dist,mesh,N,S)
      Q_list<-pre_data[[1]]
      pi_list<-pre_data[[2]]
      s0<-pre_data[[3]]
      
 
      obj <- sim_CTMC_ac(Q_list,Time,s0,S,N)
      
      
      mmpp <- sim_MMPP(obj,lambda)
      
      
      observed_idx <- which(!vapply(mmpp$tt, is.null, logical(1)))
      
      # Subset tt and ss
      mmpp$tt <- mmpp$tt[observed_idx]
      mmpp$ss <- mmpp$ss[observed_idx]
      pi_list  <- pi_list[observed_idx]
      ac_dist  <- ac_dist[observed_idx[-length(observed_idx)], , drop = FALSE]
      
      observed_ind <- length(mmpp$tt)-1
      lambda_fixed <- rep(1, S)
      lambda_fixed[traps_loc] <- 0
      lambda_template <- rep(0, S)
      lambda_template[traps_loc] <- 1
      
      
      Q_fixed <- 1-neighbour
      Q_template <- 2*neighbour  #see if Q and Lambda is the first parameter?
      
      # build the TMB data
      tmbdata <- list()
      
      tmbdata$Us <- as.integer(unlist(lapply(mmpp$tt, length)))
      tmbdata$t <- unlist(mmpp$tt)
      tmbdata$s <- as.integer(unlist(mmpp$ss)-1)
      pi_matrix <- do.call(rbind, pi_list)
      tmbdata$pi_matrix <- pi_matrix
      
      tmbdata$n_ac <- n_ac
      
      # where we have an NA in tt (no history) set the length to zero
      tmbdata$Us[unlist(lapply(mmpp$tt, \(x) all(is.na(x))))] <- 0
      
      tmbdata$t0 <- 0
      tmbdata$Time <- Time
      tmbdata$n_states <- S
      tmbdata$n_indiv <- length(mmpp$tt)
      tmbdata$camera_count <- camera_count
      
      tmbdata$n_obs <- observed_ind
      tmbdata$ac_dist <- ac_dist
      tmbdata$area_cell_ac <- a_ac
      tmbdata$area_region <- A_total
      tmbdata$ac_to_states_dist <- ac_fit_dist
      
      
      tmbdata$Q_template <- Q_template
      tmbdata$Q_fixed <- Q_fixed
      tmbdata$lambda_template <- lambda_template
      tmbdata$lambda_fixed <- lambda_fixed
      
      # Inference
      
      # starting values
      theta_init <- c( log (0.01),log(0.04),0.3)
      
      start_time <- Sys.time()
      
      
      # make the autodiff function
      obj <- MakeADFun(data=tmbdata, parameters=list(theta=theta_init), DLL="like_MMPP_ac")
      
      
      #obj$fn(theta_init)
      # do the optimization
      out <- nlminb(theta_init, obj$fn, obj$gr)
      end_time <- Sys.time()
      
      rep<-sdreport(obj)
      
      summ<-summary(rep)
      
      
      Pop<-confint_pop_mmmpp_ac(out$par, Time, observed_ind, S, neighbour, traps_loc,  camera_count,  rep$cov, n_ac,ac_fit_dist)
      nll <- obj$fn(out$par)
      AIC <- 2*3 + nll*2
      
      n_events <- sum(lengths(mmpp$tt))
      elapsed_sec <- as.numeric(difftime(end_time, start_time, units="secs"))
      
      N_hat <- Pop[1]
      SE_N  <- Pop[2]
      

      
      ## Now the Umove-SCR model
      
      
      Q_fixed <- 1 - neighbour         # positions of structural zeros
      Q_template <- neighbour          # positions to fill with transition rates (neighbours)
      
      lambda_fixed <- rep(1,S)
      lambda_fixed[traps_loc] <- 0  # free parameters at active trap locations only
      lambda_template<- rep(0,S)
      lambda_template[traps_loc] <- 2
      
      
      
      tmbdata2 <- list(
        Us = lengths(mmpp$tt),
        t = unlist(mmpp$tt),
        s = as.integer(unlist(mmpp$ss) - 1L),  # 0-based indexing for C++
        n_obs = observed_ind,
        f = rep(1/S, S),
        t0 = 0,
        Time = Time,
        n_states = S,
        n_indiv = length(mmpp$tt),
        Q_template = Q_template,
        Q_fixed = Q_fixed,
        lambda_template = lambda_template,
        lambda_fixed = lambda_fixed,
        camera_count = camera_count)
      
      # where we have an NA in tt (no history) set the length to zero (in our case that will just be the last term)
      tmbdata2$Us[unlist(lapply(mmpp$tt, \(x) all(is.na(x))))] <- 0
      
      
      theta_init <- rep(-2.5,2)
      
      start_time2 <- Sys.time()
      # make the autodiff function
      obj2 <- MakeADFun(data=tmbdata2, parameters=list(theta=theta_init), DLL="like_MMPP")
      
      # do the optimization
      invisible(capture.output({
        fit2 <- nlminb(start = theta_init, objective = obj2$fn, gradient = obj2$gr)
      }))
      
      
      end_time2 <- Sys.time()
      
      report2<- sdreport(obj2, par.fixed = fit2$par)
      param2<-summary(report2)
      
      elapsed_sec2 <- as.numeric(difftime(end_time2, start_time2, units="secs"))
      
      Pop_size2<- confint_pop_mmmpp(fit2$par, Time, observed_ind, S, neighbour, traps_loc, camera_count,  report2$cov)
      
      
      ## now the CT SCR model
      
      r2<-10
      cams<-mask[traps_loc,]
      mesh_fit = make.mask(trap,buffer=2,type="traprect",spacing=0.5)
      theta_init<-c(1,0.01,0.1)
      new_format <- mmpp_to_df(mmpp,traps_loc,Time,r2)
      ddfmat_sim <- new_format[[1]]
      dfrows_sim <- new_format[[2]]
      
      # fit CT SCR
      start_time3 <- Sys.time()
      fit_nomem <- optim(theta_init[1:2], LikelihoodCnoMem, trap = as.matrix(cams),df = ddfmat_sim, dfrows =dfrows_sim, mesh = as.matrix(mesh_fit), endt = Time, hessian=TRUE)
      end_time3 <- Sys.time()
    
      N_est_SCR<-confint_pop(fit_nomem,Time,cams,mesh_fit,observed_ind)
      SE_SCR<-sqrt(diag(solve(fit_nomem$hessian)))
      elapsed_sec3 <- as.numeric(difftime(end_time3, start_time3, units="secs"))
      
      
      results[i_sim, ] <- list(
        q, beta, N, Time, S, n_ac, observed_ind, n_events,
        summ[4,1], summ[5,1], summ[3,1],
        summ[4,2], summ[5,2], summ[3,2], N_hat, SE_N,
        param2[3,1], param2[4,1],  param2[3,2], param2[4,2], Pop_size2[1], Pop_size2[2],
        fit_nomem$par[1], fit_nomem$par[2], SE_SCR[1], SE_SCR[2], as.numeric(N_est_SCR[1]), as.numeric(N_est_SCR[2]),
        nll, AIC, elapsed_sec,elapsed_sec2,elapsed_sec3)
      
      i_sim <- i_sim+1
      
    }
  }


