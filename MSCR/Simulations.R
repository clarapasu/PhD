# ==============================================================================
# MSCR SIMULATION STUDY
# By Clara Panchaud
# ==============================================================================
# This simulation compares three spatial capture-recapture models:
# 1. MSCR (Memory Spatial Capture-Recapture) - incorporates memory effects
# 2. CT SCR (Continuous Time SCR) - no memory, continuous time
# 3. DT SCR (Discrete Time SCR) - traditional discrete time model
#
# The study evaluates how different values of the memory parameter beta 
# affect population estimation accuracy across these models.
# ==============================================================================

# Load required libraries
# ------------------------------------------------------------------------------
require(secr)
library(ggplot2)
library(numDeriv)
library(dplyr)
library(expm)
library(fields)
library(ggplot2)
library(proxy)
library(fdrtool)
library(ggpubr)
library(MASS)
library(tidyr)

# Source custom functions for likelihood calculations and simulation
# ------------------------------------------------------------------------------
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")

# set seed for reproducible results
set.seed(1)


# ==============================================================================
# STUDY SETUP AND PARAMETERS
# ==============================================================================

# Load the trap locations (from a real survey) and setup the spatial configuration
traps<-read.csv("Study2/Data/marten_traps2.csv",col.names = c("x","y"))
trap = make.poly(x=traps$x, y=traps$y)
trap <- trap[-31,] # Remove the last element as it was repeated 
cams<- read.traps(data = traps,detector="count") # we need different format of the trap locations 
                                                 # so they are named trap and cams 

T_days<-11 # Length of survey period (days)
N<-20 # True total number of individuals 
r<-1000 # Number of bins in the time discretization used to generate data
m<-1 # Memory indicator (1 = memory effects present, 0 = no memory (=SCR) )
h0<- 0.85 # Baseline hazard parameter (log-scale)
sigma<- -1.16 # Spatial scale parameter (log-scale)
beta<- c(-1.5,-1,-0.5,0,0.5,1,1.5) # Memory parameter values to test (log-scale)
theta_init<- c(1.00,  0.01 , 0.00) # Initial parameter values for optimisation 


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
n_sim <- 100

# Initialize results dataframe
# ------------------------------------------------------------------------------
results <- data.frame(matrix(ncol = 25, nrow = n_sim * length(beta)))
colnames(results) <- c(
  # MSCR model results (with memory)
  "h0", "sigma", "beta",                    # Parameter estimates
  "SE.h0", "SE.sigma", "SE.beta",          # Standard errors
  "N", "SE.N",                             # Population estimates
  
  # CT SCR model results (no memory)
  "h0.nomem", "sigma.nomem",               # Parameter estimates
  "SE.h0.nomem", "SE.sigma.nomem",         # Standard errors
  "N.nomem", "SE.N.nomem",                 # Population estimates
  
  # DT SCR model results (traditional)
  "lambda0", "sigma.secr",                      # Parameter estimates (secr package)
  "SE.lambda0", "SE.sigma.secr",                # Standard errors
  "N.secr", "SE.secr",                     # Population estimates
  
  # Additional information
  "n",                                     # Number of observed individuals
  "true.beta",                             # True beta value used
  "time.secr", "time.memory", "time.no.memory"  # Computation times
)

j<-1 # Results row counter 

# ==============================================================================
# MAIN SIMULATION LOOP
# ==============================================================================

# Loop over different beta values (memory parameter)
for (k in 1:length(beta)){
  
  # Set current parameter vector
  theta<-c(h0,sigma,beta[k])  
  
  # Run multiple simulations for each value of beta
  for (i in 1:n_sim){
    
    # ----------------------------------------------------------------------
    # DATA SIMULATION
    # ----------------------------------------------------------------------
    
    # Sample activity centers for the N individuals
    sampling<-sample(nrow(samplemask),N,replace=TRUE)
    ac<-as.matrix(samplemask[sampling,])
    
    # simulate a data set of capture histories
    df_sim<-sim_data(N,T_days,r,theta,as.matrix(cams),ac,m,samplemask)
    df_sim<-re_id(df_sim) # re-assign individual IDs since unobserved individuals were removed 
    
    # Add trap coordinates to the capture histories
    df_sim$trap_x<-trap[df_sim$y,1]
    df_sim$trap_y<-trap[df_sim$y,2]
    df_sim<-arrange(df_sim,id,Time)
    
    # ----------------------------------------------------------------------
    # DATA PREPARATION FOR MODEL FITTING
    # ----------------------------------------------------------------------
    
    # Define coarser time discretisation for model fitting 
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
    N_est<-confint_pop(fit,T_days,trap,mask,n)

    # Store MSCR results
    results[j,1:3] <- fit$par                         # Parameter estimates
    results[j,4:6] <- sqrt(diag(solve(fit$hessian)))  # Standard errors
    results[j,7:8] <- N_est[1:2]                      # Population estimate and SE
    
    
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
    fit_nomem <- optim(theta[1:2], LikelihoodCnoMem, trap = as.matrix(cams),df = ddfmat_sim, dfrows =dfrows_sim, mesh = meshmat, endt = T_days, hessian=TRUE)
    end_time <- Sys.time()
    time_nomemory<- end_time - start_time
    
    # Calculate population estimate
    N_est_nomem<-confint_pop(fit_nomem,T_days,trap,mask,n)
    
    # Store CT SCR results
    results[j, 9:10] <- fit_nomem$par                         # Parameter estimates
    results[j, 11:12] <- sqrt(diag(solve(fit_nomem$hessian))) # Standard errors
    results[j, 13:14] <- N_est_nomem[1:2]                     # Population estimate and SE
    
    # ----------------------------------------------------------------------
    # MODEL 3: DT SCR (USING SECR PACKAGE)
    # ----------------------------------------------------------------------
    
    # Prepare data in secr package format
    capt <- data.frame(session = rep("test",each=dim(df_sim)[1]),
                       ID = df_sim$id,
                       occasion = ceiling(df_sim$Time), # Convert to discrete occasions
                       trapID = as.character(df_sim$y),
                       stringsAsFactors = FALSE)
    
    # Create capture history object
    sim_ch<-make.capthist(capt,cams)
    
    
    start_time <- Sys.time()
    # Fit traditional discrete-time SCR model
    afit = secr.fit(sim_ch,mask=mask,detectfn = 14,trace = FALSE)
    end_time <- Sys.time()
    time_scr<- end_time - start_time
    
    # Calculate population estimate from density estimate
    N_secr <- exp(coef(afit)[1,1]) * A # Density * Area
    SE_N_secr<- A* sqrt(diag(solve(afit$fit$hessian)))[1]* exp(coef(afit)[1,1])

    # Store DT SCR results
    results[j, 15:16] <- afit$fit$estimate[2:3]                    # Parameter estimates
    results[j, 17:18] <- sqrt(diag(solve(afit$fit$hessian)))[2:3]  # Standard errors
    results[j, 19] <- N_secr                                       # Population estimate
    results[j, 20] <- SE_N_secr                                    # Standard error
    
    # ----------------------------------------------------------------------
    # STORE ADDITIONAL INFORMATION
    # ----------------------------------------------------------------------
    
    results[j,21]<-n # Number of observed individuals
    results[j,22]<-beta[k] # True beta value
    results[j,23:25]<-c(time_scr,time_memory,time_nomemory) # Computation times
    
    j<-j+1 # Move to next results row
  }
  
}

# ==============================================================================
# RESULTS PROCESSING AND VISUALIZATION
# ==============================================================================

# Set true population size
true_N <- 20


# Calculate performance metrics grouped by beta value
# ------------------------------------------------------------------------------
metrics_by_beta <- results %>%
  group_by(true.beta) %>%
  summarise(

    # ==================================================
    # MSCR MODEL METRICS
    # ==================================================
    
    # Calculate CI bounds from estimates and standard errors
    CI_lower = N - 1.96 * SE.N,
    CI_upper = N + 1.96 * SE.N,
    
    # Percentage Bias: ((estimated - true) / true) * 100
    bias_mscr = mean((N - true_N) / true_N * 100, na.rm = TRUE),
    
    # 95% CI Width: mean width of confidence intervals
    ci_width_mscr = mean(CI_upper - CI_lower, na.rm = TRUE),
    
    # Coverage Probability: % of CIs containing true value
    coverage_mscr = mean(CI_lower <= true_N & CI_upper >= true_N, na.rm = TRUE) * 100,
    
    # Root Mean Square Error
    rmse_mscr = sqrt(mean((N - true_N)^2, na.rm = TRUE)),
    
    # ==================================================
    # CT SCR MODEL METRICS 
    # ==================================================
    
    CI_lower.nomem = N.nomem - 1.96 * SE.N.nomem,
    CI_upper.nomem = N.nomem + 1.96 * SE.N.nomem,
    
    bias_nomem = mean((N.nomem - true_N) / true_N * 100, na.rm = TRUE),
    ci_width_nomem = mean(CI_upper.nomem - CI_lower.nomem, na.rm = TRUE),
    coverage_nomem = mean(CI_lower.nomem <= true_N & CI_upper.nomem >= true_N, na.rm = TRUE) * 100,
    rmse_nomem = sqrt(mean((N.nomem - true_N)^2, na.rm = TRUE)),
    
    # ==================================================
    # DT SCR MODEL METRICS 
    # ==================================================
    
    CI_lower.secr = N.secr - 1.96 * SE.secr,
    CI_upper.secr = N.secr + 1.96 * SE.secr,
    
    bias_secr = mean((N.secr - true_N) / true_N * 100, na.rm = TRUE),
    ci_width_secr = mean(CI_upper.secr - CI_lower.secr, na.rm = TRUE),
    coverage_secr = mean(CI_lower.secr <= true_N & CI_upper.secr >= true_N, na.rm = TRUE) * 100,
    rmse_secr = sqrt(mean((N.secr - true_N)^2, na.rm = TRUE)),
    
    .groups = 'drop'
  )


# Process the format of the results for plotting
plot_data <- results %>%
  dplyr::select(true.beta, N.secr, N, N.nomem) %>%
  pivot_longer(
    cols = c(N.secr, N, N.nomem),
    names_to = "model",
    values_to = "N"
  ) %>%
  mutate(
    model = case_when(
      model == "N.secr" ~ "DT SCR",
      model == "N" ~ "MSCR",
      model == "N.nomem" ~ "CT SCR"
    ),
    model = factor(model),
    true.beta = as.factor(true.beta)
  )

# Plot the results of the three models for the different values of beta (Figure 3.4 in thesis)
ggplot(plot_data, aes(x = factor(true.beta), y = N, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = 20, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("white", "gray80", "gray50")) +
  labs(title = "Population Estimates by Model",
       x = expression(log(beta)), y = "Population Estimate (N)") +
  theme_minimal() +
  theme(legend.position = "bottom")

