# ==============================================================================
# DISCRETISATIONS TESTS FOR SURVEY 1
# By Clara Panchaud
# ==============================================================================
# This code tests different time and space discretisations until finding
# the stable values used in Analysis1.R
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
traps_og<-read.csv("Study1/marten_traps1.csv")



## create the mask over the region; different buffers were tested 
## until the estimated parameters remained stable


 space_list <- c(0.1, 0.2,0.5,1)
 L_list <- c(10,50,100,200)
# 
 results<-data.frame(matrix(ncol = 8, nrow = 16))
# 
 colnames(results)<-c( "spacing", "L", "time.m", "time.nm", "N", "N.nm", "Likelihood.m", "Likelihood.nm")
 j<-1
# 
 for (space in space_list){
   for (L in L_list){

     
        trap = make.poly(x=traps_og$x, y=traps_og$y)
        trap <- trap[-31,]
        mask = make.mask(trap,buffer=2,spacing=space,type="trapbuffer")
        meshmat<-as.matrix(mask)
        traps<-as.matrix(trap)
        
        
        T<-12
        n<-length(unique(data$id))
        
        ## make a version of the data set with the discretized time
        ddf <- data.frame(t = as.numeric(), y = as.integer(), id = as.integer())
        for(i in unique(data$id)){
          df <- discretize(data[data$id == i, ],T, L)
          df$id <- i
          ddf <- rbind(ddf, df)
        }
        ddfmat = as.matrix(ddf)
        dfrows = as.numeric(table(ddf$id))
        
        ## inital parameters for MSCR (h0,sigma,beta)
        theta<-c(1.5, -2 ,1.5) 
        
        
        ## fit the MSCR model and save the running time 
        start_time <- Sys.time()
        fit <- optim(theta, LikelihoodC, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T,hessian=TRUE) 
        end_time <- Sys.time()
        fit_time<-end_time - start_time
        
        ## extract the estimated parameters and confidence intervals
        theta_est<-fit$par
        param<-confint_param(fit,T,traps,mask)
        Likelihood.m <- LikelihoodC(theta_est,traps,ddfmat, dfrows, meshmat, T)
        N_est<-confint_pop(fit,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
        
        ## initial parameters for SCR (h0,sigma)
        theta<-c(1.5,-2) 
        
        ## fit the SCR model and save the running time
        start_time <- Sys.time()
        fit_nomem <- optim(theta, LikelihoodCnoMem, trap = traps, df = ddfmat, dfrows = dfrows, mesh = meshmat, endt = T,hessian=TRUE) 
        fit_time_nomem<-end_time <- Sys.time()
        time_nomem<- end_time - start_time
        ## extract the estimated parameters and confidence intervals
        param_nm<-confint_param(fit_nomem,T,traps,mask)
        N_nm<-confint_pop(fit_nomem,T,traps,mask,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05)
        Likelihood.nm <- LikelihoodCnoMem(fit_nomem$par,traps, ddfmat, dfrows, meshmat, T)
        

        fit_time<-as.numeric(fit_time, units="mins")
        time_nomem<-as.numeric(time_nomem, units="mins")
  
        results[j,]<-c(space, L, fit_time, time_nomem, N_est[1], N_nm[1], Likelihood.m, Likelihood.nm)
        j <- j+1

}}
     

results1<-read.csv("results1.csv",row.names = 1)
results1<-na.omit(results1)
results3<-read.csv("results3.csv",row.names = 1)
results3<-na.omit(results3)
results2<-read.csv("results2.csv",row.names = 1)
results2<-na.omit(results2)
results<-rbind(results1,results3,results2)


results
