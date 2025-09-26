# ==============================================================================
# ACTIVITY CENTRE PROBABILITY DENSITY FUNCTION PLOTS
# by Clara Panchaud
# ==============================================================================
# This script creates spatial probability density plots showing the likely
# locations of an individual's activity centre based on its capture history.
# Compares the plots obtained from using
# MSCR (with memory) vs CT SCR (without memory).
# 
# AC PDF = Activity Center Probability Density Function
# Shows where an individual is most likely to have its activity centre
# based on the spatial pattern of its captures
# ==============================================================================

# Load required libraries
# ------------------------------------------------------------------------------
library(lubridate)
library(dplyr)
library(ggplot2)
library(secr)
library(sf)

# Load custom functions
# ------------------------------------------------------------------------------
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")

# ==============================================================================
# DATA LOADING
# ==============================================================================

# Load trap locations and capture data
traps<-read.csv("Study1/Data/marten_traps1.csv")
df<-read.csv("Study1/Data/marten_data1.csv")


# Load pre-fitted model results from previous analysis
# These .Rdata files contain the fitted model objects from "Analysis.R"
fit <- get(load("Study1/Models/fit_MSCR.Rdata")) 
fit_nomem <- get(load("Study1/Models/fit_CTSCR.Rdata"))

# Extract parameter estimates from fitted models
theta_est<-fit$par  # MSCR parameters: (h0, sigma, beta)
theta_est_nm<-fit_nomem$par # CT SCR parameters: [h0, sigma]

# ==============================================================================
# SPATIAL MESH SETUP
# ==============================================================================

# Create spatial mesh for AC PDF calculation
# Space parameter controls mesh resolution (larger = coarser, faster computation)
space <- 0.05  # 0.05 km spacing (use 0.05 for thesis plots, but increase
               # to get results faster)

# Convert traps to required formats
trap_mat<-as.matrix(traps)
trap_poly = make.poly(x=traps$x, y=traps$y)
trap_poly <- trap_poly[-31,]

# Create spatial mask around trap array
# buffer=2: extends 2km beyond trap locations
mask = make.mask(trap_poly,buffer=2,spacing=space,type="trapbuffer")
meshmat<-as.matrix(mask)

# ==============================================================================
# INDIVIDUAL DATA PREPARATION
# ==============================================================================

# Select an individual for AC PDF demonstration
# (This could be changed to any individual ID in the dataset, 
# for the thesis we run this code with a,b,c)
## select an individual of the dataset and the traps where it was seen
data<-df[df$id==8,]

# Create dataframe tracking which traps caught this individual
trap_ind<-data.frame( x = trap_poly$x,
                      y = trap_poly$y,
                      seen = rep(0,30),
                      seen_indic=rep(0,30))

# Fill in capture information for each trap
for ( i in 1:length(table(data$y))){
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen<-table(data$y)[i]
  trap_ind[as.numeric(names(table(data$y))[i]),]$seen_indic<-1
}

# Discretise the continuous capture times using 12 days divided by 100 (as in Analysis1.R)
data<-discretize(data,12,100)

# ==============================================================================
# AC PDF CALCULATION - MSCR MODEL (WITH MEMORY)
# ==============================================================================


# Calculate the likelihood of activity center being at each mesh point
# For MSCR model (memory = 1)
ac_posterior_ind<-apply(meshmat,1,Likelihood_ind,theta=theta_est,trap=trap_poly,data=data,memory=1)

# Calculate normalisation constant (integral over entire mesh)
# This ensures the AC PDF integrates to 1 (is a proper probability density)
integral<-as.double(exp(Likelihood_integrate(theta_est,trap_poly,data,mask,1)))

# Create normalised AC PDF data frame
ac_density_ind<-data.frame(x = meshmat[,1],
                           y = meshmat[,2],
                           value = unlist(ac_posterior_ind)/integral)

# ==============================================================================
# AC PDF PLOT - MSCR MODEL
# ==============================================================================

# Create AC PDF plot for MSCR model for the selected individual
cell <- attr(mask, "spacing")
plot_ind <- ggplot(ac_density_ind, aes(x = x, y = y, fill = value)) +
  geom_tile(width = cell, height = cell) +
  geom_point(
    data = trap_ind,
    aes(x = x, y = y, size = as.factor(seen)),  
    shape = 23,
    fill = ifelse(trap_ind$seen == 0, "white", "grey"),
    color = "black",
    inherit.aes = FALSE
  ) +
  coord_fixed(xlim = c(314.5, 317.5), ylim = c(4965, 4968), expand = FALSE) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 200))+
  scale_x_continuous(breaks = c(315, 316, 317))+
  labs(
    x = "Easting (km)",
    y = "Northing (km)",
    fill = "AC PDF",
    size = "Number of captures"
  ) +
  scale_size_discrete(
    range = c(5, 10),
    breaks = sort(unique(trap_ind$seen)), # only those that occur
    guide = guide_legend(
      override.aes = list(
        fill = ifelse(sort(unique(trap_ind$seen)) == 0, "white", "grey"),
        shape = 23,
        color = "black"
      )
    )
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.key = element_rect(fill = "white")
  )


#================================================================================
# AC PDF CALCULATION - CT SCR MODEL (NO MEMORY)
# ==============================================================================

# Repeat AC PDF calculation for CT SCR model (memory = 0)
ac_posterior_ind_nm<-apply(meshmat,1,Likelihood_ind,theta=theta_est_nm,trap=trap_poly,data=data,memory=0)
integral<-as.double(exp(Likelihood_integrate(theta_est_nm,trap_poly,data,mask,0)))
ac_density_ind_nm<-data.frame(x = meshmat[,1],
                              y = meshmat[,2],
                              value = unlist(ac_posterior_ind_nm)/integral)

# ==============================================================================
# AC PDF PLOT - CT SCR MODEL  
# ==============================================================================

# Create AC PDF plot for MSCR model for the selected individual
cell <- attr(mask, "spacing")
plot_ind_nm <- ggplot(ac_density_ind_nm, aes(x = x, y = y, fill = value)) +
  geom_tile(width = cell, height = cell) +
  geom_point(
    data = trap_ind,
    aes(x = x, y = y, size = as.factor(seen)),  
    shape = 23,
    fill = ifelse(trap_ind$seen == 0, "white", "grey"),
    color = "black",
    inherit.aes = FALSE
  ) +
  coord_fixed(xlim = c(314.5, 317.5), ylim = c(4965, 4968), expand = FALSE) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 200))+
  scale_x_continuous(breaks = c(315, 316, 317))+
  labs(
    x = "Easting (km)",
    y = "Northing (km)",
    fill = "AC PDF",
    size = "Number of captures"
  ) +
  scale_size_discrete(
    range = c(5, 10),
    breaks = sort(unique(trap_ind$seen)), # only those that occur
    guide = guide_legend(
      override.aes = list(
        fill = ifelse(sort(unique(trap_ind$seen)) == 0, "white", "grey"),
        shape = 23,
        color = "black"
      )
    )
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.key = element_rect(fill = "white")
  )

# ==============================================================================
# DISPLAY AND SAVE THE PLOTS 
# ==============================================================================

library(ggplot2)
library(patchwork)

# Remove legends from both plots
plot_ind_noleg <- plot_ind + theme(legend.position = "none")
plot_ind_nm_noleg <- plot_ind_nm + theme(legend.position = "none")

# Use one of the original plots to supply the legend
# and collect it below the combined plots
final_plot <- plot_ind + plot_ind_nm_noleg + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")  # horizontal legend

# Display
final_plot

# Save high-quality figure
#ggsave("Figures/AC_PDF_horizontal_legend_patchwork.png", final_plot, width = 12, height = 6, dpi = 300)
