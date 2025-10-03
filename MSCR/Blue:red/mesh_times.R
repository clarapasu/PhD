# ==============================================================================
# BLUE/RED MESH SIMULATIONS
# By Clara Panchaud
# ==============================================================================
# This code explores the computational efficiency of the MSCR model by testing a method
# combining spatial meshes, such that the mesh is finer in areas where the individual was 
# observed (and therefore where the likelihood contribution will be greater)
# ==============================================================================

# Load required libraries
# ------------------------------------------------------------------------------
library(fdrtool)
library(secr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(numDeriv)

# Source custom functions for likelihood calculations and simulation
# ------------------------------------------------------------------------------
Rcpp::sourceCpp("Functions/LikelihoodC.cpp")
source("Functions/Sim_Func.R")
source("Functions/Fit_Func.R")
source("Functions/Fit_mesh.R")

# set seed for reproducible results
set.seed(1)

# Load the data and trap locations (from a real survey)
df<-read.csv("Study2/Data/marten_data2.csv")
traps_og<-read.csv("Study2/Data/marten_traps2.csv",col.names = c("x","y"))

# ==============================================================================
# FIRST TEST AND PLOT
# Test the code to produce the coarser and finer meshes
# Plot the result
# ==============================================================================

space<- 0.5 # length (km) of one side of a coarse cell
threshold<-2 # threshold distance from observed traps that the fine mesh covers
divisor <- 2 # divisor^2 is the number of fine cells covering one coarse cell

# Create trap polygon and mask for survey area, convert them to matrices for computation
trap = make.poly(x=traps_og$x, y=traps_og$y)
trap <- trap[-31,]
mask = make.mask(trap,buffer=2,spacing=space,type="trapbuffer")
meshmat<-as.matrix(mask)
traps<-as.matrix(trap)


# Select individual for analysis and prepare the data
T_days<-12 # Survey length
r<-10 # Temporal discretisation parameter
i<-2 # Index of chosen individual to run the tests on
data<-discretize(df[df$id==i,],T_days,r) # Obtain discretised version of the individual's capture history

# Identify cameras where the individual was detected
cams_ind<-traps[unique(data$y)[!unique(data$y)==0],,drop=FALSE]

# Generate the dual-resolution mesh (blue/red meshes)
mesh<-new_mesh(data,traps,space,mask,threshold, divisor)

# ==============================================================================
# VISUALIZATION: Plot study area with both meshes and trap locations
# vary space, threshold and divisor to obtain all the plots from Figure 3.6 of the thesis. 
# ==============================================================================

# Get extent of the mask for the plot boundaries
x_range <- range(mask$x) 
y_range <- range(mask$y) 
y_range[1]<-y_range[1]-1
y_range[2]<-y_range[2]+1

#Plot the landscape and meshes
plot(NA, xlim =x_range, ylim = y_range, asp=1, 
     xlab = "Longitude", ylab = "Latitude", main = "Study Area")
box(col = NA) 
plot(mask, dots = FALSE, border = 1, ppoly = FALSE, add = TRUE)
plotMaskEdge(mask, add = TRUE, col = "black", lwd = 1.5)
points(mesh[[1]], col = "red", pch = 16, cex = 0.5) # red mesh (coarse resolution)
points(mesh[[2]], col = "blue", pch = 16, cex = 0.5) # blue mesh (fine resolution)
points(traps, col="black",pch=4, cex=1.5) # all trap locations
points(cams_ind, col="yellow",pch=8, cex=2) # traps with detections


# ==============================================================================
# TEST OF MESH CONFIGURATIONS
# ==============================================================================

# Define model parameters (we select the results from Analysis2.R)
h0<- 0.85
sigma<- -1.17
beta<- 0.74
theta<-c(h0,sigma,beta)

# Define parameter grids for the analysis
space_list <- c(1, 2)
threshold_list <- c(1.5, 2)
divisor_list <- c(2,4, 6, 8)

# Initialise the results data frame
results<-data.frame(matrix(ncol = 7, nrow = 16))
colnames(results)<-c( "dim 1", "dim 2", "Likelihood", "time", "space", "threshold", "divisor")

# Loop through all parameter combinations
j<-1
for (space in space_list){
  for (threshold in threshold_list){
    for (divisor in divisor_list){

      # Create mask with current spacing parameter
      mask = make.mask(trap,buffer=2,spacing=space,type="trapbuffer")
      
      # Generate dual-resolution mesh with current parameters
      mesh<-new_mesh(data,traps,space,mask,threshold,divisor)

      # Time the likelihood calculation
      start_time <- Sys.time()
      L <-Likelihood_integrate_mesh(theta,trap,data,mesh,1)
      end_time <- Sys.time()
      fit_time<- end_time - start_time
      fit_time <- as.numeric(fit_time, units="mins")

      # Record mesh dimensions
      dim1 <- dim(mesh[[1]])[1]
      dim2 <- dim(mesh[[2]])[1]

      results[j,]<- c(dim1, dim2, L, fit_time, space, threshold, divisor)
      j<-j+1
    }
  }
}

# Calculate the combined mesh dimension for each row of the results
results$dim <- results$`dim 1` + results$`dim 2`


#Display the results (Some are extracted to form Table 3.8 of the thesis)
print(results, row.names = FALSE)



