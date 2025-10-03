# ==============================================================================
# MESH GENERATION AND LIKELIHOOD FUNCTIONS
# By Clara Panchaud
# ==============================================================================
# These functions create dual-resolution meshes and calculate the likelihood for
# the MSCR model
# ==============================================================================

# BACKGROUND: Data discretization approach
# ------------------------------------------------------------------------------
# The dataset is assumed to come in continuous time (as in simulation), i.e. for 
# each individual we have a list of capture times and traps where the animal was seen. 
# The "discretize" function creates a discretisation of time from 0 to T in r segments.
# Segments in which the individual is unseen get the value 0, while segments where 
# the individual is observed get cut in two at the event time. The segment until the 
# event time gets the value of the trap where it is seen, and so on. 
# This discretization is done for one individual at a time.
# 
# NOTE: The Likelihood function takes "mesh" as a vector of meshes (it also works 
# for a single mesh), which enables the blue/red mesh setting


# ==============================================================================
# FUNCTION: new_mesh
# ==============================================================================
# Creates a dual-resolution mesh with fine resolution near observed locations
# and coarse resolution elsewhere
# ------------------------------------------------------------------------------
#' Generate dual-resolution spatial mesh
#' 
#' @param data data frame. Capture history of one individual in discretised format
#' @param cams Matrix with 2 columns. Each row contains a camera trap's coordinates
#' @param space numeric(1). Spacing between mesh points in coarse mesh
#' @param mask mesh object (see secr package) with 2 columns, representing potential 
#'        activity centre locations for integration (created with the spacing entered above)
#' @param threshold numeric(1). Distance threshold defining "near" observed locations
#' @param fine_divisor integer(1). Factor by which to divide spacing in fine mesh
#' 
#' @return list with two elements:
#'         coarse - mesh object for areas far from observations
#'         fine - mesh object for areas near observations
#' 
new_mesh <- function(data, cams, space, mask, threshold, fine_divisor = 2) {
  
  # Identify the cameras where the individual was detected
  cams_ind<-cams[unique(data$y)[!unique(data$y)==0], ,drop=FALSE]
  
  d1 <- nrow(cams_ind) # number of cameras with detections
  d2 <- nrow(mask) # number of mask points at the start
  
  # Calculate distances from all mask points to observed cameras
  dists <- as.matrix(dist(rbind(as.matrix(cams_ind), as.matrix(mask))))
  a <- dists[(d1 + 1):(d1 + d2), 1:d1]  # distance matrix: mask points × cameras
  
  # Identify mask points within threshold distance of observed cameras
  if (d1==1){
    close_index <-unique(which(a<threshold))
  }else {
    close_index<-unique(which(a < threshold, arr.ind = TRUE)[,1])
  }
  close<-mask[close_index,]  # mask points requiring fine resolution
  
  # Create fine mesh around each close point
  fine_mask_list <- vector("list", nrow(close))
  
  for (j in seq_len(nrow(close))) {
    
    # Define coarse and fine spacing parameters
    coarse_spacing <- space
    fine_spacing <- space/ fine_divisor
    center_x <- close[j,1]
    center_y <- close[j,2]
    
    x_seq <- seq(center_x - coarse_spacing/2 + fine_spacing/2,
                 center_x + coarse_spacing/2 - fine_spacing/2,
                 by = fine_spacing)
    y_seq <- seq(center_y - coarse_spacing/2 + fine_spacing/2,
                 center_y + coarse_spacing/2 - fine_spacing/2,
                 by = fine_spacing)
    
    fine_points <- expand.grid(x = x_seq, y = y_seq)
    fine_mask_list[[j]] <- fine_points
  }
  
  # Combine all fine mesh points and remove duplicates
  fine_mask <- unique(do.call(rbind, fine_mask_list))
  
  # Create coarse mesh by removing points replaced by fine mesh
  coarse_mask <- mask[-close_index,]
  dummy_mask<-coarse_mask[1:nrow(fine_mask),]
  dummy_mask$x <- fine_mask$x
  dummy_mask$y <- fine_mask$y
  
  # Set cell areas for integration
  attr(coarse_mask,"a")<-space^2 * 10^-4
  attr(dummy_mask,"a")<-fine_spacing^2 *10^-4
  
  # Return list with coarse and fine meshes
  list(coarse = coarse_mask, fine = dummy_mask)
}



# ==============================================================================
# FUNCTION: Likelihood_integrate
# ==============================================================================
# Integrate likelihood over activity center locations using dual-resolution mesh
# ------------------------------------------------------------------------------
#' Calculate integrated likelihood across spatial mesh(es)
#' 
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param data data frame obtained from discretize() for this individual
#' @param mesh list. One or more mesh objects for integration
#' @param memory integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' 
#' @return numeric. Log-likelihood value for the individual's capture history, with activity centre integrated out
Likelihood_integrate_mesh<-function(theta,trap,data,mesh,memory){
  
  Likeli<-0 # integrated likelihood accumulator
  A<-0  # total area accumulator
  
  # Loop through all meshes (coarse and fine)
  for (mask in mesh){
    # Skip empty meshes ()
    if (dim(mask)[1]!=0){
      a = attr(mask,"a" ) # grid cell area
      D<-dim(mask)[1]  # number of grid cells in the mesh
      A = A+a*D # total mesh area
      L<-0 # Initalised likelihood for this mesh
      mask<-matrix(c(mask[,1],mask[,2]),ncol=2,nrow=D)
      
      # Sum likelihood across all mesh points
      for (i in 1:D){
        L<-L+Likelihood_ind(theta,trap,data,mask[i,],memory)
      }
      
      # Add weighted contribution to total likelihood
      Likeli<-Likeli+(L*a)
    }
  }
  
  # Return log-likelihood normalised by the total area
  return(log(Likeli/A))
}







