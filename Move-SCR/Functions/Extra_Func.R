# ==============================================================================
# POPULATION SIZE ESTIMATION UNDER THE Move-SCR MODEL
# By Clara Panchaud
# ==============================================================================

# A helper function to start: 

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



#' Form the transition matrix Q
#'
#' @param theta the two parameters included in Q
#' @param ac_dist a matrix of the distances between the activity centre and the grid cells,
#'                such that ac_dist[i,k] is the distance between the activity centre i
#'                and the centre of grid cell k
#'                
#' @return  three lists of length equal to the number of simulated individuals, where Q_list contains
#'          the Q matrices that are now individual specific, same for pi_list containing the initial state probabilities
#'          and s0 the initial state locations. 


make_Q<-function(theta,ac_dist,mesh,N,S){
  
  Q_list <- list()
  pi_list <- list()
  s0 <- numeric(N)
  
  neighbour<-neighbour_matrix(mesh)
  
  for (i in 1:N){
    Q<-neighbour
    for (j in 1:S){
      for (k in 1:S){
        if(Q[j,k]==1){
          Q[j,k]= exp( theta[1] - theta[2] * ac_dist[i,k]) 
        }
      }
    }
    diag(Q)<- - rowSums(Q)
    
    eig <- eigen(t(Q))
    pi <- Re(eig$vectors[, which.min(abs(eig$values))])
    pi <- abs(pi) / sum(abs(pi))  # Ensure positive and normalized
    
    # Safety check
    if(any(pi < 0) || any(!is.finite(pi))) {
      pi <- rep(1/S, S)  # Fallback to uniform distribution
    }
    
    initial_state <- sample(1:S, size = 1, prob = pi)
    
    Q_list[[i]] <- Q
    pi_list[[i]] <- pi
    s0[i] <- initial_state
  }
  
  return(list(Q_list,pi_list,s0))
}



#' Estimate population size and detection probability under the UMove-SCR model
#'
#' Computes the estimated population size (N̂) and the overall detection 
#' probability (1 - P_unobserved) based on the model parameters and the survey setup.
#'
#' @param theta numeric(2). Model parameters on the log-scale: 
#'        c(q, log(lambda)).
#'         where exp(q)     : baseline transition rate between states
#                lambda     : detection rate in states with traps on
#' @param Time numeric(1). Total duration of the survey.
#' @param observed_ind integer(1). Number of individuals observed at least once.
#' @param S integer(1). Number of spatial states in the mesh.
#' @param neighbour S x S matrix. Adjacency matrix defining
#'        neighbouring relationships between spatial states.
#' @param traps_on integer vector. Indices of states where traps are active.
#' @param f numeric vector of length S. Initial probability distribution over states.
#'
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item obs_prob: overall probability of being observed at least once
#'         }
#'         
#'         
Pop_size_mmmpp <- function(theta, Time, observed_ind, S, neighbour, traps_on, f, camera_count = NULL) {
  q_hat <- exp(theta[1]) # transition rate
  l_hat <- exp(theta[2]) # detection rate
  
  # Construct generator matrix for movement
  Q_est <- neighbour * q_hat
  diag(Q_est) <- -rowSums(Q_est)
  
  # If camera_count not provided, create vector of ones
  if (is.null(camera_count)) {
    camera_count <- rep(1, S)
  }
  
  # Construct detection rate matrix (lambda_i = 0 outside trap states)
  lambda_est <- rep(0, S)
  lambda_est[traps_on] <- l_hat * camera_count[traps_on]
  Lambda_mat <- diag(lambda_est)
  
  # Compute the probability of remaining unobserved over the entire survey
  # via matrix exponential of (Q - Lambda)
  unobs_prob <- sum(f %*% expm(Time * (Q_est - Lambda_mat)))
  
  # Detection probability
  obs_prob <- 1 - unobs_prob
  
  # Population size estimated via Horvitz–Thompson-like estimator
  N_est <- observed_ind / obs_prob
  
  return(c(N_est, obs_prob))
}



#' Compute population size estimate and its standard error
#'
#' @param theta_est numeric(2). Estimated parameters.
#' @param Time numeric(1). Total survey duration.
#' @param observed_ind integer(1). Number of individuals observed in the data.
#' @param S integer(1). Number of spatial states.
#' @param neighbour S x S matrix. Adjacency matrix.
#' @param traps_on integer vector. Indices of states with traps.
#' @param f numeric vector of length S. Initial state distribution.
#' @param cov 2x2 covariance matrix of parameter estimates.
#' @param distribution character(1). "poisson" or "binomial", determines the 
#'        sampling variance component.
#'     
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item SE: standard error of N_est
#'         }
confint_pop_mmmpp <- function(theta_est, Time, observed_ind, S, neighbour, traps_on, f, cov, camera_count = NULL, distribution = "binomial") {
  
  # Compute point estimates of population size and detection probability using the previous function
  Pop_size <- Pop_size_mmmpp(theta_est, Time, observed_ind, S, neighbour, traps_on, f, camera_count)
  N_est <- Pop_size[1]
  obs_prob <- Pop_size[2]
  
  # Numerical gradient of N_est with respect to theta (for delta method)
  grad_vec <- numDeriv::grad(
    func = function(theta) Pop_size_mmmpp(theta, Time, observed_ind, S, neighbour, traps_on, f, camera_count),
    x = theta_est
  )
  
  # Sampling variance component: depends on assumed distribution
  s2 <- switch (tolower(distribution),
                poisson  = observed_ind/obs_prob^2,
                binomial = observed_ind*((1-obs_prob)/(obs_prob^2)))
  
  # Total variance
  Var <- as.numeric(t(grad_vec) %*% cov %*% grad_vec) + s2
  SE <- sqrt(Var)
  
  return(c(N_est=N_est, SE=SE))
}



# ==============================================================================
# Repeat with the ACMove-SCR model 
# ==============================================================================



#' Estimate population size and detection probability under the ACMove-SCR model
#'
#' Computes the estimated population size (N̂) and the overall detection 
#' probability (1 - P_unobserved) based on the model parameters and the survey setup.
#'
#' @param theta numeric(3). Model parameters on the log-scale: 
#'        c(log(lambda), alpha, beta).
#' @param Time numeric(1). Total duration of the survey.
#' @param observed_ind integer(1). Number of individuals observed at least once.
#' @param S integer(1). Number of spatial states in the mesh.
#' @param neighbour S x S matrix. Adjacency matrix defining
#'        neighbouring relationships between spatial states.
#' @param traps_on integer vector. Indices of states where traps are active.
#' @param camera_count integer vector. Number of cameras present in each state cell. (We allow max 1 in our example, but this vector allows this assumption to be relaxed.)
#' @param n_ac numeric(1). Number of cells in the activity centre mesh used for integration. 
#' @param ac_to_states_dist N_ac x S matrix.  Matrix defining the distances between each potential activity centre locations and the state space. 
#' 
#' 
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item obs_prob: overall probability of being observed at least once
#'         }

Pop_size_full_ac <- function(theta, Time, observed_ind, S, neighbour, traps_on, camera_count, n_ac, ac_to_states_dist) {
  l <- exp(theta[1]) # Detection rate lambda
  theta1 <- theta[2] # alpha
  theta2 <- theta[3] # beta
  
  # Construct detection rate matrix with the camera count multiplying lambda
  lambda_est <- rep(0, S)
  lambda_est[traps_on] <- l * camera_count[traps_on]
  Lambda_mat <- diag(lambda_est)
  
  unobs_prob <- 0
  
  # Construct generator matrix for movement
  for (i in 1:n_ac) {
    Q <- neighbour
    
    for (j in 1:S) {
      for (k in 1:S) {
        if (Q[j,k] == 1) {
          Q[j,k] = exp(theta1 - theta2 * ac_to_states_dist[i, k])
        }
      }
    }
    diag(Q) <- -rowSums(Q)
    
    # Use the stationary distribution as the initial distribution
    eig <- eigen(t(Q))
    pi <- Re(eig$vectors[, which.min(abs(eig$values))])
    pi <- pi / sum(pi)
    
    # Compute the probability of remaining unobserved over the entire survey
    # via matrix exponential of (Q - Lambda)
    unobs_prob <- unobs_prob + sum(pi %*% expm(Time * (Q - Lambda_mat)))
  }
  
  p_unobs = unobs_prob / n_ac
  p_obs   <- max(1 - p_unobs, 1e-12)  
  
  # Population size estimated via Horvitz–Thompson-like estimator
  N_est <- observed_ind / p_obs
  
  return(c(N_est, p_obs))  
}


#' This function returns ONLY the population estimate from the previous function, as a scalar
Pop_size_mmmpp_ac <- function(theta, Time, observed_ind, S, neighbour, traps_on, camera_count,  n_ac, ac_to_states_dist) {
  Pop_full <- Pop_size_full_ac(theta, Time, observed_ind, S, neighbour, traps_on, camera_count,  n_ac, ac_to_states_dist)
  return(Pop_full[1])  
}


#' Compute population size estimate and its standard error
#'
#' @param theta_est numeric(3). Estimated parameters.
#' @param Time numeric(1). Total survey duration.
#' @param observed_ind integer(1). Number of individuals observed in the data.
#' @param S integer(1). Number of spatial states.
#' @param neighbour S x S matrix. Adjacency matrix.
#' @param traps_on integer vector. Indices of states with traps.
#' @param camera_count integer vector. Number of cameras present in each state cell. 
#' @param cov 2x2 covariance matrix of parameter estimates.
#' @param n_ac numeric(1). Number of cells in the activity centre mesh used for integration. 
#' @param ac_to_states_dist N_ac x S matrix.  Matrix defining the distances between each potential activity centre locations and the state space. 
#' @param distribution character(1). "poisson" or "binomial", determines the 
#'        sampling variance component.
#'     
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item SE: standard error of N_est
#'         }
# Fixed confidence interval function
confint_pop_mmmpp_ac <- function(theta_est, Time, observed_ind, S, neighbour, traps_on, camera_count, cov, n_ac, ac_to_states_dist, distribution = "binomial") {
  
  # Compute point estimates of population size and detection probability using the previous function
  Pop_full <- Pop_size_full_ac(theta_est, Time, observed_ind, S, neighbour, traps_on, camera_count,  n_ac, ac_to_states_dist)
  N_est <- Pop_full[1]
  obs_prob <- Pop_full[2]
  
  # Calculate gradient using the scalar function
  grad_vec <- numDeriv::grad(
    func = function(theta) Pop_size_mmmpp_ac(theta, Time, observed_ind, S, neighbour, traps_on, camera_count,  n_ac, ac_to_states_dist),
    x = theta_est
  )
  
  # Sampling variance component: depends on assumed distribution
  s2 <- switch(tolower(distribution),
               poisson  = observed_ind/obs_prob^2,
               binomial = observed_ind*((1-obs_prob)/(obs_prob^2)))
  
  # Total variance
  Var <- as.numeric(t(grad_vec) %*% cov %*% grad_vec) + s2
  SE <- sqrt(Var)
  
  return(c(N_est=N_est, SE=SE))
}




#' Convert MMPP simulation output to discretised detection data used for CT SCR
#'
#' Turns the MMMPP data into a discretized detection history format suitable 
#' for analysis with CT SCR. The function maps 
#' detections to trap locations, removes duplicate detections within a time 
#' threshold, and discretises continuous time data into intervals.
#'
#' @param mmpp List. MMPP simulation output containing:
#'        \itemize{
#'          \item ss: spatial location index of detections
#'          \item tt: detection times
#'        }
#' @param traps_on Integer vector. Indices of active trap locations.
#' @param T Numeric(1). Total survey duration.
#' @param r2 Numeric(1). Parameter passed to discretise function.
#' @param mesh Data frame. Mesh defining the state space of the MMMPP, here defining the trap locations with columns x and y.
#'
#' @return List of length 2:
#'         \itemize{
#'           \item ddfmat_sim: matrix of discretised detection data
#'           \item dfrows_sim: numeric vector with number of detections per individual
#' 
mmpp_to_df<-function(mmpp,traps_on,T,r2,mesh){
  # Convert MMPP simulation output to observation data frame
  # Map over spatial locations (ss) and times (tt) to create tibbles
  obs_df <- map2(
    .x = mmpp$ss,
    .y = mmpp$tt,
    .f = ~ tibble(y = .x, Time = .y)
  ) %>%
    # Combine all observations into single data frame with ID column for each individual
    bind_rows(.id = "id") %>%
    mutate(id = as.integer(id))
  
  # Convert to base R data frame
  obs_df<-as.data.frame(obs_df)
  
  # Remove last row as it is an artifact 
  obs_df<-obs_df[-nrow(obs_df),]
  
  # Add trap coordinates from mesh based on trap locations
  obs_df$trap_x <- mesh[obs_df$y,]$x
  obs_df$trap_y <- mesh[obs_df$y,]$y
  
  # Remap trap IDs to match the active traps
  obs_df$y<- match(obs_df$y, traps_on)
  
  # Remove duplicate detections that occur within 1 hour (0.04166667 days)
  # This filters out rapid re-detections at the same trap
  df_sim <- obs_df %>%
    arrange(id, Time) %>% 
    filter(!(id == lag(id) & (Time - lag(Time) < 0.04166667)))
  
  # Pre-allocate data frame for discretized observations
  n_rows <- sum(table(df_sim$id))
  ddf_sim <- data.frame(
    t = numeric(n_rows),
    y = integer(n_rows),
    id = integer(n_rows)
  )
  
  row_idx <- 1
  # Loop through each individual and discretise their detection history
  for(i in unique(df_sim$id)){
    # Discretise observations for individual i into time intervals using r2
    data <- discretize(df_sim[df_sim$id == i, ],T, r2)
    data$id <- i
    
    # Determine how many rows this individual contributes
    n_new <- nrow(data)
    
    # Insert discretised data into pre-allocated data frame
    ddf_sim[row_idx:(row_idx + n_new - 1), ] <- data
    
    # Update row index for next individual
    row_idx <- row_idx + n_new
    
    # Convert to matrix format
    ddfmat_sim = as.matrix(ddf_sim)
    
    # Count detections per individual
    dfrows_sim = as.numeric(table(ddf_sim$id))
  }
  
  return(list(ddfmat_sim, dfrows_sim))
}


