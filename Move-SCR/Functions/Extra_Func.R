# ==============================================================================
# POPULATION SIZE ESTIMATION UNDER THE Move-SCR MODEL
# By Clara Panchaud
# ==============================================================================


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
Pop_size_mmmpp <- function(theta, Time, observed_ind, S, neighbour, traps_on, f) {
  q_hat <- exp(theta[1]) # transition rate
  l_hat <- exp(theta[2]) # detection rate
  
  # Construct generator matrix for movement
  Q_est <- neighbour * q_hat
  diag(Q_est) <- -rowSums(Q_est)
  
  # Construct detection rate matrix (lambda_i = 0 outside trap states)
  lambda_est <- rep(0, S)
  lambda_est[traps_on] <- l_hat 
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
confint_pop_mmmpp <- function(theta_est, Time, observed_ind, S, neighbour, traps_on, f, cov, distribution = "binomial") {
  
  # Compute point estimates of population size and detection probability using the previous function
  Pop_size <- Pop_size_mmmpp(theta_est, Time, observed_ind, S, neighbour, traps_on, f)
  N_est <- Pop_size[1]
  obs_prob <- Pop_size[2]
  
  # Numerical gradient of N_est with respect to theta (for delta method)
  grad_vec <- numDeriv::grad(
    func = function(theta) Pop_size_mmmpp(theta, Time, observed_ind, S, neighbour, traps_on, f),
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

