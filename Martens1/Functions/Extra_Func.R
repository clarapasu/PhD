

Pop_size_mmmpp <- function(theta, Time, observed_ind, S, neighbour, traps_on, f) {
  q_hat <- exp(theta[1])
  l_hat <- exp(theta[2])
  
  Q_est <- neighbour * q_hat
  diag(Q_est) <- -rowSums(Q_est)
  
  lambda_est <- rep(0, S)
  lambda_est[traps_on] <- l_hat 
  Lambda_mat <- diag(lambda_est)
  
  unobs_prob <- sum(f %*% expm(Time * (Q_est - Lambda_mat)))
  N_est <- observed_ind / (1 - unobs_prob)
  
  obs_prob <- 1 - unobs_prob
  return(c(N_est, obs_prob))
}


confint_pop_mmmpp <- function(theta_est, Time, observed_ind, S, neighbour, traps_on, f, cov, distribution = "binomial") {
  
  
  Pop_size <- Pop_size_mmmpp(theta_est, Time, observed_ind, S, neighbour, traps_on, f)
  
  N_est <- Pop_size[1]
  obs_prob <- Pop_size[2]
  
  
  grad_vec <- numDeriv::grad(
    func = function(theta) Pop_size_mmmpp(theta, Time, observed_ind, S, neighbour, traps_on, f),
    x = theta_est
  )
  
  s2 <- switch (tolower(distribution),
                poisson  = observed_ind/obs_prob^2,
                binomial = observed_ind*((1-obs_prob)/(obs_prob^2)))
  
  
  Var <- as.numeric(t(grad_vec) %*% cov %*% grad_vec) + s2
  SE <- sqrt(Var)
  
  return(c(N_est=N_est, SE=SE))
}

