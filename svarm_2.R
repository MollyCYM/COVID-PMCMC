rm(list = ls())
# Load necessary packages
library(astsa)         # For time series data simulation
library(pomp)          # Particle filtering
library(mvtnorm)       # Multivariate normal distribution
library(truncnorm)
# Model Parameters
n_iter <- 1000000      # Number of MCMC iterations
n_particles <- 100   # Number of particles for the bootstrap particle filter
# Load your data and clean it
y_data <- read.csv('svmgen_y1.csv', header = FALSE, skip = 1)  # Skip first row
y <- as.numeric(y_data[, 2]) 

# Proposal distribution for PMMH (Random walk Metropolis-Hastings) with truncated normals
propose_params <- function(theta) {
  beta_proposal <- rtruncnorm(1, a = 0, b = Inf, mean = theta[1], sd = 0.1)      # Truncated at 0 for beta
  phi_proposal <- rtruncnorm(1, a = -1, b = 1, mean = theta[2], sd = 0.1)        # Truncated between -1 and 1 for phi
  sigma_eta_proposal <- rtruncnorm(1, a = 0, b = Inf, mean = theta[3], sd = 0.1) # Truncated at 0 for sigma_eta
  
  return(c(beta_proposal, phi_proposal, sigma_eta_proposal))
}

# Prior distributions
prior_density <- function(theta) {
  dtruncnorm(theta[1], a = 0, b = Inf, mean = 0.8, sd = 0.1) +  # Truncated normal for beta (positive values)
    dtruncnorm(theta[2], a = -1, b = 1, mean = 0.9, sd = 0.1) +  # Truncated normal for phi (between -1 and 1)
    dtruncnorm(theta[3],  a = 0, b = Inf, mean = sqrt(0.05), sd = 0.1)  # Truncated normal for sigma_eta (positive values)
}

# Bootstrap particle filter for the SV model
bootstrap_filter <- function(y, theta, n_particles) {
  T <- length(y)
  beta <- theta[1]
  phi <- theta[2]
  sigma_eta <- theta[3]
  
  alpha_particles <- matrix(0, n_particles, T)
  weights <- matrix(0, n_particles, T)
  unnormalized_weights <- matrix(0, n_particles, T)  # To store unnormalized weights
  
  # Initial state distribution
  alpha_particles[, 1] <- rnorm(n_particles, mean = 0, sd = sigma_eta / sqrt(1 - phi^2))
  
  for (t in 1:T) {
    if (t > 1) {
      alpha_particles[, t] <- phi * alpha_particles[, t-1] + rnorm(n_particles, 0, sigma_eta)
    }
    
    # Calculate the unnormalized weights (without dividing by sum)
    unnormalized_weights[, t] <- dnorm(y[t], mean = 0, sd = beta * exp(alpha_particles[, t] / 2))
    
    # Normalize weights for resampling (but don't use this for log-likelihood)
    weights[, t] <- unnormalized_weights[, t] / sum(unnormalized_weights[, t]) # Normalizing for resampling only
    
    # Resampling step
    indices <- sample(1:n_particles, size = n_particles, prob = weights[, t], replace = TRUE)
    alpha_particles[, t] <- alpha_particles[indices, t]
  }
  
  # Compute the log-likelihood using unnormalized weights
  log_likelihood <- sum(log(colMeans(unnormalized_weights)))
  
  return(list(alpha_particles = alpha_particles, log_likelihood = log_likelihood))
}

# PMMH Algorithm
pmmh <- function(y, n_iter, n_particles) {
  theta_chain <- matrix(0, n_iter, 3)
  theta_chain[1, ] <- c(0.8, 0.9, 0.1)  # Initial values for beta, phi, sigma_eta
  
  log_likelihood_chain <- numeric(n_iter)
  
  res <- bootstrap_filter(y, theta_chain[1, ], n_particles)
  log_likelihood_chain[1] <- res$log_likelihood
  
  for (iter in 2:n_iter) {
    theta_proposal <- propose_params(theta_chain[iter - 1, ])
    
    # # Check prior constraints
    # if (theta_proposal[1] <= 0 || theta_proposal[3] <= 0 || abs(theta_proposal[2]) > 1) {
    #   theta_chain[iter, ] <- theta_chain[iter - 1, ]
    #   log_likelihood_chain[iter] <- log_likelihood_chain[iter - 1]
    #   next
    # }
    
    # Bootstrap particle filter for the proposed parameters
    res_proposal <- bootstrap_filter(y, theta_proposal, n_particles)
    log_likelihood_proposal <- res_proposal$log_likelihood
    
    # Calculate log-prior ratio
    log_prior_current <- log(dtruncnorm(theta_chain[iter - 1, 1], a = 0, b = Inf, mean = 0.8, sd = 0.1)) +
      log(dtruncnorm(theta_chain[iter - 1, 2], a = -1, b = 1, mean = 0.9, sd = 0.1)) +
      log(dtruncnorm(theta_chain[iter - 1, 3], a = 0, b = Inf, mean = sqrt(0.05), sd = 0.1))
    
    log_prior_proposed <- log(dtruncnorm(theta_proposal[1], a = 0, b = Inf, mean = 0.8, sd = 0.1)) +
      log(dtruncnorm(theta_proposal[2], a = -1, b = 1, mean = 0.9, sd = 0.1)) +
      log(dtruncnorm(theta_proposal[3], a = 0, b = Inf, mean = sqrt(0.05), sd = 0.1))
    
    log_prior_ratio <- log_prior_proposed - log_prior_current
    
    # Compute log likelihood ratio
    log_likelihood_ratio <- log_likelihood_proposal - log_likelihood_chain[iter - 1]
    
    # Calculate the log-proposal density for theta_proposed given theta_current
    log_proposal_current_to_proposed <- log(dtruncnorm(theta_proposal[1], a = 0, b = Inf, mean = theta_chain[iter - 1, 1], sd = 0.1)) +
      log(dtruncnorm(theta_proposal[2], a = -1, b = 1, mean = theta_chain[iter - 1, 2], sd = 0.1)) +
      log(dtruncnorm(theta_proposal[3], a = 0, b = Inf, mean = theta_chain[iter - 1, 3], sd = 0.1))
    # Calculate the log-proposal density for theta_current given theta_proposed (reverse move)
    log_proposal_proposed_to_current <- log(dtruncnorm(theta_chain[iter - 1, 1], a = 0, b = Inf, mean = theta_proposal[1], sd = 0.1)) +
      log(dtruncnorm(theta_chain[iter - 1, 2], a = -1, b = 1, mean = theta_proposal[2], sd = 0.1)) +
      log(dtruncnorm(theta_chain[iter - 1, 3], a = 0, b = Inf, mean = theta_proposal[3], sd = 0.1))
    log_proposal_ratio <- log_proposal_proposed_to_current - log_proposal_current_to_proposed
    
    # Calculate the final log acceptance ratio
    log_accept_ratio <- log_likelihood_ratio + log_prior_ratio + log_proposal_ratio
    
    # Accept or reject
    if (log(runif(1)) < log_accept_ratio) {
      theta_chain[iter, ] <- theta_proposal
      log_likelihood_chain[iter] <- log_likelihood_proposal
    } else {
      theta_chain[iter, ] <- theta_chain[iter - 1, ]
      log_likelihood_chain[iter] <- log_likelihood_chain[iter - 1]
    }
  }
  
  
  return(list(theta_chain = theta_chain, log_likelihood_chain = log_likelihood_chain, final_alpha = res_proposal$alpha_particles))
}
result <- pmmh(y, n_iter, n_particles)
write.csv(result$final_alpha, "../data/svm_alpha2.csv")
write.csv(result$log_likelihood_chain, "../data/svm_likelihood2.csv")
write.csv(result$theta_chain[, 1],"../data/svm_beta2.csv")
write.csv(result$theta_chain[, 2],"../data/svm_phi2.csv")
write.csv(result$theta_chain[, 3],"../data/svm_sigma2.csv")