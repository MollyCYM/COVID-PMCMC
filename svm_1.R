rm(list = ls())
# Load necessary packages
library(mcmcse)        # For MCMC diagnostics
library(astsa)         # For time series data simulation
library(pomp)          # Particle filtering
library(mvtnorm)       # Multivariate normal distribution

# Model Parameters
n_iter <- 1000      # Number of MCMC iterations
n_particles <- 100   # Number of particles for the bootstrap particle filter


# # Load your observation data
# y_data <- read.csv('svm_y2r.csv', header = FALSE)[,1]  # Assuming your observations are in the first column
# y <- y_data  # This replaces the simulated 'y' with your real data
# 
# set.seed(123)
# beta_true <- 0.6         # True value of beta
# phi_true <- 0.975        # True value of phi
# sigma_eta_true <- sqrt(0.02) # True value of sigma_eta
# Simulated data
set.seed(123)
T <- 2000                # Number of time steps
beta_true <- 0.6         # True value of beta
phi_true <- 0.975        # True value of phi
sigma_eta_true <- sqrt(0.02) # True value of sigma_eta

# Simulate latent states (AR(1) process)
alpha <- numeric(T)
alpha[1] <- rnorm(1, mean = 0, sd = sigma_eta_true / sqrt(1 - phi_true^2))
for (t in 2:T) {
  alpha[t] <- phi_true * alpha[t-1] + rnorm(1, 0, sigma_eta_true)
}
y <- numeric(T)
for (t in 1:T) {
# Simulate observations
y <- beta_true * exp(alpha / 2) * rnorm(1,0,1)}

# Log-likelihood of the observations given states
log_likelihood <- function(y, alpha, beta) {
  sum(dnorm(y, mean = beta * exp(alpha / 2), sd = 1, log = TRUE))
}

write.csv(x, "svmgen_alpha1.csv")
write.csv(y, "svmgen_y1.csv")
y_data <- read.csv('svmgen_y1.csv', header = FALSE)  # Assuming your observations are in the first column
# Remove the first two rows
y_data <- y_data[-c(1, 1), ]
# Keep only the second column
y <- y_data[, 2]
# Proposal distribution for PMMH (Random walk Metropolis-Hastings)
propose_params <- function(theta) {
  rnorm(3, mean = theta, sd = c(0.1, 0.01, 0.01))  # Propose new beta, phi, sigma_eta
}

# Prior distributions
prior_density <- function(theta) {
  dnorm(theta[1], 0.8, 0.1, log = TRUE) +       # Prior for beta
    dnorm(theta[2], 0.45, 1, log = TRUE) +      # Prior for phi
    dnorm(theta[3], sqrt(0.04), 0.1, log = TRUE) # Prior for sigma_eta
}

# Bootstrap particle filter for the SV model
bootstrap_filter <- function(y, theta, n_particles) {
  T <- length(y)
  beta <- theta[1]
  phi <- theta[2]
  sigma_eta <- theta[3]
  
  alpha_particles <- matrix(0, n_particles, T)
  weights <- matrix(0, n_particles, T)
  
  # Initial state distribution
  alpha_particles[, 1] <- rnorm(n_particles, mean = 0, sd = sigma_eta / sqrt(1 - phi^2))
  
  for (t in 1:T) {
    if (t > 1) {
      alpha_particles[, t] <- phi * alpha_particles[, t-1] + rnorm(n_particles, 0, sigma_eta)
    }
    weights[, t] <- dnorm(y[t], mean = beta * exp(alpha_particles[, t] / 2), sd = 1)
    weights[, t] <- weights[, t] / sum(weights[, t]) # Normalize weights
    
    # Resampling
    indices <- sample(1:n_particles, size = n_particles, prob = weights[, t], replace = TRUE)
    alpha_particles[, t] <- alpha_particles[indices, t]
  }
  
  log_likelihood <- sum(log(colMeans(weights)))
  
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
    
    # Check prior constraints
    if (theta_proposal[1] <= 0 || theta_proposal[3] <= 0 || abs(theta_proposal[2]) > 1) {
      theta_chain[iter, ] <- theta_chain[iter - 1, ]
      log_likelihood_chain[iter] <- log_likelihood_chain[iter - 1]
      next
    }
    
    # Bootstrap particle filter for the proposed parameters
    res_proposal <- bootstrap_filter(y, theta_proposal, n_particles)
    log_likelihood_proposal <- res_proposal$log_likelihood
    
    # Compute log prior ratio
    log_prior_ratio <- prior_density(theta_proposal) - prior_density(theta_chain[iter - 1, ])
    
    # Compute acceptance ratio
    log_accept_ratio <- (log_likelihood_proposal - log_likelihood_chain[iter - 1]) + log_prior_ratio
    
    # Accept or reject
    if (log(runif(1)) < log_accept_ratio) {
      theta_chain[iter, ] <- theta_proposal
      log_likelihood_chain[iter] <- log_likelihood_proposal
    } else {
      theta_chain[iter, ] <- theta_chain[iter - 1, ]
      log_likelihood_chain[iter] <- log_likelihood_chain[iter - 1]
    }
  }
  
  return(list(theta_chain = theta_chain, log_likelihood_chain = log_likelihood_chain, final_alpha = res$alpha_particles))
}

# Run PMMH algorithm
result <- pmmh(y, n_iter, n_particles)
# Load necessary libraries
library(mvtnorm)

# PMMH Result (Assume the 'pmmh' function is already executed and results are available)
# Use the final parameter estimates from the last iteration
beta_estimate <- result$theta_chain[n_iter, 1]
phi_estimate <- result$theta_chain[n_iter, 2]
sigma_eta_estimate <- result$theta_chain[n_iter, 3]

# Use the last set of latent states (alpha) from the last iteration of the particle filter
# Assuming alpha_particles stores the latent states across particles in the last iteration
alpha_estimate <- apply(result$final_alpha, 2, mean)  # Take mean of particles
# Load alpha generated data
x_data <- read.csv('svm_x2.csv', header = FALSE)  # Assuming your observations are in the first column
# Remove the first two rows
x_data <- x_data[-c(1, 2), ]
# Keep only the second column
alpha_saved <- x_data[, 2]
# Display the cleaned data (optional)
head(alpha_saved) 

# Plot the new generated observations
plot(alpha_saved, type = "l", col = "blue", lwd = 2, ylim = range(c(y, new_observations)),
     main = "Comparison of Input Observations and New Generated Observations", 
     xlab = "Time", ylab = "Observations")
lines(alpha_estimate, col = "red", lwd = 2)
legend("topright", legend = c("Input Observations", "New Generated Observations"),
       col = c("blue", "red"), lty = 1, lwd = 2)


write.csv(result$latent_states, "svm_alpha1.csv")
write.csv(result$theta_chain[, 1],"svm_beta1.csv")
write.csv(result$theta_chain[, 2],"svm_phi1.csv")
write.csv(result$theta_chain[, 3],"svm_sigma1.csv")