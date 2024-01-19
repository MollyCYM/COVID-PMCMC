### Simulate data from a stochastic volatility model
rm(list=ls())
set.seed(56798)
library("truncnorm")
### Define true parameters
beta <- 0.6
Phi <- 0.975
cQ <- sqrt(0.02) 

### Simulate SV data
T <- 2000

num=Time=T
x=arima.sim(n=num+1,list(ar=0.975),sd=sqrt(0.02))
y=ts(beta*exp(x[-1]/2)*rnorm(num,0,sqrt(1)))  

### Particle filter for SV model likelihood approximation

# particle_filter_sv <- function(y, phi, sigma, n_particles = 1000) {
#   particles <- matrix(0, nrow = n_particles, ncol = length(y))
#   weights <- rep(1 / n_particles, n_particles)
#   
#   for (t in 2:length(y)) {
#     particles[, t] <- phi * particles[, t - 1] + rnorm(n_particles, 0, sigma)
#     weights <- weights * dnorm(y[t] - 0, mean = 0, sd = true_sigma_obs)  # Assuming normal observation model
#     weights <- weights / sum(weights)
#   }
#   
#   return(list(particles = particles, weights = weights))
# }

# SIRPF<-function(y,N,beta,Phi,cQ){
#   particlesSIR0<-matrix(0,N,1)
#   weightsSIR0<-matrix(0,N,1)
#   particlesSIR <- matrix(0,N,T)
#   resampleSIR<-matrix(0,N,T)
#   weightsSIR <- matrix(0,N,T)
#   WeightsSIR <- matrix(0,N,T)
#   totalweight<-vector(length=T)
#   likelihood<-vector(length=T)
#   hSIRPF <- vector(length=T) 
#   hSIRvar<- vector(length=T)
#   for (i in 1:N) {
#     particlesSIR0[i,1]=rnorm(1, 0 , sd=sqrt((cQ^2)/(1 - Phi^2)))
#     weightsSIR0[i,1]=1/N
#   }
#   for (t in 1:T){
#     for (i in 1:N){
#       if(t==1){particlesSIR[i,t] = Phi*(particlesSIR0[i,1]) + rnorm(1,0,cQ) 
#       weightsSIR[i,t]= dnorm(y[t],mean=0,sd=sqrt(exp(particlesSIR[i,t])*beta^2))
#       }
#       if(t>1){particlesSIR[i,t]=Phi*(resampleSIR[i,(t-1)]) + rnorm(1,0,cQ) 
#       weightsSIR[i,t]=dnorm(y[t],mean=0,sd=sqrt(exp(particlesSIR[i,t])*beta^2))
#       
#       }}
#     
#     likelihood[t]=mean(weightsSIR[,t])
#     totalweight = sum(weightsSIR[,t])
#     for(i in 1:N){ 
#       WeightsSIR[i,t]=weightsSIR[i,t]/totalweight[t]}
#     
#     #Resampling step:
#     resampleSIR[,t] = sample(particlesSIR[,t], size=N, replace=TRUE, prob =WeightsSIR[,t])
#     hSIRPF[t] = mean(resampleSIR[,t])
#     hSIRvar[t]= var(resampleSIR[,t])
#   }
#   loglikelihood=sum(log(likelihood))
#   return(loglikelihood)}
SIRPF <- function(y, N, beta, Phi, cQ) {
  particlesSIR0 <- matrix(0, N, 1)
  weightsSIR0 <- matrix(0, N, 1)
  particlesSIR <- matrix(0, N, length(y))
  resampleSIR <- matrix(0, N, length(y))
  weightsSIR <- matrix(0, N, length(y))
  WeightsSIR <- matrix(0, N, length(y))
  totalweight <- vector(length = length(y))
  likelihood <- vector(length = length(y))
  hSIRPF <- vector(length = length(y))
  hSIRvar <- vector(length = length(y))
  
  for (i in 1:N) {
    particlesSIR0[i, 1] = rnorm(1, 0, sd = sqrt((cQ^2) / (1 - Phi^2)))
    weightsSIR0[i, 1] = 1 / N
  }
  
  for (t in 1:length(y)) {
    for (i in 1:N) {
      if (t == 1) {
        particlesSIR[i, t] = Phi * (particlesSIR0[i, 1]) + rnorm(1, 0, cQ)
        weightsSIR[i, t] = dnorm(y[t], mean = 0, sd = sqrt(exp(particlesSIR[i, t]) * beta^2))
      }
      if (t > 1) {
        particlesSIR[i, t] = Phi * (resampleSIR[i, (t - 1)]) + rnorm(1, 0, cQ)
        weightsSIR[i, t] = dnorm(y[t], mean = 0, sd = sqrt(exp(particlesSIR[i, t]) * beta^2))
      }
    }
    
    likelihood[t] = mean(weightsSIR[, t])
    totalweight[t] = sum(weightsSIR[, t])
    
    if (is.finite(totalweight[t]) && totalweight[t] > 0) {
      WeightsSIR[, t] = weightsSIR[, t] / totalweight[t]
    } else {
      # Handle the case where totalweight is zero or NA
      WeightsSIR[, t] = 1 / N
    }
    
    # Resampling step:
    resampleSIR[, t] = sample(particlesSIR[, t], size = N, replace = TRUE, prob = WeightsSIR[, t])
    hSIRPF[t] = mean(resampleSIR[, t])
    hSIRvar[t] = var(resampleSIR[, t])
  }
  
  loglikelihood = sum(log(likelihood))
  return(loglikelihood)
}
set.seed(93)
log_likelihood <- SIRPF(y, N = 200, beta = 0.6, Phi = 0.975, cQ = sqrt(0.02))
print(log_likelihood)

# PMMH algorithm for SV model
pmmh_sv <- function(data, n_iterations) {
  # Initialize parameters
  beta <- 0.8
  phi <- 0.45
  cQ <- sqrt(0.04)
  
  # MCMC samples
  beta_samples <- numeric(n_iterations)
  phi_samples <- numeric(n_iterations)
  cQ_samples <- numeric(n_iterations)
  
  for (iteration in 1:n_iterations) {
    # Sample new parameters using a proposal distribution
    #? not vector form? j-1...
    beta_proposal <- rtruncnorm(1, a=0, b=Inf, mean=beta, sd=0.1)
    phi_proposal <- rtruncnorm(1, a=-1, b=1, mean=phi, sd=0.01)
    cQ_proposal <- rnorm(1, mean=cQ, sd=0.001)
    
    # Compute the likelihood using particle filter
    # particle_result <- particle_filter_sv(data, phi_proposal, sigma_proposal)
    
    # Estimate the marginal likelihood
    marginal_likelihood <- SIRPF(y,N=200,beta=beta_proposal,Phi=phi_proposal,cQ=cQ_proposal)
    
    # Calculate the acceptance probability
    acceptance_prob <- min(1, marginal_likelihood)
    
    # Accept or reject the proposal
    if (runif(1) < acceptance_prob) {
      beta <- beta_proposal
      phi <- phi_proposal
      cQ <- cQ_proposal
    }
    
    # Save samples
    beta_samples[iteration] <- beta
    phi_samples[iteration] <- phi
    cQ_samples[iteration] <- cQ
  }
  
  return(list(beta_samples=beta_samples, phi_samples = phi_samples, cQ_samples = cQ_samples))
}

# Run the PMMH algorithm for SV model
result <- pmmh_sv(y,n_iterations = 10)




