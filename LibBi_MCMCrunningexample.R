#LibBi deterministic model MCMC inference part
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)

# Load the data
v <- read.csv("365days.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time", "value")
ncores <- 8 
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
  obs y

  state S
  state E
  state I
  state R
  
  input N
  param sigma
  param beta
  param gamma
  param tau
  param mu

  sub parameter {
    sigma ~ truncated_gaussian(0.20379467, 0.2, lower = 0) 
    gamma ~ truncated_gaussian(0.12460946, 0.2, lower = 0) 
    beta ~ truncated_gaussian(0.57586873, 0.3, lower = 0) 
    mu ~ truncated_gaussian(0.09454979, 0.5, lower = 0) 
    tau ~ uniform(0, 1)
  }

  sub initial {
    S <-999999
    E <- 1 
    I <-0
    R <-0
  }

  sub transition(delta = 1) {
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dS/dt = -(beta*S*I)/N
      dE/dt = (beta*S*I)/N - sigma*E
      dI/dt = sigma*E - gamma*I - mu*I
      dR/dt = (mu+gamma)*I
    }
  }

  sub observation {
    y ~ log_normal(log(max((sigma*E)/5, 0)), tau)
  }

  sub proposal_parameter {
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    beta ~ gaussian(beta, 0.01)
    mu ~ gaussian(mu,0.001)
    tau ~ gaussian(tau, 0.05)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 1000000)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(sigma =0.19, gamma =0.2, beta=0.53, mu =0.0015)

bi <- sample(bi_model, end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'model') %>% 
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 500, thin = 1)


bi_lst <- bi_read(bi %>% sample_obs)