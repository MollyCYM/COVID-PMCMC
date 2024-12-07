rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
library(readr)
# Load the data
v <- read.csv("259wk1.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()

y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
  obs y
  
  const h = 1
  state S
  state E
  state I
  state R
  state x

  state Z

  input N
  param k
  param gamma
  param sigma // Noise driver
  param theta
  param a
  param b
  param E0
  param I0
  param R0
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(1.59, 0.02, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(1.08, 0.075, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    theta ~ truncated_gaussian(0.05, 0.01, lower = 0)
    a ~ gaussian(-0.02, 0.01)
    b ~ gaussian(-0.2, 0.1)
    x0 ~ uniform(0,2)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
    tau ~ uniform(0, 1)
  }

  sub initial {
    S <- N
    R <- R0*S
    S <- S - R

    E <- exp(E0 + log(S))
    S <- S - E
    I <- exp(I0 + log(S))
    S <- S - I
    x<- x0
    Z <- 0
  }

  sub transition(delta = h) {
    Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = theta*(a+b-x)+sigma*e/h
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5.0, 0)), tau)
  }

  sub proposal_parameter {
    k ~ gaussian(k, 0.005)
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    theta ~ truncated_gaussian(theta, 0.01,lower=0)
    a ~ gaussian(a, 0.01)
    b ~ gaussian(b, 0.01)
    x0 ~ gaussian(x0, 0.05)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
    tau ~ gaussian(tau, 0.05)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 100000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)