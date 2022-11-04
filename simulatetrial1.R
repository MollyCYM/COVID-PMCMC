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
options(digits=2)
# Load the data
v <- read.csv("simulate366.csv", header=FALSE, stringsAsFactors=FALSE) 
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
  state M
  
  input N
  param sigma
  param beta
  param gamma
  param tau
  param mu

  sub parameter {
    sigma ~ truncated_gaussian(0.25, 0.2, lower = 0) 
    gamma ~ truncated_gaussian(0.2, 0.2, lower = 0) // gamma is the period, not the rate
    beta ~ truncated_gaussian(0.5, 0.3, lower = 0) 
    mu ~ truncated_gaussian(0.001, 0.25, lower = 0) 
    tau ~ uniform(0, 1)
  }

  sub initial {
    S <-999999
    E <- 1 
    I <-0
    R <-0
    M <-0
  }

  sub transition(delta = 1) {
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dS/dt = -(beta*S*I)/N
      dE/dt = (beta*S*I)/N - sigma*E
      dI/dt = sigma*E - gamma*I - mu*I
      dR/dt = gamma*I
      dM/dt = mu*I
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

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 100000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)
write.csv(bi_lst,"../data/SEIR1.csv")
fitY <- bi_lst$y %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(y %>% rename(Y = value))
write.csv(fitY,"../data/SEIRy1.csv")


write.csv(bi_lst$sigma$value,"../data/SEIRsigma1.csv")
write.csv(bi_lst$gamma$value,"../data/SEIRgamma1.csv")
write.csv(bi_lst$beta$value,"../data/SEIRbeta1.csv")
write.csv(bi_lst$mu$value,"../data/SEIRmu1.csv")


