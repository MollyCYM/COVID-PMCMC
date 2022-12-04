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
v <- read.csv("5days.csv", header=FALSE, stringsAsFactors=FALSE) 
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
    sigma ~ truncated_gaussian(0.20379467, 0.2, lower = 0) 
    gamma ~ truncated_gaussian(0.12460946, 0.2, lower = 0) // gamma is the period, not the rate
    beta ~ truncated_gaussian(0.57586873, 0.3, lower = 0) 
    mu ~ truncated_gaussian(0.001, 0.001, lower = 0) 
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
    tau ~ truncated_gaussian(tau, 0.05, lower=0)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 1000000)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(sigma =0.19, gamma =0.2, beta=0.53, mu =0.0015)

bi <- sample(bi_model, end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'model') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 100, thin = 1) %>% # burn in 
  sample(nsamples = 100000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)
write.csv(bi_lst, "../data/short1.csv")
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
write.csv(fitY,"../data/short_y1.csv")

Mmodel <- read.csv("simulatestates.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[,3]
S<-S[1:30]
E<-Mmodel[,4]
E<-E[1:30]
I<-Mmodel[,5]
I<-I[1:30]
R<-Mmodel[,6]
R<-R[1:30]
M<-Mmodel[,7]
M<-M[1:30]
S <- data.frame(value = S) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitS <-bi_lst$S %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(S %>% rename(S = value))
write.csv(fitS,"../data/short_S1.csv")

E <- data.frame(value = E) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitE <-bi_lst$E %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(E %>% rename(E = value))
write.csv(fitE,"../data/short_E1.csv")

I <- data.frame(value = I) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitI <-bi_lst$I %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(I %>% rename(I = value))
write.csv(fitI,"../data/short_I1.csv")

R <- data.frame(value = R) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitR <-bi_lst$R %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(R %>% rename(R = value))
write.csv(fitR,"../data/short_R1.csv")

M <- data.frame(value = M) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitM <-bi_lst$M %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(M %>% rename(M = value))
write.csv(fitM,"../data/short_M1.csv")

write.csv(bi_lst$sigma$value,"../data/short_sigma1.csv")
write.csv(bi_lst$gamma$value,"../data/short_gamma1.csv")
write.csv(bi_lst$beta$value,"../data/short_beta1.csv")
write.csv(bi_lst$mu$value,"../data/short_mu1.csv")




