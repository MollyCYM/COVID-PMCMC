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
v <- read.csv("SEIR_366.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time", "value")

ncores <- 8
minParticles <- max(ncores, 16)
model_str <-
"
model SEIR {
  obs y
  
  state S
  state E
  state I
  state R
  
  input N
  
  param sigma
  param beta
  param gamma
  param mu
  param tau
  
  sub parameter {
    sigma ~ truncated_gaussian(0.23, 0.05, lower = 0) //
    gamma ~ truncated_gaussian(0.12460946, 0.2, lower = 0) // gamma is the period, not the rate
    beta ~ truncated_gaussian(0.57586873, 0.3, lower = 0)
    mu ~ truncated_gaussian(0.001, 0.001, lower = 0)
    tau ~ uniform(0,1)
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
      dR/dt = (gamma+mu)*I
    }
  }
  
  sub observation {
    y ~ log_normal(log(max((sigma*E)/5, 0)), tau)
  }
  
  sub proposal_parameter {
    sigma ~ truncated_gaussian(sigma, 0.001,lower=0)
    gamma ~ truncated_gaussian(gamma, 0.001, lower=0)
    beta ~ truncated_gaussian(beta, 0.01, lower=0)
    mu ~ truncated_gaussian(mu,0.0001, lower=0)
    tau ~ truncated_gaussian(tau,0.0001, lower=0)
  }
}" 
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)


input_lst <- list(N =1000000 )
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(sigma =0.19, gamma =0.2, beta=0.53, mu =0.0015)
bi <- sample(bi_model, end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 2000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>%                                 
adapt_particles(min = minParticles, max = minParticles*500) %>%
adapt_proposal(min = 0.05, max = 0.4) %>%
sample(nsamples = 1000, thin = 1)

bi_lst <- bi_read(bi %>% sample_obs)
write.csv(bi_lst, "../data/para0_model1.csv")
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
write.csv(fitY,"../data/para0_y1.csv")

Mmodel <- read.csv("simulatestates4.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[,3]
E<-Mmodel[,4]
I<-Mmodel[,5]
R<-Mmodel[,6]
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
write.csv(fitS,"../data/para0_S1.csv") 

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
write.csv(fitE,"../data/para0_E1.csv")

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
write.csv(fitI,"../data/para0_I1.csv")

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
write.csv(fitR,"../data/para0_R1.csv")

write.csv(bi_lst$sigma$value,"../data/para0_sigma1.csv")
write.csv(bi_lst$gamma$value,"../data/para0_gamma1.csv") 
write.csv(bi_lst$beta$value,"../data/para0_beta1.csv")
write.csv(bi_lst$mu$value,"../data/para0_mu1.csv")
    