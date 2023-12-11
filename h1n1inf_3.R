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
v <- read.csv("h1n1bm_ydaily1.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time","value")

ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model h1n1bm {
  obs y

  state S
  state E
  state I
  state R
  state x

  input N
  
  param k
  param gamma
  param sigma // Noise driver
  param tau

  sub parameter {
    k ~ truncated_gaussian(2, 3, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(5, 5, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(0.1, 0.1, lower = 0)
    tau ~ truncated_gaussian(0.5, 0.1, lower = 0)
  }
  
  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.01, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.01, lower = 0) 
    sigma ~ truncated_gaussian(sigma, 0.001, lower = 0)
    tau ~ gaussian(tau, 0.001)
  }

  sub initial {
    S <- N-1
    E <- 0
    I <- 1
    R <- 0
  }

  sub transition(delta = 1) {
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
    }
  }

  sub observation {
    y ~ log_normal(log(max((E/k)/10, 0)), tau)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k=2, gamma=8, sigma=0.07, tau=0.1)
#LibBi wrapper 
#run launches LibBi with a particular set of command line arguments
bi_model <- libbi(model,end_time = end_time, input = input_lst, 
                  init=init_list,obs = obs_lst)
#RBi.helpers adapt_particle
particles_adapted <- bi_model %>%
  sample(nsamples = 2000, nparticles = minParticles, 
         nthreads = ncores, proposal = 'prior') %>%
  adapt_particles(min = minParticles, max = minParticles*500)

#RBi.helpers adapt_proposal
proposal_adapted <- particles_adapted %>%
  sample(target = "posterior", nsamples = 2000, 
         nthreads = ncores, proposal = 'model') %>%
  adapt_proposal(min = 0.1, max = 0.4)

#Running pMCMC with burn-in
bi <- proposal_adapted %>%
  sample(nsamples = 5000, thin = 1,init=init_list)

bi_lst <- bi_read(bi %>% sample_obs)

#Data save command
write.csv(bi_lst,"../data/h1n1inf_model3.csv")
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
write.csv(fitY,"../data/h1n1inf_y3.csv")

plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df,"../data/h1n1inf_beta3.csv")

Mmodel <- read.csv("h1n1bm_model1.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[-1,7]
E<-Mmodel[-1,9]
I<-Mmodel[-1,11]
R<-Mmodel[-1,13]

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
write.csv(fitS,"../data/h1n1inf_S3.csv")

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
write.csv(fitE,"../data/h1n1inf_E3.csv")

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
write.csv(fitI,"../data/h1n1inf_I3.csv")

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
write.csv(fitR,"../data/h1n1inf_R3.csv")


write.csv(bi_lst$k$value,"../data/h1n1inf_alpha3.csv")
write.csv(bi_lst$gamma$value,"../data/h1n1inf_gamma3.csv")
write.csv(bi_lst$sigma$value,"../data/h1n1inf_sigma3.csv")
write.csv(bi_lst$tau$value,"../data/h1n1inf_tau3.csv")


