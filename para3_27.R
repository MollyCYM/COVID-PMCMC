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
v <- read.csv("H1N1BM_wk3.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()

y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
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

  state Z

  input N
  param k
  param gamma
  param sigma // Noise driver
  param E0
  param I0
  param R0
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(1.59, 0.2, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(1.08, 0.75, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    x0 ~ uniform(-5,2)
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
    x <- x0
    Z <- 0
  }

  sub transition(delta = 1) {
  Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

  sub proposal_parameter {
    k ~ gaussian(k, 0.005)
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
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
init_list <- list(sigma=0.07, gamma=1.08, k=1.59, tau=0.1)


bi <- sample(bi_model, end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 2000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=123) %>% 
  adapt_particles(min = minParticles, max = minParticles*500) %>%
  adapt_proposal(min = 0.1, max = 0.4) %>%
  sample(nsamples =10000, thin = 1)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/para3_model27.csv")
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
write.csv(fitY,"../data/para3_y27.csv")

plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df,"../data/para3_beta27.csv")

Mmodel <- read.csv("simulateh1n1states1.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[,4]
E<-Mmodel[,5]
I<-Mmodel[,6]
R<-Mmodel[,7]

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
write.csv(fitS,"../data/para3_S27.csv")

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
write.csv(fitE,"../data/para3_E27.csv")

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
write.csv(fitI,"../data/para3_I27.csv")

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
write.csv(fitR,"../data/para3_R27.csv")


write.csv(1/bi_lst$k$value,"../data/para3_alpha27.csv")
write.csv(1/bi_lst$gamma$value,"../data/para3_gamma27.csv")
write.csv(bi_lst$sigma$value,"../data/para3_sigma27.csv")
write.csv(bi_lst$tau$value,"../data/para3_tau27.csv")


