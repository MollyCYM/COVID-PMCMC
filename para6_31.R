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
v <- read.csv("covidoudg2_y4w.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()

y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
L <- read.csv("Forcing.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(Forcing) <- c("time","value")

x_read <- read.csv("covidoudg2_x41.csv", header=FALSE, stringsAsFactors=FALSE)
x <- data.frame(value = x_read) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(x) <- c("time","value")

ncores <- 15
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
  obs y
  
  state S
  state E
  state I
  state R
  state mu
  state x

  state Z

  input N
  input Forcing
  input x
  
  const sigma = 0.000001
  
  param k
  param gamma
  param theta
  param a
  param b
  param tau
  param x0
  
  sub parameter {
    k ~ truncated_gaussian(5, 1, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(9, 1, lower = 0) // gamma is the period, not the rate
    theta ~ truncated_gaussian(0.05, 0.1, lower = 0)
    tau ~ truncated_gaussian(0.1, 0.05, lower = 0)
    a ~ truncated_gaussian(-0.02, 0.05, upper = 0)
    b ~ truncated_gaussian(-0.2, 0.1, upper = 0)
    x0 ~ gaussian(-0.02, 0.2) //long-term ou mean & sd
  }

  sub initial {
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
    Z <- 0
    x <- x0
  }

  sub transition(delta = 1) {
  Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    mu <- a+b*Forcing
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dS/dt = -exp(x)*S*(0.1*I+E)/N
      dE/dt = exp(x)*S*(0.1*I+E)/N - E*(1/k+1/gamma)
      dI/dt = E/k-I*(1/gamma+0.0087)
      dR/dt = (I+E)/gamma+0.0087*I
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.01, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.01, lower = 0) 
    theta ~ truncated_gaussian(theta, 0.001, lower = 0)
    tau ~ gaussian(tau, 0.001)
    a ~ gaussian(a, 0.001)
    b ~ gaussian(b, 0.001)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 52196381,Forcing=Forcing)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k=5, gamma=9,theta=0.05,tau=0.1,a=-0.02,b=-0.2)

bi <- sample(bi_model,target = "posterior", end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 2000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=0066661) %>% 
  adapt_particles(min = minParticles, max = minParticles*500) %>%
  adapt_proposal(min = 0.1, max = 0.4) %>%
  sample(nsamples = 10000, thin = 1)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/para6_model31.csv")
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
write.csv(fitY,"../data/para6_y31.csv")

plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df,"../data/para6_beta31.csv")

Mmodel <- read.csv("covidoudg2_model4.csv", header=TRUE, stringsAsFactors=FALSE)
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
write.csv(fitS,"../data/para6_S31.csv")

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
write.csv(fitE,"../data/para6_E31.csv")

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
write.csv(fitI,"../data/para6_I31.csv")

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
write.csv(fitR,"../data/para6_R31.csv")


write.csv(1/bi_lst$k$value,"../data/para6_alpha31.csv")
write.csv(1/bi_lst$gamma$value,"../data/para6_gamma31.csv")
write.csv(bi_lst$tau$value,"../data/para6_tau31.csv")
write.csv(bi_lst$theta$value,"../data/para6_theta31.csv")
write.csv(bi_lst$a$value,"../data/para6_a31.csv")
write.csv(bi_lst$b$value,"../data/para6_b31.csv")







