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
v <- read.csv("covidinter1.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time", "value")
L <- read.csv("Forcing.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) %>%
mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(Forcing) <- c("time","value")

# Forcing2 <- data.frame(value = L) %>%
#   mutate(time = seq(1, by = 1, length.out = n())) %>%
#   dplyr::select(time, F)
# colnames(Forcing2) <- c("time","value")
#Forcing <- as.numeric(unlist(Forcing))
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
  state mu
  state x

  state Z
  
  input N
  input Forcing
  
  param k
  param gamma
  param sigma // Noise driver
  param theta
  param a
  param b
  param tau

  sub parameter {
    k ~ truncated_gaussian(7, 0.1, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(5, 0.1, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(sqrt(0.004), 0.0001, lower = 0)
    theta ~ truncated_gaussian(0.05, 0.001, lower = 0)
    a ~ gaussian(-0.02, 0.01)
    b ~ gaussian(-0.2, 0.01)
    tau ~ uniform(0,1)
  }

  sub initial {
    S <- N
    R <- 0.03*S
    S <- S - R

    E <- exp(-15 + log(S))
    S <- S - E
    I <- exp(-10 + log(S))
    S <- S - I
    x<- 1
    mu<- 1
    Z <- 0
  }

  sub transition(delta = h) {
    noise e
    e ~ wiener()
    mu <- a+b*Forcing
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
    
      dx/dt = theta*(mu-x)+sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I*(1/gamma+0.0087)
      dR/dt = I/gamma+0.0087*I
      dZ/dt = E/k

    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

  sub proposal_parameter {
    k ~ gaussian(k, 0.01)
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    theta ~ truncated_gaussian(theta, 0.001,lower=0)
    a ~ gaussian(a, 0.01)
    b ~ gaussian(b, 0.01)
    tau ~ gaussian(tau, 0.01)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N=52196381,Forcing=Forcing)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k =7, gamma =5, sigma=sqrt(0.004), theta =0.05, a=-0.02, b=-0.2, tau=0.8)

bi <- sample(bi_model, end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=12) %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 2000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)
#dx/dt = theta*(a+b*Forcing-x)+sigma*e/h


write.csv(bi_lst, "../data/covid365_1.csv")
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
write.csv(fitY,"../data/covid365_y1.csv")


plot_df <- bi_lst$x %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) 
write.csv(plot_df,"../data/covid365_beta1.csv")

fitmu <-bi_lst$mu %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) 
write.csv(fitmu,"../data/covid365_mu1.csv")

Mmodel <- read.csv("simcovidinter1.csv", header=TRUE, stringsAsFactors=FALSE)
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
write.csv(fitS,"../data/covid365_S1.csv")

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
write.csv(fitE,"../data/covid365_E1.csv")

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
write.csv(fitI,"../data/covid365_I1.csv")

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
write.csv(fitR,"../data/covid365_R1.csv")



write.csv(bi_lst$k$value,"../data/covid365_alpha1.csv")
write.csv(bi_lst$gamma$value,"../data/covid365_gamma1.csv")
write.csv(bi_lst$sigma$value,"../data/covid365_sigma1.csv")
write.csv(bi_lst$theta$value,"../data/covid365_theta1.csv")
