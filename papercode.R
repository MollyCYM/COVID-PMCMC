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
v <- read.csv("covidposter_wk.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()

y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
L <- read.csv("Forcing_364.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(Forcing) <- c("time","value")

ncores <- 8
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
  
  param k
  param gamma
  param sigma // Noise driver
  param theta
  param a
  param b
  param E0
  param I0
  param R0
  param tau
  
  sub parameter {
    k ~ truncated_gaussian(5, 0.01, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(5, 0.01, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(sqrt(0.004), 0.00001, lower = 0)
    theta ~ truncated_gaussian(0.05, 0.0001, lower = 0)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
    tau ~ uniform(0, 1)
    a ~ gaussian(-0.02, 0.01)
    b ~ gaussian(-0.2, 0.01)
  }

  sub initial {
    S <- N
    R <- R0*S
    S <- S - R

    E <- exp(E0 + log(S))
    S <- S - E
    I <- exp(I0 + log(S))
    S <- S - I
    x <- log(0.8)
    Z <- 0
  }

  sub transition(delta = 1) {
  Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    mu <- a+b*Forcing
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = theta*(mu-x)+sigma*e
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
    k ~ truncated_gaussian(k, 0.001, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.001, lower = 0) 
    sigma ~ truncated_gaussian(sigma, 0.001, lower = 0)
    theta ~ truncated_gaussian(theta, 0.001, lower = 0)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
    tau ~ gaussian(tau, 0.05)
    a ~ gaussian(a, 0.01)
    b ~ gaussian(b, 0.01)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=1111) %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 5000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/Exp6_model2.csv")
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
write.csv(fitY,"../data/Exp6_y2.csv")

plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df,"../data/Exp6_beta2.csv")

plot_df1 <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(np) %>% mutate(value = value - value[1]) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df1,"../data/Exp6_beta02.csv")

fitmu <-bi_lst$mu %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) 
write.csv(fitmu,"../data/Exp6_mu2.csv")

Mmodel <- read.csv("Covidou1.csv", header=TRUE, stringsAsFactors=FALSE)
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
write.csv(fitS,"../data/Exp6_S2.csv")

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
write.csv(fitE,"../data/Exp6_E2.csv")

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
write.csv(fitI,"../data/Exp6_I2.csv")

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
write.csv(fitR,"../data/Exp6_R2.csv")

write.csv(1/bi_lst$k$value,"../data/Exp6_alpha2.csv")
write.csv(1/bi_lst$gamma$value,"../data/Exp6_gamma2.csv")
write.csv(bi_lst$sigma$value,"../data/Exp6_sigma2.csv")
write.csv(bi_lst$theta$value,"../data/Exp6_theta2.csv")


