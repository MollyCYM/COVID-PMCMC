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
v <- read.csv("covidoudg_y1w.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()

y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
L <- read.csv("Forcing.csv", header=FALSE, stringsAsFactors=FALSE)
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
  param tau
  
  sub parameter {
    k ~ truncated_gaussian(5, 0.05, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(9, 0.09, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(sqrt(0.004), 0.001, lower = 0)
    theta ~ truncated_gaussian(0.05, 0.001, lower = 0)
    tau ~ truncated_gaussian(0.1, 0.001, lower = 0)
    a ~ gaussian(-0.02, 0.001)
    b ~ gaussian(-0.2, 0.01)
  }

  sub initial {
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
    Z <- 0
    x <- 0
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
    k ~ truncated_gaussian(k, 0.01, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.01, lower = 0) 
    sigma ~ truncated_gaussian(sigma, 0.001, lower = 0)
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
init_list <- list(k=5, gamma=9, sigma=sqrt(0.004),theta=0.05,tau=0.1,a=-0.02,b=-0.2)

bi <- sample(bi_model,target = "posterior", end_time = end_time, input = input_lst, init=init_list,
             obs = obs_lst, nsamples = 1000,nparticles = minParticles,
             nthreads = ncores, proposal = 'model',seed=0066661) %>% 
  adapt_particles(min = minParticles, max = minParticles*500) %>%
  adapt_proposal(min = 0.1, max = 0.4)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/para6_model1.csv")

write.csv(1/bi_lst$k$value,"../data/para6_alpha1.csv")
write.csv(1/bi_lst$gamma$value,"../data/para6_gamma1.csv")
write.csv(bi_lst$sigma$value,"../data/para6_sigma1.csv")
write.csv(bi_lst$tau$value,"../data/para6_tau1.csv")
write.csv(bi_lst$theta$value,"../data/para6_theta1.csv")
write.csv(bi_lst$a$value,"../data/para6_a1.csv")
write.csv(bi_lst$b$value,"../data/para6_b1.csv")







