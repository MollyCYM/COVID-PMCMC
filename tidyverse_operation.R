rm(list=ls())
# library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
# Load the data
v <- read.csv("covidoudg2_y2w.csv", header=FALSE, stringsAsFactors=FALSE) |>
  rowSums()

y <- data.frame(value = v) |>
  # mutate(time = seq(7, by = 7, length.out = n())) |>
  mutate(time = seq(7, by = 7, length.out = 52)) |>
  dplyr::select(time, value)
L <- read.csv("Forcing.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) |>
  mutate(time = seq(1, by = 1, length.out = 365)) |>
  dplyr::select(time,V1 )
colnames(Forcing) <- c("time","value")

ncores <- 12
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
  
  sub parameter {
    k ~ truncated_gaussian(5, 1, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(9, 1, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(sqrt(0.004), 0.1, lower = 0)
    theta ~ truncated_gaussian(0.05, 0.1, lower = 0)
    a ~ truncated_gaussian(-0.02, 0.05, upper = 0)
    b ~ truncated_gaussian(-0.2, 0.1, upper = 0)
  }

  sub initial {
    x ~ gaussian(-0.02, 0.2)
    
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
    Z <- 1/k
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
    y ~ poisson(Z/5)
  }

  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.01, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.01, lower = 0) 
    sigma ~ truncated_gaussian(sigma, 0.001, lower = 0)
    theta ~ truncated_gaussian(theta, 0.001, lower = 0)
    a ~ gaussian(a, 0.001)
    b ~ gaussian(b, 0.001)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 52196381,Forcing=Forcing)
end_time <- max(y$time)
obs_lst <- list(y = y |> dplyr::filter(time <= end_time))
init_list <- list(k=6, gamma=10, sigma=0.07,theta=0.07,a=-0.04,b=-0.4)

bi <- sample(bi_model,target = "posterior", end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 2000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=00000888) |> 
  adapt_particles(min = minParticles, max = minParticles*4) |>
  adapt_proposal(min = 0.05, max = 0.4) |>
  sample(nsamples = 20, thin = 1, init=init_list)

bi_lst <- bi_read(bi |> sample_obs)
bi_lst <- bi_read(bi |> sample_obs())

Mmodel <- read.csv("covidoudg2_model2.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[-1,7]
E<-Mmodel[-1,9]
I<-Mmodel[-1,11]
R<-Mmodel[-1,13]

S <- data.frame(value = S) |>
  mutate(time = seq(1, by = 1, length.out = 365)) |>
  dplyr::select(time, value)
fitS <-bi_lst$S |>
  group_by(time) |>
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) |> ungroup() |>
  left_join(S |> rename(S = value))
write.csv(fitS,"para5_S3121211.csv")

