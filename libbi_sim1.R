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
v <- read.csv("covidou_wk3.csv", header=FALSE, stringsAsFactors=FALSE) %>%
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
  state x

  state Z

  input N
  
  param k
  param E0
  param Forcing
  param mu
  
  sub parameter {
    k ~ truncated_gaussian(5, 0.01, lower = 0) 
    E0 ~ truncated_gaussian(1, 0.001, lower = 0) 
  }

  sub initial {
    S <- N-E0
    E <- E0
    I <- 0
    R <- 0
    x <- log(0.8)
    Z <- 0
  }

  sub transition(delta = 1) {
  Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    mu <- -0.02-0.2*Forcing
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = 0.05*(mu-x)+sqrt(0.004)*e
      dS/dt = -exp(x)*S*(0.1*I+E)/N
      dE/dt = exp(x)*S*(0.1*I+E)/N - E*(1/k+1/5)
      dI/dt = E/k-I*(1/5+0.0087)
      dR/dt = (I+E)/5+0.0087*I
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), 0.8)
  }

  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.0001, lower = 0) 
    E0 ~ gaussian(E0, 0.0001)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k=5)
bi_model <- libbi(model, obs=obs_lst, end_time = end_time, input = input_lst, init=init_list )
sim <- rbi::simulate(bi_model,nsamples=10)
sim_res <- bi_read(sim, type="state")
ggplot(sim_res$x, aes(x=time, group=np))+
  geom_line(aes(y=value)) +
  ylab("x")
ggplot(sim_res$S, aes(x=time, group=np))+
  geom_line(aes(y=value)) +
  ylab("S")

