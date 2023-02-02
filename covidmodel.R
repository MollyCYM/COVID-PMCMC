rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
times <- 1:365
N=56000000     #England population size in 2023


mu<-vector(length=365) 
for (t in 1:365) {
  if (t <= 60){
    mu[t]=exp(-0.92)      }                    #No lock-down policy
  else if (t <=280){                         
    mu[t]=exp(-0.92-1.38)     }                #Strict national lock-down 
  else {mu[t]=exp(-0.92-0.69) }                #Mitigate Tier System
}
mu<-ts(mu)
# Load the data
v <- read.csv("simcoviddata.csv", header=FALSE, stringsAsFactors=FALSE) 
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
  state x
  
  state Z

  input N
  input I
  
  param alpha
  param gamma
  param sigma // Noise driver
  param E0
  param I0
  param R0
  param x0
  param tau

  sub parameter {
    alpha ~ truncated_gaussian(0.125, 0.05, lower = 0) 
    gamma ~ truncated_gaussian(0.07, 0.01, lower = 0) 
    sigma ~ uniform(0,1) // or uniform(0,0.1)
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
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = theta*(exp(b0+b1*I[t])-x) + sigma*e
      dS/dt = -x*S*(0.1*I+E)/N
      dE/dt = x*S*(0.1*I+E)/N - E*(alpha+gamma)
      dI/dt = E*alpha-I*(gamma+0.0087)
      dR/dt = (I+E)*gamma+0.0087*I
      dZ/dt = E*alpha
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

  sub proposal_parameter {
    alpha ~ gaussian(alpha, 0.05)
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
input_lst <- list(N = 56000000, mu)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 10, thin = 5) %>% # burn in 
  sample(nsamples = 200, thin = 5) #20000, 40000

bi_lst <- bi_read(bi %>% sample_obs)
write.csv(bi_lst,"60w2.csv")
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
write.csv(fitY,"60wy2.csv")

plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df,"60wbeta2.csv")

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
write.csv(plot_df1,"60wbeta12.csv")

write.csv(bi_lst$k$value,"60wpalpha2.csv")
write.csv(1/bi_lst$k$value,"60walpha2.csv")
write.csv(bi_lst$gamma$value,"60wpgamma2.csv")
write.csv(1/bi_lst$gamma$value,"60wgamma2.csv")

