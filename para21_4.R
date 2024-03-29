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
xread <- read.csv("para21_4x.csv", header=FALSE, stringsAsFactors=FALSE)
x <- data.frame(value = xread) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(x) <- c("time","value")

# L <- read.csv("poster1_mu1.csv", header=TRUE, stringsAsFactors=FALSE)
# mu <- data.frame(value = L) %>%
#   mutate(time = seq(1, by = 1, length.out = n())) %>%
#   dplyr::select(time,value.x )
# colnames(mu) <- c("time","value")
# 
# ornstein_uhlenbeck <- function(n,theta,sigma,x0){
#   dt  <- 1
#   x<-vector(length=366)
#   for (i in 1:(n+1)) {
#     if (i==1){x[i]=x0}
#     else{
#       x[i]  <-  x[i-1] + theta*(mu[i-1]-x[i-1])*dt + sigma*e[i-1]}
#   }
#   return(x);
# }

ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
  obs y
  
  state S
  state E
  state I
  state R

  state Z

  input N
  
  param k
  param I0
  param x
  
  sub parameter {
    k ~ truncated_gaussian(5, 0.01, lower = 0) 
    I0 ~ truncated_gaussian(-11.3657, 0.001, lower = -11.38) 
  }

  sub initial {
    S <- N
    R <- 0.2439*S
    S <- S - R

    E <- exp(-14.7816 + log(S))
    S <- S - E
    I <- exp(I0 + log(S))
    S <- S - I
    Z <- 0
  }

  sub transition(delta = 1) {
  Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
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
    I0 ~ gaussian(I0, 0.005)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k=3)

bi <- sample(bi_model, end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=111123) %>% 
  adapt_particles(min = minParticles, max = minParticles*500) %>%
  adapt_proposal(min = 0.1, max = 0.4) %>%
  sample(nsamples = 1, thin = 1) %>% # burn in 
  sample(nsamples = 5000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/para21_model4.csv")
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
write.csv(fitY,"../data/para21_y4.csv")

# plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
#   group_by(time) %>%
#   mutate(
#     q025 = quantile(value, 0.025),
#     q25 = quantile(value, 0.25),
#     q50 = quantile(value, 0.5),
#     q75 = quantile(value, 0.75),
#     q975 = quantile(value, 0.975)
#   ) %>% ungroup()
# write.csv(plot_df,"../data/para21_beta4.csv")

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
write.csv(fitS,"../data/para21_S4.csv")

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
write.csv(fitE,"../data/para21_E4.csv")

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
write.csv(fitI,"../data/para21_I4.csv")

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
write.csv(fitR,"../data/para21_R4.csv")

write.csv(1/bi_lst$k$value,"../data/para21_alpha4.csv")


