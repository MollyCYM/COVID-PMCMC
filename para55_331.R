rm(list=ls())
set.seed(0066666662)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
# Load the data
v <- read.csv("covidoudg2_y121d.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time","value")

ncores <- 12
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
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
    k ~ truncated_gaussian(5, 2, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(9, 2, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(sqrt(0.004), 0.01, lower = 0)
    tau ~ truncated_gaussian(0.1, 0.01, lower = 0)
  }
  
  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.005, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.005, lower = 0) 
    sigma ~ truncated_gaussian(sigma, 0.0001, lower = 0)
    tau ~ truncated_gaussian(tau, 0.0001,lower = 0)
  }
  
  sub initial {
    x ~ gaussian(0, 0.3) 
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
  }

  sub transition(delta = 1) {
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I*(1/gamma+0.0087)
      dR/dt = I*(1/gamma+0.0087)
    }
  }

  sub observation {
    y ~ log_normal(log(max((E/k)/5, 0)), tau)
  }
}"

model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])

input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k=5, gamma=9, sigma=sqrt(0.004),tau=0.1)
#LibBi wrapper 
#run launches LibBi with a particular set of command line arguments
bi_model <- libbi(model,end_time = end_time, input = input_lst, 
                  init=init_list, obs = obs_lst)
#RBi.helpers adapt_particle
particles_adapted <- bi_model %>%
  sample(nsamples = 1000, nparticles = minParticles, 
         nthreads = ncores, proposal = 'prior') %>%
  adapt_particles(min = minParticles, max = minParticles*500)

#RBi.helpers adapt_proposal
proposal_adapted <- particles_adapted %>%
  sample(target = "posterior", nsamples = 1000, 
         nthreads = ncores, proposal = 'model') %>%
  adapt_proposal(min = 0.2, max = 0.4)

#Running pMCMC with burn-in
bi <- proposal_adapted %>%
  sample(nsamples = 1000, thin = 1,init=init_list) %>%
  sample(nsamples = 10000, thin = 1)
bi_lst <- bi_read(bi %>% sample_obs)

fitY <- bi_lst$y %>%
  group_by(time) %>% 
  ungroup() %>%
  left_join(y %>% rename(Y = value))
write.csv(fitY,"../data/para55_y331.csv")

plot_df <- bi_lst$x  %>%
  group_by(time) 
write.csv(plot_df,"../data/para55_x331.csv")

Mmodel <- read.csv("covidoudg2_model121.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[-1,7]
E<-Mmodel[-1,9]
I<-Mmodel[-1,11]
R<-Mmodel[-1,13]

S <- data.frame(value = S) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitS <-bi_lst$S %>%
  group_by(time) %>%
  ungroup() %>%
  left_join(S %>% rename(S = value))
write.csv(fitS,"../data/para55_S331.csv")

E <- data.frame(value = E) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitE <-bi_lst$E %>%
  group_by(time) %>%
  ungroup() %>%
  left_join(E %>% rename(E = value))
write.csv(fitE,"../data/para55_E331.csv")

I <- data.frame(value = I) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitI <-bi_lst$I %>%
  group_by(time) %>%
  ungroup() %>%
  left_join(I %>% rename(I = value))
write.csv(fitI,"../data/para55_I331.csv")

R <- data.frame(value = R) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
fitR <-bi_lst$R %>%
  group_by(time) %>%
  ungroup() %>%
  left_join(R %>% rename(R = value))
write.csv(fitR,"../data/para55_R331.csv")

write.csv(bi_lst$loglikelihood$value,"../data/para55_loglik331.csv")
write.csv(bi_lst$k$value,"../data/para55_alpha331.csv")
write.csv(bi_lst$gamma$value,"../data/para55_gamma331.csv")
write.csv(bi_lst$sigma$value,"../data/para55_sigma331.csv")
write.csv(bi_lst$tau$value,"../data/para55_tau331.csv")






