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
v <- read.csv("Eng_308dy.csv", header=FALSE, stringsAsFactors=FALSE) 

y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time","value")

L <- read.csv("Eng_308forcing.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
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
  state x

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
    theta ~ truncated_gaussian(0.05, 0.2, lower = 0)
    a ~ truncated_gaussian(-0.02, 0.2, upper = 0)
    b ~ truncated_gaussian(-0.2, 0.5, upper = 0)
  }
  
 sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.001, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.001, lower = 0) 
    sigma ~ truncated_gaussian(sigma, 0.0001, lower = 0)
    theta ~ truncated_gaussian(theta, 0.0001, lower = 0)
    a ~ gaussian(a, 0.0001)
    b ~ gaussian(b, 0.0001)
 }
  
  sub initial {
    x ~ gaussian(a, sigma/sqrt(2*theta) )
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
  }

  sub transition(delta = 1) {
    noise e
    e ~ wiener()
    inline mu = a+b*Forcing
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = theta*(mu-x)+sigma*e
      dS/dt = -exp(x)*S*(0.1*I+E)/N
      dE/dt = exp(x)*S*(0.1*I+E)/N - E*(1/k+1/gamma)
      dI/dt = E/k-I*(1/gamma+0.0087)
      dR/dt = (I+E)/gamma+0.0087*I
    }
  }

  sub observation {
    y ~ poisson((E/k)/5)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])

input_lst <- list(N = 52196380,Forcing=Forcing)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k=5, gamma=9, sigma=sqrt(0.004), theta=0.05, a=-0.02, b=-0.2)
#LibBi wrapper 
#run launches LibBi with a particular set of command line arguments
bi_model <- libbi(model,end_time = end_time, input = input_lst, 
                  init=init_list,obs = obs_lst)
#RBi.helpers adapt_particle
particles_adapted <- bi_model %>%
  sample(nsamples = 2000, nparticles = minParticles, 
         nthreads = ncores, proposal = 'prior') %>%
  adapt_particles(min = minParticles, max = minParticles*500)

#RBi.helpers adapt_proposal
proposal_adapted <- particles_adapted %>%
  sample(target = "posterior", nsamples = 2000, 
         nthreads = ncores, proposal = 'model') %>%
  adapt_proposal(min = 0.1, max = 0.4)

#Running pMCMC with burn-in
bi <- proposal_adapted %>%
  sample(nsamples = 8000, thin = 1)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/para10_model2.csv")
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
write.csv(fitY,"../data/para10_y2.csv")

plot_df <- bi_lst$x |> mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df,"../data/para10_beta2.csv")

fitS <-bi_lst$S %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() 
write.csv(fitS,"../data/para10_S2.csv")

fitE <-bi_lst$E %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() 
write.csv(fitE,"../data/para10_E2.csv")

fitI <-bi_lst$I %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(fitI,"../data/para10_I2.csv")

fitR <-bi_lst$R %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() 
write.csv(fitR,"../data/para10_R2.csv")


write.csv(bi_lst$k$value,"../data/para10_alpha2.csv")
write.csv(bi_lst$gamma$value,"../data/para10_gamma2.csv")
write.csv(bi_lst$sigma$value,"../data/para10_sigma2.csv")
write.csv(bi_lst$theta$value,"../data/para10_theta2.csv")
write.csv(bi_lst$a$value,"../data/para10_a2.csv")
write.csv(bi_lst$b$value,"../data/para10_b2.csv")







