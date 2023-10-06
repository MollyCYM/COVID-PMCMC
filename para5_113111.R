rm(list=ls())
set.seed(07811690)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
# Load the data
L <- read.csv("forcing30.csv", header=FALSE, stringsAsFactors=FALSE)
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
    k ~ truncated_gaussian(5, 1, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(9, 1, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(sqrt(0.004), 0.1, lower = 0)
    theta ~ truncated_gaussian(0.05, 0.2, lower = 0)
    tau ~ truncated_gaussian(0.1, 0.1, lower = 0)
    a ~ truncated_gaussian(-0.02, 0.1, upper = 0)
    b ~ truncated_gaussian(-0.2, 0.2, upper = 0)
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
  
  sub initial {
    x ~ gaussian(a, sigma/sqrt(2*theta) ) 
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
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }
}"

model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
input_lst <- list(N = 52196381,Forcing=Forcing)
init_list <- list(k=3, gamma=7, sigma=0.04, theta=0.03, tau=0.08, b_0=-0.01, b_1=-0.1)

#One. Generate 10 weeks data:10*7 
covid_data <- generate_dataset(model, end_time = 52 * 7, noutputs = 52,input = input_lst, init=init_list)
#Save and read the libbi output as a dataframe
dataset <- bi_read(covid_data)
covid_model <- libbi(model,input = input_lst, init=init_list,obs = covid_data)
#RBi.helpers adapt_particle
particles_adapted <- covid_model %>%
  sample(nsamples = 2000, nparticles = minParticles, 
         nthreads = ncores, proposal = "prior") %>%
  adapt_particles(min = minParticles, max = minParticles*500)

#RBi.helpers adapt_proposal
proposal_adapted <- particles_adapted %>%
  sample(target = "posterior", proposal = "model") %>%
  adapt_proposal(min = 0.1, max = 0.4,adapt="both")

#Running pMCMC with burn-in
bi <- proposal_adapted %>%
  sample(nsamples = 10000, thin = 1, init=init_list)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/para5_model113111.csv")
fitY <- bi_lst$y %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(dataset$y %>% rename(Y = value))
write.csv(fitY,"../data/para5_y113111.csv")

plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
left_join(dataset$x %>% rename(x = value))
write.csv(plot_df,"../data/para5_beta113111.csv")

fitS <-bi_lst$S %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(dataset$S %>% rename(S = value))
write.csv(fitS,"../data/para5_S113111.csv")

fitE <-bi_lst$E %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(dataset$E %>% rename(E = value))
write.csv(fitE,"../data/para5_E113111.csv")

fitI <-bi_lst$I %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(dataset$I %>% rename(I = value))
write.csv(fitI,"../data/para5_I113111.csv")

fitR <-bi_lst$R %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(dataset$R %>% rename(R = value))
write.csv(fitR,"../data/para5_R113111.csv")


write.csv(bi_lst$k$value,"../data/para5_alpha113111.csv")
write.csv(bi_lst$gamma$value,"../data/para5_gamma113111.csv")
write.csv(bi_lst$sigma$value,"../data/para5_sigma113111.csv")
write.csv(bi_lst$tau$value,"../data/para5_tau113111.csv")
write.csv(bi_lst$theta$value,"../data/para5_theta113111.csv")
write.csv(bi_lst$a$value,"../data/para5_a113111.csv")
write.csv(bi_lst$b$value,"../data/para5_b113111.csv")







