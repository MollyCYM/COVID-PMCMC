rm(list=ls())
set.seed(0066661)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
# Load the data
v <- read.csv("covidoudg2_y821w.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()
y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)

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
  state x
  state Z

  input N
  input Forcing
  
  const k = 5
  const gamma = 9
  const theta = 0.2
  const a = -0.02
  const b = -0.2
  
  param tau

  
  sub parameter {
    tau ~ truncated_gaussian(0.1, 0.1, lower = 0)
  }
  
  sub proposal_parameter {
    tau ~ gaussian(tau, 0.001)
  }
  
  sub initial {
    x <- a
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
    Z <- 1/k
  }

  sub transition(delta = 1) {
  Z <- ((t_now) % 7 == 0 ? 0 : Z)
    inline mu = a+b*Forcing
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = theta*(mu-x)
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
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(tau=0.1)
#LibBi wrapper 
#run launches LibBi with a particular set of command line arguments
bi_model <- libbi(model,end_time = end_time, input = input_lst, 
                  init=init_list, obs = obs_lst)
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
  sample(nsamples = 10000, thin = 1, init=init_list)

bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/para5_model1131433.csv")
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
write.csv(fitY,"../data/para5_y1131433.csv")

plot_df <- bi_lst$x |> mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(plot_df,"../data/para5_beta1131433.csv")

Mmodel <- read.csv("covidoudg2_model821.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[-1,5]
E<-Mmodel[-1,7]
I<-Mmodel[-1,9]
R<-Mmodel[-1,11]

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
write.csv(fitS,"../data/para5_S1131433.csv")

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
write.csv(fitE,"../data/para5_E1131433.csv")

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
write.csv(fitI,"../data/para5_I1131433.csv")

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
write.csv(fitR,"../data/para5_R1131433.csv")



write.csv(bi_lst$tau$value,"../data/para5_tau1131433.csv")








