rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
library(readr)
# Load the data
# v <- read.csv("covid259days2.csv", header=FALSE, stringsAsFactors=FALSE) %>%
#   rowSums()
# 
# y <- data.frame(value = v) %>%
#   mutate(time = seq(7, by = 7, length.out = n())) %>%
#   dplyr::select(time, value)
v <- read.csv("covid259days3.csv", header=FALSE, stringsAsFactors=FALSE) 
L <- read.csv("F.csv", header=TRUE, stringsAsFactors=FALSE)
Input <- data.frame(value = L)
colnames(Input) <- c("time", "value")
Input<-data.frame(value=Input$value)
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time", "value")
ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
dim i(size = 259)

  obs y[i]

  state S[i]
  state E[i]
  state I[i]
  state R[i]
  state x[i]
  state mu[i]

  state Z[i]

  input N
  
  param k
  param gamma
  param sigma // Noise driver
  param theta
  param a
  param b
  param E0
  param I0
  
  param R0
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(0.2, 0.01, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(5, 0.05, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    theta ~ uniform(0,1)
    a ~ truncated_gaussian(1, 0.03, lower = 0)
    b ~ uniform(-2,0)
    x0 ~ uniform(1,2)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
    tau ~ uniform(0, 1)
  }

  sub initial {
    S[i] <- N
    R[i] <- R0*S[i]
    S[i] <- S[i] - R[i]

    E[i] <- exp(E0 + log(S[i]))
    S[i] <- S[i] - E[i]
    I[i] <- exp(I0 + log(S[i]))
    S[i] <- S[i] - I[i]
    x[i] <- x0
    mu[i] <- 3
    Z[i] <- 0
  }

  sub transition(delta = 1) {
    noise e[i]
    e[i] ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dmu[i]/dt = exp(a + b)
      dx[i]/dt = x[i-1] + theta*(mu[i-1]-x[i-1]) + sigma*e[i-1]
      dS[i]/dt = -x[i]*S[i]*(E[i]+0.1*I[i])/N
      dE[i]/dt = x[i]*S[i]*(E[i]+0.1*I[i])/N - E[i]*(1/k+1/gamma)
      dI[i]/dt = E[i]/k-I[i]*(1/gamma+0.0087)
      dR[i]/dt = (I[i]+E[i])/gamma+0.0087*I[i]
      dZ[i]/dt = E[i]/k
    }
  }

  sub observation {
    y[i] ~ log_normal(log(Z[i]), tau)
  }

  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.005,lower=0)
    sigma ~ truncated_gaussian(sigma, 0.01,lower=0)
    theta ~ truncated_gaussian(theta, 0.01,lower=0)
    gamma ~ truncated_gaussian(gamma, 0.01,lower=0)
    a ~ truncated_gaussian(a, 0.01,lower=0)
    b ~ gaussian(b, 0.01)
    x0 ~ gaussian(x0, 0.05)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
    tau ~ truncated_gaussian(tau, 0.05,lower=0)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 10000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)

####################################################################################
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
library(readr)
library(ncdf4) 
bi_open("/Users/mollycui/Desktop/R script/Research 10-Epid PMCMC/input.nc",file="input")
# Load the data
v <- read.csv("covid259days3.csv", header=FALSE, stringsAsFactors=FALSE) 
L <- read.csv("F.csv", header=TRUE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L)
colnames(Forcing) <- c("time", "value")
Forcing<-data.frame(value=Forcing$value)
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
  state mu

  state Z
  
  param k
  param gamma
  param sigma // Noise driver
  param theta
  param a
  param b
  param E0
  param I0
  
  param R0
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(0.2, 0.01, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(5, 0.05, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    theta ~ uniform(0,1)
    a ~ truncated_gaussian(1, 0.03, lower = 0)
    b ~ uniform(-2,0)
    x0 ~ uniform(1,2)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
    tau ~ uniform(0, 1)
  }

  sub initial {
    S <- 52196381
    R <- R0*S
    S <- S - R

    E <- exp(E0 + log(S))
    S <- S - E
    I <- exp(I0 + log(S))
    S <- S - I
    x <- x0
    mu <- 3
    Z <- 0
  }

  sub transition(delta = 1) {
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dmu/dt = exp(a + b)
      dx/dt =  theta*(mu-x) + sigma*e
      dS/dt = -x*S*(E+0.1*I)/52196381
      dE/dt = x*S*(E+0.1*I)/52196381 - E*(1/k+1/gamma)
      dI/dt = E/k-I*(1/gamma+0.0087)
      dR/dt = (I+E)/gamma+0.0087*I
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(Z), tau)
  }

  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.005,lower=0)
    sigma ~ truncated_gaussian(sigma, 0.01,lower=0)
    theta ~ truncated_gaussian(theta, 0.01,lower=0)
    gamma ~ truncated_gaussian(gamma, 0.01,lower=0)
    a ~ truncated_gaussian(a, 0.01,lower=0)
    b ~ gaussian(b, 0.01)
    x0 ~ gaussian(x0, 0.05)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
    tau ~ truncated_gaussian(tau, 0.05,lower=0)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
end_time <- max(y$time)
#input_lst <- list(Forcing = Forcing)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time,  obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 10000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)
write.csv(bi_lst, "covid2592.csv")
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
write.csv(fitY,"covid259_y2.csv")

btr <- read.csv("covid259beta3.csv", header=TRUE, stringsAsFactors=FALSE)
bt <- btr[,2]
B <- data.frame(value = bt) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
plot_df <- bi_lst$x %>% mutate(value = exp(value)+1) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(B %>% rename(B = value))
write.csv(plot_df,"covid259_beta2.csv")

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
write.csv(plot_df1,"covid259_beta02.csv")

Mmodel <- read.csv("simulatecovid12.csv", header=TRUE, stringsAsFactors=FALSE)
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
write.csv(fitS,"covid259_S2.csv")

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
write.csv(fitE,"covid259_E2.csv")

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
write.csv(fitI,"covid259_I2.csv")

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
write.csv(fitR,"covid259_R2.csv")



write.csv(bi_lst$k$value,"covid259_alpha2.csv")
write.csv(bi_lst$gamma$value,"covid259_gamma2.csv")
###############################################################################
rm(list=ls())
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
library(readr)
library(ncdf4) 
#bi_open("/Users/mollycui/Desktop/R script/Research 10-Epid PMCMC/input.nc",file="input")
# Load the data
v <- read.csv("covid259days3.csv", header=FALSE, stringsAsFactors=FALSE) 
L <- read.csv("F.csv", header=TRUE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L)
colnames(Forcing) <- c("value")
#Forcing<-data.frame(value=Forcing$value)
# Creates a 2D array
df<- data.frame(time = 1:259, values = L$F)
df_new <- tidyr::pivot_wider(df, names_from = time, values_from = L$F)
#Forcing <- tidyr::pivot_wider(Forcing, names_from =time, values_from = value)
Forcing <-as.matrix(Forcing)
#Forcing <-t(Forcing)
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time", "value")
ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
dim i(size = 259)
  obs y[i]

  state S[i]
  state E[i]
  state I[i]
  state R[i]
  state x[i]
  state mu[i]

  state Z[i]

  input Forcing[i]
  
  param k
  param gamma
  param sigma // Noise driver
  param theta
  param a
  param b
  param E0
  param I0
  
  param R0
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(0.2, 0.01, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(5, 0.05, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    theta ~ uniform(0,1)
    a ~ truncated_gaussian(1, 0.03, lower = 0)
    b ~ uniform(-2,0)
    x0 ~ uniform(1,2)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
    tau ~ uniform(0, 1)
  }

  sub initial {
    S[i] <- 52196381
    R[i] <- R0*S[i]
    S[i] <- S[i] - R[i]

    E[i] <- exp(E0 + log(S[i]))
    S[i] <- S[i] - E[i]
    I[i] <- exp(I0 + log(S[i]))
    S[i] <- S[i] - I[i]
    x[i] <- x0
    mu[i] <- 3
    Z[i] <- 0
  }

  sub transition(delta = 1) {
    noise e[i]
    e[i] ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dmu[i]/dt = exp(a + b*Forcing[i])
      dx[i]/dt =  theta*(mu[i]-x[i]) + sigma*e[i]
      dS[i]/dt = -x[i]*S[i]*(E[i]+0.1*I[i])/52196381
      dE[i]/dt = x[i]*S[i]*(E[i]+0.1*I[i])/52196381 - E[i]*(1/k+1/gamma)
      dI[i]/dt = E[i]/k-I[i]*(1/gamma+0.0087)
      dR[i]/dt = (I[i]+E[i])/gamma+0.0087*I[i]
      dZ[i]/dt = E[i]/k
    }
  }

  sub observation {
    y[i] ~ log_normal(log(Z[i]), tau)
  }

  sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.005,lower=0)
    sigma ~ truncated_gaussian(sigma, 0.01,lower=0)
    theta ~ truncated_gaussian(theta, 0.01,lower=0)
    gamma ~ truncated_gaussian(gamma, 0.01,lower=0)
    a ~ truncated_gaussian(a, 0.01,lower=0)
    b ~ gaussian(b, 0.01)
    x0 ~ gaussian(x0, 0.05)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
    tau ~ truncated_gaussian(tau, 0.05,lower=0)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
end_time <- max(y$time)
input_lst <- list(Forcing = Forcing)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 10000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)
###########################################################################
