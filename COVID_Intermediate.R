rm(list=ls())
require(deSolve)
set.seed(34)   #random seed 32 34 36 37 38 39
times <- 1:365
N=52196381
#Simulate a Brownian Motion Path
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
e <- c(0, cumsum(e))

mu<-vector(length=365) 
for (t in 1:365) {
  if (t <= 120){
    mu[t]=-0.02      }                    #No lock-down policy
  else {mu[t]=-0.2-0.02 }                 #Lock-down policy
}
mu<-ts(mu)

plot(mu,type='l')
#Simulate an O-U Process
ornstein_uhlenbeck <- function(n,theta,sigma,x0){
  dt  <- 1
  x<-vector(length=366)
  for (i in 1:(n+1)) {
    if (i==1){x[i]=x0}
    else{
      x[i]  <-  x[i-1] + theta*(mu[i-1]-x[i-1])*dt + sigma*e[i-1]}
  }
  return(x);
}

x<- ornstein_uhlenbeck(365,0.05,sqrt(0.004),1)
plot(x,type='l',xlab="time")
abline(v=121, col="red")
lines(exp(mu+0.02),col="blue")
write.csv(x,"covidintex1.csv")
#Main ODE Model
COVID_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
    dt <- 1
    dS <- -exp(x[t+1])*S*I/N
    dE <- exp(x[t+1])*S*I/N - E/k
    dI <- E/k-I*(1/gamma+0.0087)
    dR <- I/gamma+0.0087*I

    return(list(c( dt, dS, dE, dI, dR)))
  })
}
params <- c(k=7, gamma=5)   #Trial: k=5 & gamma=5; k=7 & gamma=5

#library('truncnorm')
# R0 <- rtruncnorm(1, a=0, b=1, mean = 0.15, sd = 0.15)
# E0 <-runif(1,-16, -9)
# I0 <-runif(1,-16, -9)
# #x <-runif(1, -5,2)
# S <- N
# R <- R0*S
# S <- S - R
# 
# E <- exp(E0 + log(S))
# S <- S - E
# I <- exp(I0 + log(S))
# S <- S - I

R0 <-0.03
S <- N
R <- S*R0
S <- S - R

E0 <- -15
E <- exp(E0 + log(S))
S <- S - E

I0 <- -10
I <- exp(I0 + log(S))
S <- S - I
initial_state<- c(t=1, S=S, E=E, I=I, R=R)
#initial_state <- c(t=1,S=52196380, E=1, I=0, R=0)
model3 <- ode(initial_state, times, COVID_OU, params)

summary(model3)

matplot(model3, type="l", lty=1, main="COVID—OU Model",ylab="counts", xlab="Time")
legend <- colnames(model3)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z1 <-model3[,4]/7

tau1 <- 0.8#runif(1,0,1)
Y1 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}
plot(Y1,type='l',ylab="Obs Y",xlab="time",col="blue")
write.csv(Y1,"covidinterY1.csv")
write.csv(model3,"simcovidinter1.csv")

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
v <- read.csv("covidinter1.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time", "value")
L <- read.csv("Forcing.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(Forcing) <- c("time","value")

# Forcing2 <- data.frame(value = L) %>%
#   mutate(time = seq(1, by = 1, length.out = n())) %>%
#   dplyr::select(time, F)
# colnames(Forcing2) <- c("time","value")
#Forcing <- as.numeric(unlist(Forcing))
ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model dureau {
  obs y
  
  const h = 1
  state S
  state E
  state I
  state R
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
  param E0
  param I0
  param R0
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(7, 0.1, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(5, 0.1, lower = 0) // gamma is the period, not the rate
    sigma ~ truncated_gaussian(sqrt(0.004), 0.0001, lower = 0)
    theta ~ truncated_gaussian(0.05, 0.01, lower = 0)
    a ~ gaussian(-0.02, 0.01)
    b ~ gaussian(-0.2, 0.1)
    x0 ~ uniform(0,2)
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
    x<- x0
    Z <- 0
  }

  sub transition(delta = h) {
    noise e
    e ~ wiener()
    
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = theta*(a+b*Forcing-x)+sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I*(1/gamma+0.0087)
      dR/dt = I/gamma+0.0087*I
      dZ/dt = E/k

    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

  sub proposal_parameter {
    k ~ gaussian(k, 0.01)
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    theta ~ truncated_gaussian(theta, 0.01,lower=0)
    a ~ gaussian(a, 0.01)
    b ~ gaussian(b, 0.01)
    x0 ~ gaussian(x0, 0.05)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
    tau ~ gaussian(tau, 0.05)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N=52196381,Forcing=Forcing)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=1234) %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 100, thin = 1) %>% # burn in 
  sample(nsamples = 2000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)
#dx/dt = theta*(a+b*Forcing-x)+sigma*e/h


write.csv(bi_lst, "../data/covid365_1.csv")
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
write.csv(fitY,"../data/covid365_y1.csv")

btr <- read.csv("covidintex1.csv", header=TRUE, stringsAsFactors=FALSE)
bt <- exp(btr[,2])
B <- data.frame(value = bt) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(B %>% rename(B = value))
write.csv(plot_df,"../data/covid365_beta1.csv")

# plot_df1 <- bi_lst$x %>% mutate(value = exp(value)) %>%
#   group_by(np) %>% mutate(value = value - value[1]) %>%
#   group_by(time) %>%
#   mutate(
#     q025 = quantile(value, 0.025),
#     q25 = quantile(value, 0.25),
#     q50 = quantile(value, 0.5),
#     q75 = quantile(value, 0.75),
#     q975 = quantile(value, 0.975)
#   ) %>% ungroup()
# write.csv(plot_df1,"covid259_beta02.csv")

Mmodel <- read.csv("simcovidinter1.csv", header=TRUE, stringsAsFactors=FALSE)
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
write.csv(fitS,"../data/covid365_S1.csv")

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
write.csv(fitE,"../data/covid365_E1.csv")

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
write.csv(fitI,"../data/covid365_I1.csv")

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
write.csv(fitR,"../data/covid365_R1.csv")



write.csv(bi_lst$k$value,"../data/covid365_alpha1.csv")
write.csv(bi_lst$gamma$value,"../data/covid365_gamma1.csv")
write.csv(bi_lst$sigma$value,"../data/covid365_sigma1.csv")
write.csv(bi_lst$theta$value,"../data/covid365_theta1.csv")
