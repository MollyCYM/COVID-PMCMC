rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)

L <- read.csv("Forcing.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(Forcing) <- c("time","value")
e1 <- read.csv("covidoudg_e1.csv", header=FALSE, stringsAsFactors=FALSE)
e <- data.frame(value = e1) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(e) <- c("time","value")
ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model covidou {
  obs y
  
  state S
  state E
  state I
  state R
  state mu
  state x

  state Z

  input Forcing
  input e
  
  const N = 52196381
  const k = 5
  const gamma = 9
  const sigma = sqrt(0.004)
  const theta = 0.05
  const a = -0.02
  const b = -0.2
  const tau = 0.1

  sub initial {
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
    Z <- 0
    x <- 0
  }

  sub transition(delta = 1) {
  Z <- ((t_now) % 7 == 0 ? 0 : Z)
    mu <- a+b*Forcing
    ode{
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
covid <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])

T <- 365
nObs <- 365
# init_parameters <- list(S = 50318612.9, E=257.6, I= 4484.3, R=1873026.1, x=0)

synthetic_dataset <- bi_generate_dataset(model=covid, end_time=T,
                                         noutputs = nObs, seed="2")
synthetic_data <- bi_read(synthetic_dataset)

synthetic_df <- as.data.frame(synthetic_data)

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = y.value, colour="y.value")) +
  theme(legend.position="bottom") +
  ggtitle("Covid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = S.value, colour="S.value")) +
  geom_path(aes(y = E.value, colour="E.value")) +
  geom_path(aes(y = I.value, colour="I.value")) +
  geom_path(aes(y = R.value, colour="R.value")) +
  theme(legend.position="bottom") +
  ggtitle("Covid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")

#change seed, still be the same ODE solution of SEIR states 
#but generated y will change due to random dis generation
#Question remains: why ode solution in LibBi and R different?
#in LibBi, defaulted solver:
#’RK4(3)’ : An order 4(3) low-storage Runge-Kutta with adaptive step size h=1.

#in R deSolve, defaulted solver:
#The optional argument `method` selects the numerical integrator among several. 
#The default integrator is `lsoda`.
#LSODA is an integrator for solving stiff and non-stiff systems of ordinary differential
#equations. It was written in FORTRAN by Linda Petzold and Alan Hindmarsh. It can solve systems with dense or banded Jacobian when the problem is stiff. 
