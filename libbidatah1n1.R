rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)


model_str <- "
model h1n1bm {
  
  obs y
  const k = 1.59
  const gamma = 1.08
  const sigma = 0.07
  const tau = 0.1
  const N = 52196381
  const h = 1
  
  param E0
  param I0
  param R0
  param x0
  param eta
  
  state S
  state E
  state I
  state R
  state x
  state Z
  state e
  
  sub parameter {
    x0 ~ uniform(-5,2)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
  }

  sub initial {
    S <- N
    R <- R0*S
    S <- S - R

    E <- exp(E0 + log(S))
    S <- S - E
    I <- exp(I0 + log(S))
    S <- S - I
    //x <- x0
    Z <- 0
    e <- 0
  }

  sub transition(delta = h) {
    
    eta ~ gaussian(0,0.1*sqrt(h))
    e <- e+eta
    ode {
      dx/dt = sigma*eta
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

}"

h1n1 <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
rewrite(h1n1)
T <- 365
nObs <- 365
synthetic_dataset <- bi_generate_dataset(model=h1n1, end_time=T,
                                         noutputs = nObs, seed="53")
##Good seeds: 06 08 13 15 24-25 28 32 36-37 46
synthetic_data <- bi_read(synthetic_dataset)
synthetic_df <- as.data.frame(synthetic_data)
# write.csv(synthetic_df$y.value,"libbih1n1_Y1.csv")
# write.csv(synthetic_df$x.value,"libbih1n1_x1.csv")
# write.csv(synthetic_df,"libbih1n1_model1.csv")

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = exp(x.value), colour="x.value")) +
  theme(legend.position="bottom") +
  ggtitle("H1N1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = y.value, colour="y.value")) +
  #geom_path(aes(y = Z.value, colour="Z.value")) +
  theme(legend.position="bottom") +
  ggtitle("H1N1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = S.value, colour="S.value")) +
  geom_path(aes(y = E.value, colour="E.value")) +
  geom_path(aes(y = I.value, colour="I.value")) +
  geom_path(aes(y = R.value, colour="R.value")) +
  theme(legend.position="bottom") +
  ggtitle("H1N1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")
###################################Fix initial state values#############################
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)

model_str <- "
model h1n1bm {
  
  obs y
  const k = 1.59
  const gamma = 1.08
  const sigma = 0.07
  const tau = 0.1
  const N = 52196381
  
  state S
  state E
  state I
  state R
  state x
  state Z

  sub initial {
    S <- 50318612.9
    E <- 257.6
    I <- 4484.3
    R <- 1873026.1
    x <- 0
  }

  sub transition(delta = 1) {
    noise e
    e ~ wiener()
    ode {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

}"

h1n1 <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])

T <- 365
nObs <- 365
# init_parameters <- list(S = 50318612.9, E=257.6, I= 4484.3, R=1873026.1, x=0)

synthetic_dataset <- bi_generate_dataset(model=h1n1, end_time=T,
                                         noutputs = nObs, seed="11")
#Is that same?
synthetic_dataset <- bi_generate_dataset(model=h1n1, end_time=T,
                                         init=init_parameters,
                                         output_every = 1)
##
synthetic_data <- bi_read(synthetic_dataset)


synthetic_df <- as.data.frame(synthetic_data)

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = y.value, colour="y.value")) +
  geom_path(aes(y = Z.value, colour="Z.value")) +
  theme(legend.position="bottom") +
  ggtitle("H1N1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")
ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = y.value, colour="y.value")) +
  theme(legend.position="bottom") +
  ggtitle("H1N1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")
#################################Change S0 E0 I0 R0####################################
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)

model_str <- "
model h1n1bm {
  
  obs y
  const k = 1.59
  const gamma = 1.08
  const sigma = 0.07
  const tau = 0.1
  const N = 52196381
  
  state S
  state E
  state I
  state R
  state x
  state Z

  sub initial {
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
    x <- 0
  }

  sub transition(delta = 1) {
    noise e
    e ~ wiener()
    ode {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/5, 0)), tau)
  }

}"

h1n1 <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])

T <- 365
nObs <- 365
# init_parameters <- list(S = 50318612.9, E=257.6, I= 4484.3, R=1873026.1, x=0)

synthetic_dataset <- bi_generate_dataset(model=h1n1, end_time=T,
                                         noutputs = nObs, seed="1")
synthetic_data <- bi_read(synthetic_dataset)

synthetic_df <- as.data.frame(synthetic_data)

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = y.value, colour="y.value")) +
  theme(legend.position="bottom") +
  ggtitle("H1N1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")

ggplot(synthetic_df, aes(y.time)) +
  geom_path(aes(y = y.value, colour="y.value")) +
  geom_path(aes(y = S.value, colour="S.value")) +
  geom_path(aes(y = E.value, colour="E.value")) +
  geom_path(aes(y = I.value, colour="I.value")) +
  geom_path(aes(y = R.value, colour="R.value")) +
  theme(legend.position="bottom") +
  ggtitle("H1N1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")