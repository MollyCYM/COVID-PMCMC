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
model Lorenz {
  const rho = 45.92
  const beta = 4.0
  const alpha = 16.0
  
  
  state X, Y, Z
  obs X_obs
  
  
  sub initial {
    X ~ log_normal(log(1.0), 0.00002)
    Y ~ log_normal(log(1.0), 0.00002)
    Z ~ log_normal(log(1.0), 0.00002)
  }
  
  
  sub transition(delta = 0.0001) {
    ode {
      dX/dt = alpha * (Y - X)
      dY/dt = X * (rho - Z) - Y
      dZ/dt = X * Y - beta * Z
    }
  }
  
  
  sub observation {
    X_obs ~ normal(X, 0.2)
  }
}
"

Lorenz <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])

T <- 0.2
nObs <- 100
init_parameters <- list(X = 1, Y = 1, Z = 1)


synthetic_dataset <- bi_generate_dataset(model=Lorenz, end_time=T,
                                         init=init_parameters,
                                         noutputs = nObs)

synthetic_data <- bi_read(synthetic_dataset)


synthetic_df <- as.data.frame(synthetic_data)


ggplot(synthetic_df, aes(X.time)) +
  geom_path(aes(y = X.value, colour="X.value")) +
  geom_path(aes(y = X_obs.value, colour="X_obs.value")) +
  theme(legend.position="bottom") +
  ggtitle("Lorenz") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")

#####################################################################################
rm(list=ls())
library('rbi')
library(ggplot2)


model_file_name <- "Lorenz.bi"


Lorenz <- bi_model(model_file_name)
##
synthetic_dataset <- bi_generate_dataset(model=Lorenz,end_time=0.2,noutputs=100)
synthetic_data <- bi_read(synthetic_dataset)
synthetic_df <- as.data.frame(synthetic_data)
#Or
T <- 0.2
nObs <- 100
init_parameters <- list(X = 1, Y = 1, Z = 1)


synthetic_dataset <- bi_generate_dataset(model=Lorenz,end_time=T,
                                         init=init_parameters,
                                         noutputs = nObs)
synthetic_data <- bi_read(synthetic_dataset)
synthetic_df <- as.data.frame(synthetic_data)
##

ggplot(synthetic_df, aes(X.time)) +
  geom_path(aes(y = X.value, colour="X.value")) +
  geom_path(aes(y = X_obs.value, colour="X_obs.value")) +
  theme(legend.position="bottom") +
  ggtitle("Lorenz") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Time") +
  ylab("Value")