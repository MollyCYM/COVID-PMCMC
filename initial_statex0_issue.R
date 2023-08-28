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
v <- read.csv("covidoudg2_y1w.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()

y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
L <- read.csv("Forcing.csv", header=FALSE, stringsAsFactors=FALSE)
Forcing <- data.frame(value = L) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time,V1 )
colnames(Forcing) <- c("time","value")

ncores <- 8
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
  param tau
  param a
  param b
  
  
  
  sub parameter {
    k ~ truncated_gaussian(5, 0.05, lower = 0) 
    gamma ~ truncated_gaussian(9, 0.09, lower = 0) 
    sigma ~ truncated_gaussian(sqrt(0.004), 0.001, lower = 0)
    theta ~ truncated_gaussian(0.05, 0.001, lower = 0)
    tau ~ truncated_gaussian(0.1, 0.001, lower = 0)
    a ~ gaussian(-0.02, 0.001)
    b ~ gaussian(-0.2, 0.01)
  }
  
  sub initial {
    S <- N-1
    E <- 1
    I <- 0
    R <- 0
    Z <- 1/k
    
    x ~ gaussian(-0.02,0.2) //x0 draw from unconditional distribution
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

sub proposal_parameter {
    k ~ truncated_gaussian(k, 0.01, lower = 0) 
    gamma ~ truncated_gaussian(gamma, 0.01, lower = 0) 
    sigma ~ truncated_gaussian(sigma, 0.001, lower = 0)
    theta ~ truncated_gaussian(theta, 0.001, lower = 0)
    tau ~ gaussian(tau, 0.001)
    a ~ gaussian(a, 0.001)
    b ~ gaussian(b, 0.001)
  }


}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
rewrite(bi_model)
input_lst <- list(N = 52196381,Forcing=Forcing)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(k=5, gamma=9, sigma=sqrt(0.004),theta=0.05,tau=0.1,a=-0.02,b=-0.2)

bi <- sample(bi_model,target = "posterior", end_time = end_time, init=init_list, input = input_lst, obs = obs_lst, nsamples = 2000, nparticles = minParticles, nthreads = ncores, proposal = 'model',seed=0066661) %>% 
  adapt_particles(min = minParticles, max = minParticles*500) %>%
  adapt_proposal(min = 0.1, max = 0.4) %>%
  sample(nsamples = 100, thin = 1)

bi_lst <- bi_read(bi %>% sample_obs)

###########################################PZ Model  example
rm(list=ls())
# library(tidyverse)
# library(ggplot2)
# library(ggpubr)
# library(pander)
# library(lubridate)
# library(latex2exp)
library(rbi)
library(rbi.helpers)
# Load the data
v <- read.csv("obs_P.csv", header=FALSE, stringsAsFactors=FALSE)
P_obs <- data.frame(value = v) |>
  mutate(time = seq(1, by = 1, length.out = n())) |>
  dplyr::select(time, V1)
colnames(P_obs) <- c("time", "value")

ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model PZ {
  const c = 0.25   // zooplankton clearance rate
  const e = 0.3    // zooplankton growth efficiency
  const m_l = 0.1  // zooplankton linear mortality
  const m_q = 0.1  // zooplankton quadratic mortality

  param mu, sigma  // mean and standard deviation of phytoplankton growth
  state P, Z       // phytoplankton, zooplankton
  noise alpha      // stochastic phytoplankton growth rate
  obs P_obs        // observations of phytoplankton
  
  sub parameter {
    mu ~ uniform(0.0, 1.0)
    sigma ~ uniform(0.0, 0.5)
  }
  
  sub proposal_parameter {
    mu ~ truncated_gaussian(mu, 0.02, 0.0, 1.0);
    sigma ~ truncated_gaussian(sigma, 0.01, 0.0, 0.5);
  }

  sub initial {
    P <- log(2.0)
    Z ~ log_normal(log(2.0), 0.1)
  }

  sub transition(delta = 1.0) {
    alpha ~ gaussian(mu, sigma)
    ode {
      dP/dt = alpha*P - c*P*Z
      dZ/dt = e*c*P*Z - m_l*Z - m_q*Z*Z
    }
  }

  sub observation {
    P_obs ~ log_normal(log(P), 0.2)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
# end_time <- max(P_obs$time)
end_time <- 101
# obs_lst <- list(P_obs = P_obs %>% dplyr::filter(time <= end_time))
obs_lst <- list(P_obs = P_obs %>% dplyr::filter(time <= 100))

# sim <- bi_generate_dataset(bi_model, end_time = 100)
# obs_lst <- bi_read(sim)["P_obs"]
bi <- sample(bi_model, end_time = end_time, obs = obs_lst, nsamples = 200, nparticles = minParticles, nthreads = ncores, proposal = 'model') |>
  adapt_particles(min = minParticles, max = minParticles*500) |>
  adapt_proposal(min = 0.1, max = 0.4) |>
  sample(nsamples = 100, thin = 1)

bi_lst <- bi_read(bi %>% sample_obs)

####################################################Seb example code####################################################
library(rbi)
library(rbi.helpers)

set.seed(1234)
end_time <- 101

ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model PZ {
  const c = 0.25   // zooplankton clearance rate
  const e = 0.3    // zooplankton growth efficiency
  const m_l = 0.1  // zooplankton linear mortality
  const m_q = 0.1  // zooplankton quadratic mortality

  param mu, sigma  // mean and standard deviation of phytoplankton growth
  state P, Z       // phytoplankton, zooplankton
  noise alpha      // stochastic phytoplankton growth rate
  obs P_obs        // observations of phytoplankton

  sub parameter {
    mu ~ uniform(0.0, 1.0)
    sigma ~ uniform(0.0, 0.5)
  }

  sub proposal_parameter {
    mu ~ truncated_gaussian(mu, 0.02, 0.0, 1.0);
    sigma ~ truncated_gaussian(sigma, 0.01, 0.0, 0.5);
  }

  sub initial {
    P <- log(2.0)  // Change P_0 starts from a fixed number
    Z ~ log_normal(log(2.0), 0.1)
  }

  sub transition(delta = 1.0) {
    alpha ~ gaussian(mu, sigma)
    ode {
      dP/dt = alpha*P - c*P*Z
      dZ/dt = e*c*P*Z - m_l*Z - m_q*Z*Z
    }
  }

  sub observation {
    P_obs ~ log_normal(log(P), 0.2)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
# create simulated data set
sim <- bi_generate_dataset(bi_model, end_time = 100)

obs_lst <- bi_read(sim)["P_obs"]

bi <- sample(
  bi_model, end_time = end_time, obs = obs_lst, nsamples = 200,
  nparticles = minParticles, nthreads = ncores, proposal = 'model'
) |>
  adapt_particles(min = minParticles, max = minParticles * 500) |>
  adapt_proposal(min = 0.1, max = 0.4) |>
  sample(nsamples = 100, thin = 1)




