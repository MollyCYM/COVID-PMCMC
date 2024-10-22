rm(list=ls())
set.seed(0000889977)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
ncores <- 12
# Load the data
v <- read.csv("svm_y2r.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time","value")
model_str <- "
model StochasticVolatility {

  obs y
  
  state alpha
  
  noise epsilon
  
  param beta
  param phi
  param sigma_eta
  
  sub parameter {
    beta ~ truncated_gaussian(0.8, 0.1, lower = 0)
    phi ~ truncated_gaussian(0.45, 1, lower = -1, upper=1)
    sigma_eta ~ truncated_gaussian(sqrt(0.04), 0.1, lower=0)
  }
  
  sub proposal_parameter {
    beta ~ truncated_gaussian(beta, 0.1, lower = 0) 
    phi ~ truncated_gaussian(phi, 0.01, lower = -1, upper=1) 
    sigma_eta ~ truncated_gaussian(sigma_eta, 0.001)
  }
  
  sub initial {
    /* stationary distribution for alpha */
    alpha ~ normal(0.0, sqrt(sigma_eta**2/(1.0 - phi**2)));
  }

  sub transition {
    epsilon ~ normal(0,1);
    alpha <- phi*alpha + sigma_eta*epsilon;
  }

  sub observation {
    y ~ normal(0, beta*exp(0.5*alpha));
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(beta=0.8,phi=0.45,sigma_eta=sqrt(0.04))
#LibBi wrapper 
#run launches LibBi with a particular set of command line arguments
bi_model <- libbi(model,end_time = end_time, 
                  init=init_list, obs = obs_lst)
#RBi.helpers adapt_particle
particles_adapted <- bi_model %>%
  sample(nsamples = 2000, nparticles = 16, 
         nthreads = ncores, proposal = 'prior') %>%
  adapt_particles(min = 16, max = 16*500)

#RBi.helpers adapt_proposal
proposal_adapted <- particles_adapted %>%
  sample(target = "posterior", nsamples = 2000, 
         nthreads = ncores, proposal = 'model') %>%
  adapt_proposal(min = 0.1, max = 0.4)

#Running pMCMC with burn-in
bi <- proposal_adapted %>%
  sample(nsamples = 5000, thin = 1,init=init_list) %>%
  sample(nsamples = 15000, thin = 1)
bi_lst <- bi_read(bi %>% sample_obs)

write.csv(bi_lst,"../data/svminf_model2.csv")
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
write.csv(fitY,"../data/svminf_y2.csv")

fitalpha <- bi_lst$alpha %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()
write.csv(fitalpha,"../data/svminf_alpha2.csv")

write.csv(bi_lst$beta$value,"../data/svminf_beta2.csv")
write.csv(bi_lst$phi$value,"../data/svminf_phi2.csv")
write.csv(bi_lst$sigma_eta$value,"../data/svminf_sigma2.csv")