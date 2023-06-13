rm(list=ls())
library(rbi)
library(rbi.helpers)
library(stringi) ## for reading the model from a string
library(pomp) ## for the bsflu data set
library(tidyverse)
set.seed(296825852)
head(bsflu)
model_str <- '
model bsflu {
 const N = 763
 const timestep = 1/12
 param Beta, mu_I, mu_R1,  rho
 state S, I, R1
 noise infection, recovery, leave_bed
 obs Incidence
sub parameter {
 Beta ~ uniform(1, 5)
 mu_I ~ uniform(0.5, 3)
 rho ~ uniform(0.5, 1)
}
sub initial {
 S <- N - 1 // susceptibles
 I <- 1 // infectious
 R1 <- 1 // recovered but bed-confined
}
sub transition (delta = timestep) {
 infection ~ binomial(S, 1 - exp(-Beta * I/N * timestep))
 recovery ~ binomial(I, 1 - exp(-mu_I * timestep))
 leave_bed ~ binomial(R1, 1 - exp(-mu_R1 * timestep))
 S <- S - infection
 I <- I + infection - recovery
 R1 <- R1 + recovery - leave_bed
}
sub observation {
 Incidence ~ poisson(rho * R1 + 1e-6)
}
}
'
flu_model <- bi_model(lines = stri_split_lines(model_str)[[1]]) %>%
  fix(mu_R1 = 1/(sum(bsflu$B)/512))
obs <- bsflu %>%
  select(time=date, value=B) %>%
  list(Incidence=.) %>%
  time_to_numeric(origin=as.Date("1978-01-21"), unit="day")
bi <- libbi(model=flu_model, obs=obs, end_time=nrow(bsflu))
rewrite(bi)
flu_model

params <- c(Beta=2,mu_I=1,rho=0.9,mu_R1=1/3,mu_R2=1/2)
sim <- rbi::simulate(bi, init=as.list(params), nsamples=10)
sim_res <- bi_read(sim, type="state")
ggplot(sim_res$R1, aes(x=time, group=np))+
  geom_line(aes(y=value)) +
  ylab("R1")

particles_adapted <- bi %>%
  sample(proposal="prior", nsamples=2000, nparticles=1024) %>%
  adapt_particles(max=2**20)
nparticles <- particles_adapted$options$nparticles
nparticles

proposal_adapted <- particles_adapted %>%
  sample(proposal="prior", nsamples=2000) %>%
  adapt_proposal(min=0.1, max=0.3, adapt="both")

posterior <- proposal_adapted %>%
  sample(nsamples=10000)

#########################################################
rm(list=ls())
library(rbi)
library(rbi.helpers)
library(stringi) ## for reading the model from a string
library(pomp) ## for the bsflu data set
library(tidyverse)
set.seed(296825852)
head(bsflu)
mu_R1=1/(sum(bsflu$B)/512)
model_str <- '
model bsflu {
 const N = 763
 const timestep = 1/12
 input mu_R1
 param Beta, mu_I, rho
 state S, I, R1
 noise infection, recovery, leave_bed
 obs Incidence
sub parameter {
 Beta ~ uniform(1, 5)
 mu_I ~ uniform(0.5, 3)
 rho ~ uniform(0.5, 1)
}
sub initial {
 S <- N - 1 // susceptibles
 I <- 1 // infectious
 R1 <- 1 // recovered but bed-confined
}
sub transition (delta = timestep) {
 infection ~ binomial(S, 1 - exp(-Beta * I/N * timestep))
 recovery ~ binomial(I, 1 - exp(-mu_I * timestep))
 leave_bed ~ binomial(R1, 1 - exp(-mu_R1 * timestep))
 S <- S - infection
 I <- I + infection - recovery
 R1 <- R1 + recovery - leave_bed
}
sub observation {
 Incidence ~ poisson(rho * R1 + 1e-6)
}
}
'
flu_model <- bi_model(lines = stri_split_lines(model_str)[[1]])
obs <- bsflu %>%
  select(time=date, value=B) %>%
  list(Incidence=.) %>%
  time_to_numeric(origin=as.Date("1978-01-21"), unit="day")
input_lst <- list(mu_R1)
bi <- libbi(model=flu_model, obs=obs, end_time=nrow(bsflu))
bi <- attach_data(bi, "input", input_lst)
rewrite(bi)
flu_model

params <- c(Beta=2,mu_I=1,rho=0.9,mu_R1=1/3,mu_R2=1/2)
sim <- rbi::simulate(bi, init=as.list(params), nsamples=10)
sim_res <- bi_read(sim, type="state")
ggplot(sim_res$R1, aes(x=time, group=np))+
  geom_line(aes(y=value)) +
  ylab("R1")
ggplot(sim_res$S, aes(x=time, group=np))+
  geom_line(aes(y=value)) +
  ylab("S")
ggplot(sim_res$I, aes(x=time, group=np))+
  geom_line(aes(y=value)) +
  ylab("I")
