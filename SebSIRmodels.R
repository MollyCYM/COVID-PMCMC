rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)

model_file <- system.file(package="rbi", "SIR.bi")
SIRmodel <- bi_model(model_file) # load model
SIRmodel

#Generating a dataset
set.seed(1001912)
SIRdata <- bi_generate_dataset(SIRmodel, end_time=16*7, noutputs=16)
#The generated dataset can be viewed and/or stored in a variable using bi_read:
dataset <- bi_read(SIRdata)
#We can visualise the generated incidence data with
plot(dataset$Incidence$time, dataset$Incidence$value)
lines(dataset$Incidence$time, dataset$Incidence$value)

class(SIRdata) #identify SIRdata is a libbi object
#The standard way of creating a libbi object for Bayesian inference is using the libbi command
#Th bi_generate_dataset is one particular way of generating a libbi object, 
#used only to generate test data from a model.
bi <- libbi(SIRmodel)
class(bi)
#Sample from the prior of the SIR model:
bi_prior <- sample(bi, target="prior", nsamples=1000, end_time=16*7, noutputs=16)
bi_prior
str(bi_prior)
bi_prior$options
bi_prior$output_file_name
#We can get the results of the sampling run using bi_read
prior <- bi_read(bi_prior$output_file_name)
#Or shortversion
prior <- bi_read(bi_prior)
str(prior)


#Fitting a model to data using PMCMC
#The following command will generate 16 * 10,000 = 160,000 simulations 
#and therefore may take a little while to run (if you want to see the samples progress, 
#use verbose=TRUE in the sample call).
bi <- sample(bi_prior, target="posterior", nparticles=32, obs=SIRdata)
#This samples from the posterior distribution. 
#Remember that options are preserved from previous runs (because we passed the bi) as first argument, 
#so we don't need to specify nsamples, end_time and noutputs again, unless we want to change them. 
#The nparticles option specifies the number of particles

#You can also pass a list of data frames (each element of the list corresponding to 
#one observed variable as the obs argument, for example
df <- data.frame(time = c(0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112),
                 value = c(1,6,2,26,99,57,78,57,15,9,4,1,1,1,0,2,0))
bi_df <- sample(bi_prior, target="posterior", nparticles=32, obs=list(Incidence=df))
#Is this same as?
bi_df <- sample(bi, target="posterior", nsamples=1000, end_time=16*7, noutputs=16,
                nparticles=32, obs=list(Incidence=df))
###Input, init and observation files (see the LibBi manual for details) can be specified using the init, 
###input, obs options, respectively. They can each be specified either as the name of a NetCDF file 
###containing the data, or a libbi object (in which case the output file will be taken) or directly via
###an appropriate R object containing the data (e.g., a character vector of length one, 
###or a list of data frames or numeric vectors). In the case of the command above, init is specified as 
###a list, and obs as a libbi object. The Incidence variable of the SIRdata object will be taken as 
###observations.

###The time dimension (or column, if a data frame) in the passed init, input and/or obs files 
###can be specified using the time_dim option. If this is not given, it will be assumed to be time, 
###if such a dimension exists or, if not, any numeric column not called value (or the contents of 
###the value_column option). If this does not produce a unique column name, an error will be thrown. 
###All other dimensions/columns in the passed options will be interpreted as additional dimensions 
###in the data, and stored in the dims field of the libbi object.

#Analysing an MCMC run
#Get the results of the preceding sample command:
bi_contents(bi)
posterior <- bi_read(bi) #Compared with previous data generated bi_read:
#We can see that this has two more objects than previously when we specified target="prior": loglikelihood 
#(the estimated log-likelihood of the parameters at each MCMC step) and 
#logprior (the estimated log-prior density of the parameters at each MCMC step).
#Get a summary of the parameters sampled
summary(bi)
#A summary of sampled trajectories can be obtained using
summary(bi, type="state")
#Any particular posterior sample can be viewed with extract_sample 
#(with indices running from 0 to (nsamples-1)):
extract_sample(bi, 314)

#To analyse MCMC outputs, we can use the coda package and the get_traces function of RBi. 
#Note that, to get exactly the same traces, you would have to set the seed as above.
library('coda')
traces <- mcmc(get_traces(bi))

#visualise parameter traces and densities with
plot(traces)
bi_read(SIRdata, type="param") #Compared to true generated parameter values
#Plot & Inference of MCMC part: For more details on using coda to further analyse the chains, 
#see the website of the coda package. For more plotting functionality, the ggmcmc package 
#is also worth considering.

#Predictions
#We can use the predict function to re-simulate the fitted model using the estimated parameters
pred_bi <- predict(bi, start_time=0, end_time=20*7, output_every=7, with=c("transform-obs-to-state"))
#with=c("transform-obs-to-state") tells LibBi to treat observations as a state variable, 
#that is to randomly generate observations
ps <- summary(pred_bi, type="obs")
library('ggplot2')
ggplot(ps, aes(x=time))+
  geom_line(aes(y=Median)) +
  geom_ribbon(aes(ymin=`1st Qu.`, ymax=`3rd Qu.`), alpha=0.5) +
  geom_point(aes(y=value), dataset$Incidence, color="darkred") +
  ylab("cases")
#Sample observations
obs_bi <- sample_obs(bi)
summary(obs_bi, type="obs")
dataset$Incidence
os <- summary(obs_bi, type="obs")
ggplot(os, aes(x=time))+
  geom_line(aes(y=Median)) +
  geom_ribbon(aes(ymin=`1st Qu.`, ymax=`3rd Qu.`), alpha=0.5) +
  geom_point(aes(y=value), dataset$Incidence, color="darkred") +
  ylab("cases")