#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#} else if (length(args)==1) {
#  # default output file
#  args[2] = "out.txt"
#}


library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pander)
library(lubridate)
library(latex2exp)
library(rbi)
library(rbi.helpers)
library(readr)
options(digits=2)
# Load the data
daily_data <-as.logical(args[1])
print(daily_data)
name_data <-args[2]
v <- read.csv(name_data, header=FALSE, stringsAsFactors=FALSE)

if (daily_data == TRUE) {
  print("here")
  y <- data.frame(value = v) %>%
    mutate(time = seq(1, by = 1, length.out = n())) %>%
    dplyr::select(time, V1)
  colnames(y) <- c("time", "value")
  print(y)
} else {
  y <- data.frame(value = v) %>%
    mutate(time = seq(7, by = 7, length.out = n())) %>%
    dplyr::select(time, value)
}
ncores <- as.integer(args[3])
print(ncores)
minParticles <- max(ncores, 16)
model_str <- toString(read.table(file = args[4], header = FALSE))
print(model_str)
print(typeof(model_str))
print(typeof(toString(model_str)))

model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
print(bi_model)


input_lst <- list(N = as.integer(args[5]))
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))
init_list <- list(sigma =as.double(args[6]), gamma =as.double(args[7]), beta=as.double(args[8]), mu =as.double(args[9]))
print(init_list)


bi <- sample(bi_model, end_time = end_time, input = input_lst, init=init_list, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'model') %>%
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = as.integer(args[10]), thin = 1) %>% # burn in
  sample(nsamples = as.integer(args[11]), thin = 5)





bi_lst <- bi_read(bi %>% sample_obs)
