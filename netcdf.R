library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
nc_data <- nc_open('input.nc')
# Save the print(nc) dump to a text file
{
  sink('input.txt')
  print(nc_data)
  sink()
}
variable_F <- ncvar_get(nc_data,"F")
variable_time <- ncvar_get(nc_data,"time")
plo
nlon <- dim(lon)
variable=nc_data[15]
print(typeof(variable))
print(variable)
variable_time=variable[var]

# set path and filename
ncpath <- "/Users/mollycui/Desktop/R script/Research 10-Epid PMCMC/"
ncname <- "input"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
#dname <- "tmp"  # note: tmp means temperature (not temporary)
# open a netCDF file
ncin <- nc_open(ncfname)
print(ncin)
# get longitude and latitude
var_F <- ncvar_get(ncin,"F")
nF <- dim(var_F)
head(var_F)

var_time <- ncvar_get(ncin,"time")
ntime <- dim(var_time)
head(var_time)


#############################bi_open function
#' @rdname bi_open
#' @name bi_open
#' @title Bi open
#' @description
#' This function opens an NetCDF file
#' The file can be specified as a string to the filepath, in which
#' case a NetCDF connection is opened, or directly as a NetCDF connection.
#'
#' @param x either a path to a NetCDF file, or a NetCDF connection created using \code{nc_open}, or a \code{\link{libbi}} object from which to read the output
#' @param file file to open (out of "input", "init", "obs", "output"), if \code{x} is given as a \code{libbi} object; by default, will read output file
#' @return an open NetCDF connection
#' @importFrom ncdf4 nc_open
bi_open <- function(x, file = "output")
{
  if (!missing(file) && class(x) != "libbi") {
    warning("'file' given to 'bi_open' although 'x' is not a 'libbi' object; will be ignored")
  }
  
  if (typeof(x) == "character"){
    nc <- nc_open(tools::file_path_as_absolute(x))
  } else if (class(x) == "ncdf4") {
    nc <- x
  } else if (class(x) == "libbi"){
    if (!(missing(file) || file == "output")) {
      opt_name <- paste(file, "file", sep="-")
      if (!(opt_name %in% names(x$options))) {
        stop("libbi object does not contain an '", file, "' file")
      }
      filename <- x$options[[opt_name]]
    } else {
      assert_files(x)
      filename <- x$output_file_name
    }
    nc <- nc_open(filename)
  } else {
    stop("'x' must be a 'character', 'ncdf4' or 'libbi' object.")
  }
  
  return(nc)
}
###############################################################################
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
filename <- tempfile(pattern="dummy", fileext=".nc")
a <- 3
b <- data.frame(dim_a = rep(1:3, time = 2), dim_b = rep(1:2, each = 3), value = 1:6)
variables <- list(a=a, b=b)
bi_write(filename, variables)
bi_file_summary(filename)

File <- read.csv("F.csv", header=TRUE, stringsAsFactors=FALSE)
filename <- tempfile(pattern="file", fileext=".nc")
time <-data.frame(value=File$time)
value <-data.frame(value=File$F)
variables <-list(time=time, value=value)
bi_write(filename,variables,dim_factors=259)
bi_file_summary(filename)

library(tidyr)
# Example df with same names
df <- data.frame(long = rep(1:10,2), lat = 1:40, elev = 1:40)
# Creates a 2D array
df_new <- tidyr::pivot_wider(df, names_from = long, values_from = elev)
df_new <- df_new[,-1]
#######################################################################
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
#bi_open("/Users/mollycui/Desktop/R script/Research 10-Epid PMCMC/input.nc",file="input")
# Load the data
v <- read.csv("covid259days3.csv", header=FALSE, stringsAsFactors=FALSE) 
y <- data.frame(value = v) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, V1)
colnames(y) <- c("time", "value")
var_F <- bi_open("/Users/mollycui/Desktop/R script/Research 10-Epid PMCMC/input.nc",file="input")
Forcing <- ncvar_get(var_F,"F")
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
      dmu/dt = exp(a + b*Forcing)
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
input_lst <- list(Forcing = Forcing)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 1000, thin = 1) %>% # burn in 
  sample(nsamples = 10000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)
##################################################################