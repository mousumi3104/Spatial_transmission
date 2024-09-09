library(cmdstanr)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(stringr)
library(abind)
library(scales)
library(bayesplot)
library(ggplot2)
library(ISOweek)
library(this.path)
library(ISOweek)
library(here)
# 
script_directory <- this.path::this.dir()
setwd(script_directory)
#-------------- loading data -------------------------------------------------------------------------

source("data/stan_data_arrangements.R")
stan_data <- stan_data_arrangements(death_threshold = 10, script_directory)

m <- cmdstan_model("fitting_regions.stan")

fit <- m$sample(
  data=stan_data,
  iter_sampling = 10, 
  iter_warmup =100, 
  parallel_chains = 4,
  chains=4,
  thin=1, 
  seed=1234,
  refresh = 20,
  adapt_delta = 0.9, 
  max_treedepth = 10)     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
                                                              # opposite for tigher adapt_delta)
                                                              # default adapt_delta=0.80, max.treedepth = 10                                                              # max_treedepth is for time efficiency concern (higher -> more time)  
out <- fit$draws(format = "matrix")
summary_fit <- fit$summary()
save(fit,stan_data,file=paste0('region_fitting.Rdata'))

bayesplot::mcmc_trace(fit$draws(c("mu[1]","mu[2]","mu[3]")))
bayesplot::mcmc_scatter(fit$draws(c("mu[1]","mu[2]")))

death_data_length <- stan_data$death_data_length
final_time <- stan_data$final_time

weekly_death <- apply(fit$draws("weekly_deaths",format="matrix"),2,mean)
plot(weekly_death[1:death_data_length],type="l")
points(stan_data$death[,1],col="red")

Rt <- apply(fit$draws("Rt",format ="matrix"),2,mean)
plot(Rt[1:final_time])
points(Rt[(final_time+1):(2*final_time)],col="red")
plot(Rt[((7*final_time)+1):(8*final_time)],col="red4")




