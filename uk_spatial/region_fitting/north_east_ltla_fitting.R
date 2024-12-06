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
library(readxl)
# 
script_directory <- this.path::this.dir()
setwd(script_directory)
#-------------- loading data -------------------------------------------------------------------------

source("data/stan_data_arrangements_daily.R")
stan_data_connected <- stan_data_arrangements(death_threshold = 10, script_directory)
M_regions <- stan_data_connected$M_regions
m <- cmdstan_model("fitting_region_daily.stan")

fit_connected <- m$sample(
  data=stan_data_connected,
  iter_sampling =200,
  iter_warmup =800,
  parallel_chains = 4,
  chains=4, 
  thin=1,
  seed=1234,
  refresh = 40,
  adapt_delta = 0.9,
  max_treedepth = 12,
  init = \() list(mu = rep(3.28,M_regions), 
                  initial_seeding = rep(10,M_regions),
                  tau = 0.01,
                  phi = 20,
                  gamma = 0.5))     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
# opposite for tigher adapt_delta)
# default adapt_delta=0.80, max.treedepth = 10                                                              # max_treedepth is for time efficiency concern (higher -> more time)
# out <- fit_connected$draws(format = "matrix")
summary_fit_connected <- fit_connected$summary()
# save(fit_connected,stan_data_connected,file=paste0('region_connceted_rt.Rdata'))
save(fit_connected,stan_data_connected,file=paste0('results/region_connected_rt.Rdata'))

# #-------- disconnected_rt ---------------------------------------------------------------------
stan_data_disconnected <- stan_data_connected
stan_data_disconnected$C_base = stan_data_connected$C_lockdown
# 
# # Sys.setenv(STAN_NUM_THREADS = 2)
# m <- cmdstan_model("fitting_regions.stan")#,cpp_options = list(stan_threads = TRUE))
# 
fit_disconnected <- m$sample(
  data=stan_data_disconnected,
  iter_sampling = 50,
  iter_warmup =300,
  parallel_chains = 4,
  # threads_per_chain = 2,
  chains=4,
  thin=1,
  seed=1234,
  refresh = 40,
  adapt_delta = 0.9,
  max_treedepth = 12,init = \() list(mu = rep(3.28,M_regions), 
                                     initial_seeding = rep(5,M_regions),
                                     tau = 0.01,
                                     phi = 20,
                                     gamma = 0.5))

# # out <-  fit_disconnected$draws(format = "matrix")
summary_fit_disconnected <- fit_disconnected$summary()
# # save(fit_disconnected,stan_data_disconnected,file=paste0('region_disconnceted_rt.Rdata'))
save(fit_disconnected,stan_data_disconnected,file=paste0('results/region_disconnected_rt.Rdata'))
# 
# source("plot_region_fitting.R")
# 
# 
# bayesplot::mcmc_trace(fit$draws(c("mu[1]","mu[2]","mu[3]")))
# # bayesplot::mcmc_scatter(fit$draws(c("mu[1]","mu[2]")))
# # 
# death_data_length <- stan_data_connected$death_data_length
# #final_time <- stan_data$final_time
# 
# weekly_death <- apply(fit$draws("weekly_deaths",format="matrix"),2,mean)
# plot(weekly_death[((2*death_data_length)+1) : (3*death_data_length)],type="l")
# points(stan_data_connected$death[,3],col="red")

