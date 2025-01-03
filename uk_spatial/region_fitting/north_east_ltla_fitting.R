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

source("data/data_arrangements_ne.R")
stan_data_connected <- stan_data_arrangements(death_threshold = 4, script_directory)

plot_required_date <- list(inf_start_date = stan_data_connected$inf_start_date,
                           fitting_start_date = stan_data_connected$fitting_start_date,
                           end_date = stan_data_connected$end_date)

stan_data_connected$inf_start_date <- NULL
stan_data_connected$fitting_start_date <- NULL
stan_data_connected$end_date <- NULL

M_regions <- stan_data_connected$M_regions
m <- cmdstan_model("ltla_fitting_ne.stan")

fit_connected <- m$sample(
  data=stan_data_connected,
  iter_sampling =500,
  iter_warmup =1200,
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

summary_fit_connected <- fit_connected$summary()
out_rt <- fit_connected$draws("Rt",format = "matrix")
Rt <- apply(out_rt,2,function(col) mean(col))
out_inf <- fit_connected$draws("infection",format = "matrix")
inf <- apply(out_inf,2,function(col) mean(col))
out_death <- fit_connected$draws("weekly_deaths",format = "matrix")
death <- apply(out_inf,2,function(col) mean(col))

# save(fit_connected,stan_data_connected,file=paste0('region_connceted_rt.Rdata'))
# save(fit_connected,stan_data_connected,plot_required_date, file=paste0('results/ltla_connected_ne.Rdata'))
save(fit_connected,stan_data_connected,plot_required_date, file=paste0('results/ltla_connected_ne_jan.Rdata'))

# #-------- disconnected_rt ---------------------------------------------------------------------
stan_data_disconnected <- stan_data_connected
stan_data_disconnected$C_base = stan_data_connected$C_lockdown
# 
# # Sys.setenv(STAN_NUM_THREADS = 2)
# m <- cmdstan_model("fitting_regions.stan")#,cpp_options = list(stan_threads = TRUE))
# 
fit_disconnected <- m$sample(
  data=stan_data_disconnected,
  iter_sampling = 500,
  iter_warmup =1200,
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
# save(fit_disconnected,stan_data_disconnected,plot_required_date,file=paste0('results/ltla_disconnected_ne.Rdata'))
save(fit_disconnected,stan_data_disconnected,plot_required_date,file=paste0('results/ltla_disconnected_ne_jan.Rdata'))

source("plot_ltla_ne.R")
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

