# estimation of Rt for regions of England (Fig 4)

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
library(readxl)
library(this.path)
# 
script_directory <- this.path::this.dir()
setwd(script_directory)
#-------------- loading data -------------------------------------------------------------------------

source("data/stan_data_arrangements.R")
stan_data_connected <- stan_data_arrangements(death_threshold = 10, script_directory)

plot_required_date <- list(inf_start_date = stan_data_connected$inf_start_date,
                           fitting_start_date = stan_data_connected$fitting_start_date,
                           end_date = stan_data_connected$end_date)

stan_data_connected$inf_start_date <- NULL
stan_data_connected$fitting_start_date <- NULL
stan_data_connected$end_date <- NULL

M_regions <- stan_data_connected$M_regions
m <- cmdstan_model("fitting_regions.stan")

fit_connected <- m$sample(
  data=stan_data_connected,
  iter_sampling = 800,
  iter_warmup =1000,
  parallel_chains = 4,
  chains=4, 
  thin=1,
  seed=1234,
  refresh = 40,
  adapt_delta = 0.97,
  max_treedepth = 13)
  # init = \() list(mu = rep(3.28,M_regions),
  #                 initial_seeding = rep(5,M_regions),
  #                 tau = 0.01,
  #                 phi = 20) )    # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
                                                              # opposite for tigher adapt_delta)
                                                              # default adapt_delta=0.80, max.treedepth = 10                                                              # max_treedepth is for time efficiency concern (higher -> more time)
# out <- fit_connected$draws(format = "matrix")
summary_fit_connected <- fit_connected$summary()
diagnostics <- fit_connected$diagnostic_summary()
print(diagnostics)
log_lik_summary <- fit_connected$draws("log_lik",format="matrix")
loo_result <- loo(log_lik_summary)
print(loo_result)
# save(fit_connected,stan_data_connected,file=paste0('region_connceted_rt.Rdata'))
save(fit_connected,stan_data_connected,plot_required_date,file=paste0('results/updated_si/gm_region_connected_rt_jan.Rdata'))

# #-------- disconnected_rt ---------------------------------------------------------------------
stan_data_disconnected <- stan_data_connected
stan_data_disconnected$C_base = diag(M_regions)
stan_data_disconnected$C_lockdown = diag(M_regions)


fit_disconnected <- m$sample(
  data=stan_data_disconnected,
  iter_sampling = 800,
  iter_warmup =1000,
  parallel_chains = 4,
  # threads_per_chain = 2,
  chains=4,
  thin=1,
  seed=1234,
  refresh = 40,
  adapt_delta = 0.97,
  max_treedepth = 13)#init = \() list(mu = rep(3.28,M_regions), 
   #                                    initial_seeding = rep(5,M_regions),
  #                                    tau = 0.01,
  #                                    phi = 20,
  #                                    gamma = 0.5))

# # out <-  fit_disconnected$draws(format = "matrix")
summary_fit_disconnected <- fit_disconnected$summary()
diagnostics <- fit_disconnected$diagnostic_summary()
print(diagnostics)
log_lik_summary <- fit_disconnected$draws("log_lik",format="matrix")
loo_result <- loo(log_lik_summary)
print(loo_result)
# # save(fit_disconnected,stan_data_disconnected,file=paste0('region_disconnceted_rt.Rdata'))
save(fit_disconnected,stan_data_disconnected,plot_required_date,file=paste0('results/updated_si/region_disconnected_rt_jan.Rdata'))
# 
source("plot_region_fitting.R")
