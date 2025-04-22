# # estimation of Rt for regions of England (Fig 5)


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
library(loo)
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

start_time <- Sys.time()
fit_connected <- m$sample(
  data=stan_data_connected,
  iter_sampling = 800,
  iter_warmup = 1000,
  parallel_chains = 4,
  chains=4, 
  thin=1,
  seed=1234,
  refresh = 40,
  adapt_delta = 0.97,
  max_treedepth = 13)
  #init = \() list(mu = rep(3,M_regions), 
                  # initial_seeding = rep(2,M_regions),
                  # tau = 2,
                  # phi = 20,
                  # gamma = 0.5))     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
# opposite for higher adapt_delta)
# default adapt_delta=0.80, max.treedepth = 10                                                              # max_treedepth is for time efficiency concern (higher -> more time)
end_time <- Sys.time()
total_time <- end_time-start_time
summary_fit_connected <- fit_connected$summary()
diagnostics <- fit_connected$diagnostic_summary()
print(diagnostics)

# calculate ELPD score
log_lik_summary <- fit_connected$draws("log_lik",format="matrix")
# log_lik <- matrix(log_lik_summary,ncol=624)
# log_lik1 <- log_lik[, -c(sapply(0:11, function(x) (52*x + 1):(52*x + 6)))]
loo_result <- loo(log_lik_summary)
print(loo_result)

out_rt <- fit_connected$draws("Rt",format = "matrix")
Rt <- apply(out_rt,2,function(col) mean(col))
out_inf <- fit_connected$draws("infection",format = "matrix")
inf <- apply(out_inf,2,function(col) mean(col))
out_death <- fit_connected$draws("weekly_deaths",format = "matrix")
death <- apply(out_inf,2,function(col) mean(col))

# save(fit_connected,stan_data_connected,file=paste0('region_connceted_rt.Rdata'))
# save(fit_connected,stan_data_connected,plot_required_date, file=paste0('results/ltla_connected_ne.Rdata'))
save(fit_connected,stan_data_connected,plot_required_date, file=paste0('results/updated_si/gm_ltla_connected_ne_jan.Rdata'))

# #-------- disconnected_rt ---------------------------------------------------------------------
stan_data_disconnected <- stan_data_connected
stan_data_disconnected$C_base = diag(M_regions)
stan_data_disconnected$C_lockdown = diag(M_regions)
#stan_data_connected$C_lockdown

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
  max_treedepth = 13)#,init = \() list(mu = rep(3.28,M_regions), 
                                     # initial_seeding = rep(5,M_regions),
                                     # tau = 0.01,
                                     # phi = 20,
                                     # gamma = 0.5))

   # # out <-  fit_disconnected$draws(format = "matrix")
summary_fit_disconnected <- fit_disconnected$summary()
diagnostics <- fit_disconnected$diagnostic_summary()
print(diagnostics)
log_lik_summary <- fit_disconnected$draws("log_lik",format="matrix")
loo_result <- loo(log_lik_summary)
print(loo_result)
# # save(fit_disconnected,stan_data_disconnected,file=paste0('region_disconnceted_rt.Rdata'))
# save(fit_disconnected,stan_data_disconnected,plot_required_date,file=paste0('results/ltla_disconnected_ne.Rdata'))
save(fit_disconnected,stan_data_disconnected,plot_required_date,file=paste0('results/updated_si/gm_ltla_disconnected_ne_jan.Rdata'))

source("plot_ltla_ne.R")
