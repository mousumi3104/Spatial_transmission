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
  iter_sampling = 500,
  iter_warmup =1200,
  parallel_chains = 4,
  chains=4, 
  thin=1,
  seed=1234,
  refresh = 40,
  adapt_delta = 0.9,
  max_treedepth = 12,
  init = \() list(mu = rep(3.28,M_regions), 
                  initial_seeding = rep(5,M_regions),
                  tau = 0.01,
                  phi = 20,
                  gamma = 0.5))     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
                                                              # opposite for tigher adapt_delta)
                                                              # default adapt_delta=0.80, max.treedepth = 10                                                              # max_treedepth is for time efficiency concern (higher -> more time)
# out <- fit_connected$draws(format = "matrix")
summary_fit_connected <- fit_connected$summary()
# save(fit_connected,stan_data_connected,file=paste0('region_connceted_rt.Rdata'))
save(fit_connected,stan_data_connected,plot_required_date,file=paste0('results/region_connected_rt_including_jan.Rdata'))

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
save(fit_disconnected,stan_data_disconnected,plot_required_date,file=paste0('results/region_disconnected_rt_including_jan.Rdata'))
# 
source("plot_region_fitting.R")
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

# fit <-  fit_disconnected
# Rt <- apply(fit$draws("Rt",format ="matrix"),2,mean)
# plot(Rt[1:final_time],type="l",col = "black")
# lines(Rt[(final_time+1):(2*final_time)],col="black",type="l")
# plot(Rt[((2*final_time)+1):(3*final_time)],col="red",type="l")
# plot(Rt[((3*final_time)+1):(4*final_time)],col="red",type="l")
# plot(Rt[((4*final_time)+1):(5*final_time)],col="red",type="l")
# plot(Rt[((5*final_time)+1):(6*final_time)],col="red",type="l")
# plot(Rt[((6*final_time)+1):(7*final_time)],col="red",type="l")
# plot(Rt[((7*final_time)+1):(8*final_time)],col="red",type="l")
# plot(Rt[((8*final_time)+1):(9*final_time)],col="red",type="l")
# # 
# 
# 
# fit <-  fit_connected
# Rt <- apply(fit$draws("Rt",format ="matrix"),2,mean)
# plot(Rt[((8*final_time)+1):(9*final_time)],type="l",col = "red")
# load("results/estimated_separate_rt_1.Rdata")
# out <- rstan::extract(fit)
# s <- summary(fit)
# estimated_rt <- out$Rt
# Rt = colMeans(estimated_rt)
# lines(Rt,col="red")
# # [(length(Rt)-340):length(Rt)]
# abline(v=c(a,b))


