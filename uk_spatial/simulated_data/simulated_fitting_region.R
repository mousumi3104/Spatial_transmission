# estimation of Rt with the simulated data, connected and disconnected (Fig 3 (a-c)) 

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
library(ggpubr)
library(matrixStats)
library(cowplot)
library(this.path)
library(loo)

rm(list = ls())
script_directory <- this.path::this.dir()
setwd(script_directory)


#------- simulated data extraction -------------------------------------------------------------------------------------------------
simulated_data <- readRDS("data/updated_si/simulated_infection.rds")

M_regions <- simulated_data$No_of_regions
sim_week <- simulated_data$week
sim_final_time <- simulated_data$final_time
init_seed <- simulated_data$init_seed
initial_seeding_day <- simulated_data$initial_seeding_day

C <- simulated_data$mobility
sim_Rt <- simulated_data$Rt
sim_inf <- simulated_data$infection
pop <- simulated_data$population

daily_infection_data <- data.frame(region1 = sim_inf[1:sim_final_time],
                                   region2 = sim_inf[(sim_final_time+1):(2*sim_final_time)],
                                   region3 = sim_inf[((2*sim_final_time)+1):(3*sim_final_time)],
                                   time = 1:sim_final_time)

#---- data arrangements for fitting -------------------------------------------------------------------------------------------------------
fitting_final_time <- 350
prediction_horizon <- 0
N <- fitting_final_time + prediction_horizon

#--------required distribution ------------------------------------------------------------------------------
source("data/distributions.R")
dist <- distributions(fitting_final_time)
si <- dist$si
f_death <- dist$f
f_case <- dist$f_case

fitting_start <- initial_seeding_day+1
data_infection <- floor(as.matrix(daily_infection_data %>% filter(time <= N) %>% select(starts_with("Region"))))
len_data <- nrow(data_infection)

stan_data_connected <- list(M_regions = M_regions,
                          C = C,
                          final_time = fitting_final_time,
                          N = N,
                          initial_seeding_day = initial_seeding_day,
                          data_length = len_data,
                          data_inf = data_infection,
                          SI = si,
                          pop = pop,
                          fitting_start = fitting_start)

m <- cmdstan_model("simulated_fitting_region.stan")
start_time <- Sys.time()

fit_connected <-  m$sample(
  data=stan_data_connected,
  iter_sampling = 200,
  iter_warmup = 300,
  parallel_chains = 4,
  chains=4,
  thin=1,
  seed=12345,
  refresh = 40,
  adapt_delta = 0.9,
  max_treedepth = 13,
  init = \() list(init_R = c(1,1,1.5), rw_sd = 0.01))
end_time <- Sys.time()
completion_time <- end_time - start_time
print(completion_time)
# 
summary_fit_connected <- fit_connected$summary()
log_lik_summary <- fit_connected$draws("log_lik")
loo_result <- loo(log_lik_summary)
save(fit_connected,stan_data_connected, file=paste0("results/updated_si/connected_region_fitting.Rdata"))

#----------------------------------------------------------------------------------------------------------------------
#--------------- fitting separately (without mobility) ----------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

C <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow=3, ncol=3)
len_data <- nrow(data_infection)


M_regions <- 3
start_time <- Sys.time()
stan_data_disconnected <- list(M_regions = M_regions,
                               C = C,
                               final_time = fitting_final_time,
                               N = N,
                               initial_seeding_day = initial_seeding_day,
                               data_length = len_data,                        # data_deaths = data_deaths,
                               data_inf = data_infection,                         # data_cases = data_cases,
                               SI = si,
                               pop = pop,
                               fitting_start = fitting_start)

# for (i in 1:M_regions){
  # stan_data_separate <- list(final_time = final_time,
  #                         initial_seeding_day = initial_seeding_day,
  #                         data_length = len_data,
  #                         data_fit = data_infection[,i],
  #                         SI = si,
  #                         f1 = f1, f2 = f2,
  #                         pop = pop[i],
  #                         fitting_start = fitting_start)

  fit_disconnected <-  m$sample(data=stan_data_disconnected,
                                iter_sampling = 200,
                                iter_warmup = 300,
                                parallel_chains = 4,
                                chains=4,
                                thin=1,
                                seed=12345,
                                refresh = 40,
                                adapt_delta = 0.9,
                                max_treedepth = 13,init = \() list(init_R = c(1,1,1.5), rw_sd = 0.01))
  
  end_time <- Sys.time()
  completion_time <- end_time - start_time
  print(completion_time)

  summary_fit_disconnected <- fit_disconnected$summary()
  log_lik_summary <- fit_disconnected$draws("log_lik")
  loo_result <- loo(log_lik_summary)

  save( fit_disconnected,stan_data_disconnected, file=paste0("results/updated_si/disconnected_region_fitting.Rdata"))


#-----------------------------------------------------------------------------------------------------------------------
#------------ national model fitting -----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# aggregated_data <- rowSums(data_infection)
# m <- cmdstan_model("simulated_fitting_national.stan")
# 
# 
# for (i in 1:M_regions){
# stan_data_national <- list(final_time = final_time,
#                             N=N,
#                             initial_seeding_day = initial_seeding_day,
#                             data_length=len_data,
#                             data_fit = data_infection[,i],
#                             SI= si,
#                             pop = pop[i],
#                             fitting_start = fitting_start,
#                             prediction_horizon = prediction_horizon)
# # 
# fit_national <-  m$sample(
#   data=stan_data_national,
#   iter_sampling = 500,
#   iter_warmup = 1200,
#   parallel_chains = 4,
#   chains=4,
#   thin=1,
#   seed=12345,
#   refresh = 40,
#   adapt_delta = 0.8,
#   max_treedepth = 14,init = \() list(init_R = 0,rw_sd = 0.01))
# 
#   summary_fit_national <- fit_national$summary()
#   save(fit_national, stan_data_national, file= paste0("data/fitting_national_forecast",i,".Rdata"))
# 
# }

