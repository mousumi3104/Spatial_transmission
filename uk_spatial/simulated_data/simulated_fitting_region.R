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

script_directory <- this.path::this.dir()
setwd(script_directory)

#--------------------------------------------------------------------------------------------------------
M_regions <- 3 # nolint
week <- 7
final_time <- 70 * week
pop <- c(20000000, 10000000, 15000000)

a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/GBR-estimates.csv")
Rt_1 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_1 <- Rt_1 %>% arrange(date)
Rt_1 <- Rt_1 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2021-12-31")) %>% distinct(date, .keep_all = TRUE) 
# choose only the first option
Rt_1 <- Rt_1 %>% group_by(date) %>% slice_min(row_number())

a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/IND-estimates.csv")
Rt_2 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_2 <- Rt_2 %>% arrange(date)
Rt_2 <- Rt_2 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2021-12-31")) %>% distinct(date, .keep_all = TRUE)
Rt_2 <- Rt_2 %>% group_by(date) %>% slice_min(row_number())

a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/ITA-estimates.csv")
Rt_3 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_3 <- Rt_3 %>% arrange(date)
Rt_3 <- Rt_3 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2021-12-31")) %>% distinct(date, .keep_all = TRUE)
Rt_3 <- Rt_3 %>% group_by(date) %>% slice_min(row_number())

Rt <- cbind(Rt_1 = 0.1 + Rt_1$Rt,Rt_2 = Rt_2$Rt ,Rt_3=Rt_3$Rt)
Rt <- data.frame(Rt[1:final_time,])


C <- matrix(c(0.7, 0.15, 0.15, 0.07, 0.85, 0.08, 0.21, 0.09, 0.7), nrow=3, ncol=3) 
initial_seeding_day <- 6

si <- rep(0,final_time)
si[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
for (i in 2:final_time){
  si[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
}

mean1 <- 5.1; cv1 <- 0.86; mean2 <- 17.8 ; cv2 <- 0.45;
x1 <- rgammaAlt(1e6, mean1, cv1)
x2 <- rgammaAlt(1e6, mean2, cv2)

f1 <- rep(0, final_time)
f1_cached <- ecdf(x1 + x2)
ifr <- 1
iar <- 1

convolution <- function(u) (ifr * f1_cached(u))
f1[1] <- (convolution(1.5) - convolution(0))
for (i in 2:final_time) {
  f1[i] <- (convolution(i + .5) - convolution(i - .5))
}

f2 <- rep(0, final_time)
f2_cached <- ecdf(x1)
convolution <- function(u) (f2_cached(u))
f2[1] <- (convolution(1.5) - convolution(0))
for(i in 2:final_time) {
  f2[i] <- (convolution(i + .5) - convolution(i - .5))
}

#-------- data generation ------------------------------------------------------------------------------------------------------------------
stan_data_simulation <- list(M = M_regions,
                            pop=pop,
                            final_time=final_time,
                            initial_seeding_day = initial_seeding_day,
                            init_seed= c(30,30,30),
                            SI=si,
                            f1=f1,f2=f2,
                            Rt=Rt[,1:M_regions],
                            C=C,iar=iar,
                            week=week)

m <- cmdstan_model("simulated_region.stan")
simulated_data <- m$sample(data =stan_data_simulation,
                           iter_sampling = 1,
                           chains = 1,
                           thin = 1, 
                           fixed_param = TRUE)  

daily_infection <- apply(simulated_data$draws("infection",format = "matrix"),2,mean)
daily_deaths <- apply(simulated_data$draws("daily_deaths", format = "matrix"),2,mean)
daily_cases <- apply(simulated_data$draws("daily_cases", format = "matrix"),2,mean)
weekly_death <- apply(simulated_data$draws("weekly_deaths",format= "matrix"),2,mean)

daily_deaths_data <- data.frame(region1 = daily_deaths[1:final_time],
                                region2 = daily_deaths[(final_time+1):(2*final_time)],
                                region3 = daily_deaths[((2*final_time)+1):(3*final_time)],
                                time = 1:final_time)

daily_infection_data <- data.frame(region1 = daily_infection[1:final_time],
                                   region2 = daily_infection[(final_time+1):(2*final_time)],
                                   region3 = daily_infection[((2*final_time)+1):(3*final_time)],
                                   time = 1:final_time)

daily_cases_data <- data.frame(region1 = daily_cases[1:final_time],
                               region2 = daily_cases[(final_time+1):(2*final_time)],
                               region3 = daily_cases[((2*final_time)+1):(3*final_time)],
                               time = 1:final_time)

weekly_death_data <- data.frame(region1 = weekly_death[1:week],
                                region2 = weekly_death[(week+1):(2*week)],
                                region3 = weekly_death[((2*week)+1):(3*week)])

# ---------------- plot true Rt / true daily infection------------------------------------------------------

Rt$index <- 1:nrow(Rt)
Rt_long <- Rt %>%
  pivot_longer(
    cols = starts_with("Rt"),  # Select columns that start with 'Rt'
    names_to = "variable",      # Name for the new 'variable' column
    values_to = "value"         # Name for the new 'value' column
  )

ggplot(Rt_long, aes(x = index, y = value, color = variable)) +
geom_line() +
geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Add horizontal line at y=1
labs(title = "Rt Values Over Time",
      x = "Time (day)",
      y = "Rt Value",
      color = "Rt Type") +
theme_minimal()

inf_data <- daily_infection_data
inf_data$index <- 1:nrow(inf_data)
inf_datalong <- inf_data %>% pivot_longer(cols = starts_with("Region"), names_to = "variable", values_to = "value")

ggplot(inf_datalong, aes(x = index, y = value, color = variable)) +
geom_point() +
labs(title = "Daily_infection over time",
      x = "Time (day)",
      y = "Daily infection",
      color = "Region Type") +
theme_minimal()

save(Rt,daily_infection_data,file = "data/simulated_data.Rdata")

#----------------------------------------------------------------------------------------------------------
#------------ fitting with connected model ----------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------

fitting_start <- 10
data_infection <- floor(as.matrix(daily_infection_data[, 1:M_regions]))
len_data <- nrow(daily_infection_data)

# day_week_index <- c(rep(1,(fitting_case_start -1)*7), rep(2:((final_time/7)),each=7))
# day_week_index <- day_week_index[1:final_time]
# week <- day_week_index[length(day_week_index)]

stan_data_connected <- list(M_regions = M_regions,
                          C = C,
                          final_time = final_time,
                          initial_seeding_day = initial_seeding_day,
                          data_length = len_data,
#                           # data_deaths = data_deaths,
                          data_inf = data_infection,
#                           # data_cases = data_cases,
                          SI = si,
                          f1 = f1, f2 = f2,
                          pop = pop,
                          fitting_start = fitting_start)

m <- cmdstan_model("simulated_fitting_region.stan")

fit_connected <-  m$sample(
  data=stan_data_connected,
  iter_sampling = 500,
  iter_warmup = 1500,
  parallel_chains = 4,
  chains=4,
  thin=1,
  seed=12345,
  refresh = 40,
  adapt_delta = 0.95,
  max_treedepth = 10)

summary_fit_connected <- fit_connected$summary()

# save(fit_connected,stan_data_connected, file=paste0("data/",'connected_region_fitting.Rdata'))

#----------------------------------------------------------------------------------------------------------------------
#--------------- fitting separately (without mobility) ----------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

C <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow=3, ncol=3) 
len_data <- nrow(data_infection)


M_regions <- 3
stan_data_connected <- list(M_regions = M_regions,
                            C = C,
                            final_time = final_time,
                            initial_seeding_day = initial_seeding_day,
                            data_length = len_data,
                            #                           # data_deaths = data_deaths,
                            data_inf = data_infection,
                            #                           # data_cases = data_cases,
                            SI = si,
                            f1 = f1, f2 = f2,
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

  fit_disconnected <-  m$sample(data=stan_data_connected,
                                iter_sampling = 500, 
                                iter_warmup = 1500, 
                                parallel_chains = 4,
                                chains=4,
                                thin=1,
                                seed=12345,
                                refresh = 40,
                                adapt_delta = 0.95,
                                max_treedepth = 10)
  
  est_Rt <- fit_disconnected$draws("Rt",format = "matrix")
  est_inf <- fit_disconnected$draws("infection",format = "matrix")

  # assign(paste0("stan_data_separate_",i), stan_data_separate)
  # assign(paste0("fit_disconnected_",i), fit_disconnected)
  # save( fit_disconnected, file=paste0("data/disconnected_region_fitting.Rdata"))


#-----------------------------------------------------------------------------------------------------------------------
#------------ national model fitting -----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

aggregated_data <- rowSums(data_infection)
m <- cmdstan_model("simulated_fitting_national.stan")

stan_data_national <- list(final_time = final_time,
                            initial_seeding_day = initial_seeding_day,
                            data_length=len_data,
                            data_fit = aggregated_data,
                            SI= si,
                            f1=f1,
                            f2=f2,
                            pop = sum(pop),
                            fitting_start = fitting_start)
# 
fit_national <-  m$sample(
  data=stan_data_national,
  iter_sampling = 500,
  iter_warmup = 1500,
  parallel_chains = 4,
  chains=4,
  thin=1,
  seed=12345,
  refresh = 40,
  adapt_delta = 0.95,
  max_treedepth = 10)

  summary_fit_national <- fit_national$summary()
  # save(fit_national, stan_data_national, file= paste0("data/fitting_national.Rdata"))


