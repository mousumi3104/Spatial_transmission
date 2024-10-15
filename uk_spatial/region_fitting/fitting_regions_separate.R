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
#library(rstanarm)
library(this.path)
#
script_directory <- this.path::this.dir()
setwd(script_directory)

# load("data/final_pop_2020_ltla.Rdata")
# load("data/england_death_2020.Rdata")       ## weekly data
# load("data/uk_regions_mobility_matrix.Rdata")
# mobility_change <- readRDS("data/mobility_change.rds")
# 
# # --------------------------------------------------------------------------------
# pop_2020$region <- sapply(pop_2020$region, function(x){
#   paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
# 
# death_data <- death_data %>% select(all_of(pop_2020$area_name))
# 
# #------- region_pop -----------------------------------------------------------------
# 
# pop_region <- data.frame(regions = unique(pop_2020$region))
# 
# population_region <- rep(0,nrow(pop_region))
# for (i in 1:nrow(pop_region)){
#   population_region[i] = sum(pop_2020$population[pop_2020$region == pop_region$regions[i]])
# }
# 
# pop_region$population <- population_region
# 
# #------- regions index -----------------------------------------------------------------
# region_name <- c("North east","North west","Yorkshire and the humber","East midlands","West midlands","East","London","South east","South west")
# 
# north_east_index <- which(pop_2020$region == "North east")
# north_west_index <- which(pop_2020$region == "North west")
# yorkshire_index <- which(pop_2020$region == "Yorkshire and the humber")
# east_midlands_index <- which(pop_2020$region == "East midlands")
# west_midlands_index <- which(pop_2020$region == "West midlands")
# east_index <- which(pop_2020$region == "East")
# london_index <- which(pop_2020$region == "London")
# south_east_index <- which(pop_2020$region == "South east")
# south_west_index <- which(pop_2020$region == "South west")
# 
# #-------- region wise death --------------------------------------------------------------
# 
# death_regions <- data.frame(north_east = apply(death_data[,north_east_index],1,sum),
#                             north_west = apply(death_data[,north_west_index],1,sum),
#                             yorkshire = apply(death_data[,yorkshire_index],1,sum),
#                             east_midlands = apply(death_data[,east_midlands_index],1,sum),
#                             west_midlands = apply(death_data[,west_midlands_index],1,sum),
#                             east = apply(death_data[,east_index],1,sum),
#                             london = apply(death_data[,london_index],1,sum),
#                             south_east = apply(death_data[,south_east_index],1,sum),
#                             south_west = apply(death_data[,south_west_index],1,sum))
# 
# #--------- start date of the epidemic ----------------------------------------------------------
# 
# death_regions$total_death <- apply(death_data,1, sum)
# fitting_start <- which(cumsum(death_regions$total_death) > 10 )[1]       # week when cumulative death exceeds above 10
# infection_gen_time <- 6
# epidemic_start <- fitting_death_start - infection_gen_time
# inf_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", epidemic_start), "-1"))
# end_date <- as.Date("31-12-2020",format = "%d-%m-%Y")
# week_inf_end <- as.integer(strftime(end_date, format = "%U"))+1
# 
# #--------- data arrangement ------------------------------------------------------------------
# 
# final_time <- as.numeric(end_date - inf_start_date +1)
# 
# M_regions <-ncol(death_regions)-1     # number of region
# 
# death_regions <- death_regions[epidemic_start:week_inf_end,1:M_regions]   #(adjust accordingly)
# len_data <- nrow(death_regions)
# 
# day_week_index <- array(0,final_time)
# for (t in 1:final_time){
#   day_week_index[t] = ceiling(t/7) -(fitting_start-epidemic_start-1)
# }
# day_week_index[1:((fitting_start-epidemic_start)*7)] = 1
# week <- day_week_index[length(day_week_index)]
# 
# initial_seeding_day = 6
# pop = pop_region$population[1:M_regions]          #pop_2020$population[1:M_regions]
# 
# si <- rep(0,final_time)
# si[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
# for (i in 2:final_time){
#   si[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
# }
# 
# 
# mean1 <- 5.1; cv1 <- 0.86; mean2 <-17.8 ; cv2 <- 0.45;
# x1 <- rgammaAlt(1e6,mean1,cv1)
# x2 <- rgammaAlt(1e6,mean2,cv2)
# f <- rep(0,final_time)
# 
# f_cached <- ecdf(x1+x2)
# 
# convolution <- function(u) (0.0103*f_cached(u))
# f[1] = (convolution(1.5) - convolution(0))
# for(i in 2:final_time) {
#   f[i] = (convolution(i+.5) - convolution(i-.5))
# }
# 
# #-------- arranging google mobility ---------------------------------------------------------------------------
# 
# mobility_change <- mobility_change %>%
#   group_by(area_name, date, region) %>%
#   summarize(avg_retail_and_recreation_percent_change_from_baseline = mean(retail_and_recreation_percent_change_from_baseline, na.rm = TRUE),
#             ave_grocery_and_pharmacy_percent_change_from_baseline = mean(grocery_and_pharmacy_percent_change_from_baseline, na.rm = TRUE),
#             ave_parks_percent_change_from_baseline = mean(parks_percent_change_from_baseline, na.rm = TRUE),
#             ave_transit_stations_percent_change_from_baseline = mean(transit_stations_percent_change_from_baseline, na.rm = TRUE),
#             ave_workplaces_percent_change_from_baseline = mean(workplaces_percent_change_from_baseline, na.rm = TRUE))
# 
# mobility_change_region <- mobility_change %>%
#   group_by(region, date) %>%
#   summarize(avg_retail_and_recreation_percent_change_from_baseline = mean(avg_retail_and_recreation_percent_change_from_baseline, na.rm = TRUE),
#             ave_grocery_and_pharmacy_percent_change_from_baseline = mean(ave_grocery_and_pharmacy_percent_change_from_baseline, na.rm = TRUE),
#             ave_parks_percent_change_from_baseline = mean(ave_parks_percent_change_from_baseline, na.rm = TRUE),
#             ave_transit_stations_percent_change_from_baseline = mean(ave_transit_stations_percent_change_from_baseline, na.rm = TRUE),
#             ave_workplaces_percent_change_from_baseline = mean(ave_workplaces_percent_change_from_baseline, na.rm = TRUE),.groups = "drop")
# 
# mobility_change_region$date <- as.Date(mobility_change_region$date, format = "%Y-%m-%d")
# 
# if (min(mobility_change_region$date) > inf_start_date){
#   new_dates <- seq.Date(from = inf_start_date, to = min(mobility_change_region$date), by = "day")
#   new_rows <- expand.grid(region = unique(mobility_change_region$region), date = new_dates)
#   new_rows <- new_rows %>% mutate(avg_retail_and_recreation_percent_change_from_baseline = 0,
#                                   ave_grocery_and_pharmacy_percent_change_from_baseline = 0,
#                                   ave_parks_percent_change_from_baseline = 0,
#                                   ave_transit_stations_percent_change_from_baseline = 0,
#                                   ave_workplaces_percent_change_from_baseline = 0)
# }
# mobility_change_region <- rbind(new_rows, mobility_change_region)
# desired_order_region <- pop_region$regions
# mobility_change_region$region <- factor(mobility_change_region$region, levels = desired_order_region)
# mobility_change_region <- mobility_change_region %>% arrange(region)
# 
# gmobility <- array(NA, dim = c(final_time,5,M_regions))
# for (i in 1:M_regions){
#   gmobility[,,i] <- as.matrix(mobility_change_region[(((i-1)*(final_time))+1):(i*(final_time)),3:7])
# }
# 
# for (m in 1:M_regions) {
#   for (i in 1:4){
#     for (t in 1:(final_time)){
#       if (is.nan(gmobility[t,i,m])){
# 
#         prev_val <- if (t > 1) gmobility[max(which(!is.nan(gmobility[1:(t-1),i,m] ))),i,m]
#         next_val <- if (t < final_time) gmobility[(min(which(!is.nan(gmobility[(t+1):final_time,i,m]))) + t),i,m]
#         gmobility[t,i,m] <- (prev_val + next_val)/2
#       }
#     }
#   }
# }
# 
# # ---------lockdown effect ----------------------------------------
# lockdown1_started <- as.Date("2020-03-23", format = "%Y-%m-%d")
# lockdown1_lifted <- as.Date("2020-05-10", format = "%Y-%m-%d")
# 
# lockdown2_started <- as.Date("2020-11-05", format = "%Y-%m-%d")
# lockdown2_lifted <- as.Date("2020-12-02", format = "%Y-%m-%d")
# 
# lockdown_index <-  data.frame(date = seq.Date(from = inf_start_date,as.Date("2020-12-31",format = "%Y-%m-%d"),by="day"))
# lockdown_index$r_mobility_index <- rep(0,nrow(lockdown_index))
# lockdown_index$r_mobility_index[lockdown_index$date >= lockdown1_started & lockdown_index$date <= lockdown1_lifted] <- 1
# lockdown_index$r_mobility_index[lockdown_index$date >= lockdown2_started & lockdown_index$date <= lockdown2_lifted] <- 1
# 
# lockdown_index$g_mobility_index <-  rep(3,nrow(lockdown_index))
# lockdown_index$g_mobility_index[lockdown_index$date <= lockdown1_lifted] <- 1
# lockdown_index$g_mobility_index[lockdown_index$date >= lockdown1_lifted & lockdown_index$date <= lockdown2_lifted] <- 2
# 
# positions <- which(lockdown_index$r_mobility_index == 1)
# first_in_groups <- positions[c(TRUE, diff(positions) != 1)]
# first_lockdown_start <- first_in_groups[1]
# second_lockdown_start <- first_in_groups[2]
# first_lockdown_end <- which(lockdown_index$g_mobility_index == 2)[1]
# second_lockdown_end <- which(lockdown_index$g_mobility_index == 3)[1]
# 
# lockdown_index$L1 <- ifelse(lockdown_index$g_mobility_index == 1,1,0)
# lockdown_index$L2 <- ifelse(lockdown_index$g_mobility_index == 2,1,0)
# lockdown_index$L3 <- ifelse(lockdown_index$g_mobility_index == 3,1,0)
# 
# lockdown_index <- lockdown_index %>% select(- g_mobility_index, -date)
# lockdown_index <- lockdown_index[1:final_time,]

#------------------------------------------------

source("data/stan_data_arrangements.R")
stan_data_connected <- stan_data_arrangements(death_threshold = 10, script_directory)
M_regions <- stan_data_connected$M_regions

m <- cmdstan_model("fitting_regions_separate.stan")


for (i in 1:M_regions){
 stan_data <- list(final_time= stan_data_connected$final_time,
                   W = stan_data_connected$W,
                   initial_seeding_day=stan_data_connected$initial_seeding_day,
                   data_length = stan_data_connected$death_data_length,
                   death= stan_data_connected$death[,i],
                   SI=stan_data_connected$SI,
                   f=stan_data_connected$f,
                   pop=stan_data_connected$pop[i],
                   day_week_index = stan_data_connected$day_week_index,
                   fitting_start = stan_data_connected$fitting_death_start,
                   gmobility = stan_data_connected$gmobility[,i],
                   II = as.matrix(stan_data_connected$I) )

 fit = m$sample(
   data=stan_data,
   iter_sampling = 500,
   iter_warmup =1200, 
   parallel_chains = 4,
   # threads_per_chain = 2,
   chains=4,
   thin=1, 
   seed=1234,
   refresh = 40,
   adapt_delta = 0.99, 
   max_treedepth = 13,init = \() list(mu = 3.28, 
                                      initial_seeding = 5,
                                      tau = 0.01,
                                      phi = 20,
                                      gamma = 0.5))       # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution

   Rt <- fit$draws("Rt",format="matrix")
   inf <- fit$draws("infection",format="matrix")
   death <- fit$draws("weekly_deaths",format="matrix")
   summary_fit <- fit$summary()
 
 save(fit,Rt,inf,death,summary_fit,file = paste0("results/copy/region_disconnected_xyz",i,".rds"))
   
 }

final_time <- stan_data$final_time
load("results/region_connected_rt_xyz.Rdata")
final_time <- stan_data_connected$final_time
Rt <- fit_connected$draws("Rt",format = "matrix")
mean_Rt <- colMeans(Rt)
plot(mean_Rt[(((i-1)*final_time)+1):(i*final_time)],col="red")
# # Rt_min_1 = colQuantiles(Rt,prob=0.025)
# # Rt_max_1 = colQuantiles(Rt,prob=0.975)
# # lines(Rt_min_1,col="red")
# lines(Rt_max_1,col="red")
# # 
# # 
# # 
load(paste0("results/copy/region_disconnected_xyz",i,".rds"))
Rt <- fit$draws("Rt",format = "matrix")
mean_Rt <- colMeans(Rt)
points(mean_Rt,col="black")
# # Rt_min_1 = colQuantiles(Rt,prob=0.025)
# # Rt_max_1 = colQuantiles(Rt,prob=0.975)
# # lines(Rt_min_1)
# # lines(Rt_max_1)
# # 
# # abline(v=c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end))
# # abline(h=1)
# 
# z <- fit$draws("z")
# hist(z)
# plot(x[,1,1])
# points(x[,2,1],col="red")
# points(x[,3,1],col="green")
# points(x[,4,1],col="blue")
# 
# w_e <- fit$draws("weekly_effect")
# hist(w_e[,,1:36],breaks = 20)








