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
#library(rstanarm)
library(this.path)
# 
script_directory <- this.path::this.dir()
setwd(script_directory)


#--------------------------------------------------------------------------------------------------------
M <- 3
week <- 7
final_time <- 70*week
pop <- c(20000000,10000000,15000000)

a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/GBR-estimates.csv")
Rt_1 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_1 <- Rt_1 %>% filter(date >= as.Date("2020-01-01") & date <= as.Date("2021-12-01")) %>% distinct(date, .keep_all = TRUE) # choose only the first option
a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/IND-estimates.csv")
Rt_2 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_2 <- Rt_2 %>% filter(date >= as.Date("2020-01-01") & date <= as.Date("2021-12-01")) %>% distinct(date, .keep_all = TRUE)
a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/ITA-estimates.csv")
Rt_3 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_3 <- Rt_3 %>% filter(date >= as.Date("2020-01-01") & date <= as.Date("2021-12-01")) %>% distinct(date, .keep_all = TRUE)
Rt <- cbind(Rt_1 = 0.15+(Rt_1$Rt),Rt_2 = 0.2+(Rt_2$Rt) ,Rt_3=Rt_3$Rt)
Rt <- data.frame(Rt[1:final_time,])

C <- matrix(c(0.9,0.05,0.05,0.07,0.85,0.08,0.21,0.09,0.7),nrow=3,ncol=3) 

seed_time <- 1
initial_seeding_day <-6 #c(10,5,15)
tau <- rexp(3,0.02)
init_seed <- floor(rexp(3,1/tau))

#--------------------------------------------------------------------------------------------------------
# serial interval : infection to infection (gamma distribution with shape 6.5 and rate 0.62)

inf_to_onset <- si <- rep(0,final_time)
si[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
for (i in 2:final_time){
  si[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
}

# infection to death distribution : infection to onset + onset to death (gamma distribution)
# infection to onset (gamma distribution to mean 5.1 and coefficient of variation 0.86)
#onset to death (gamma distribution to mean 17.8 and coefficient of variation 0.45)

mean1 <- 5.1; cv1 <- 0.86; mean2 <-17.8 ; cv2 <- 0.45;
x1 <- rgammaAlt(1e6,mean1,cv1)
x2 <- rgammaAlt(1e6,mean2,cv2)
f1 <- rep(0,final_time)
f2 <- rep(0,final_time)

f1_cached <- ecdf(x1+x2)  # infection to death    
f2_cached <- ecdf(x1)    # infection to onset
ifr <- 1 #0.0103
iar <- 1 #0.01        #(0.8 high reporting rate, 0.3 low reporting rate)

convolution1 <- function(u) (ifr * f1_cached(u))
f1[1] = (convolution1(1.5) - convolution1(0))
for(i in 2:final_time) {
  f1[i] = (convolution1(i+.5) - convolution1(i-.5)) 
}
convolution2 <- function(u) (iar * f2_cached(u))
f2[1] = (convolution2(1.5) - convolution1(0))
for(i in 2:final_time) {
  f2[i] = (convolution2(i+.5) - convolution2(i-.5)) 
}

day_week_index <- array(0,final_time)
for (t in 1:final_time){ 
  day_week_index[t] = ceiling(t/7) 
}
week <- day_week_index[length(day_week_index)]  

#----------------------------------------------------------------------------------------------------------------
stan_data <- list(M = M,
pop=pop,
final_time=final_time,
initial_seeding_day = initial_seeding_day,
init_seed= c(30,30,30),
SI=si,
f1=f1,f2=f2,
Rt=Rt[,1:M],
C=C,iar=iar,
day_week_index = day_week_index,
week=week)
# # 
stan_data_single_region <- list(pop=pop[1],
                             final_time = final_time,
                             initial_seeding_day = initial_seeding_day,
                             init_seed = 30,
                             SI=si,
                             f1=f1,
                             f2=f2,
                             Rt=Rt$Rt_1,
                             iar=iar)

#----------------------------------------------------------------------------------------------------------------
# m <- cmdstan_model("simulated_data.stan")
m_single_region <- cmdstan_model("simulated_data_single_region.stan")

fit <- m_single_region$sample(
  data =stan_data_single_region,# stan_data,
  iter_sampling = 1,
  chains = 1,
  thin = 1, 
  fixed_param =TRUE)  

# out <- fit$draws(format = "matrix")


daily_infection <- apply(fit$draws("infection",format = "matrix"),2,mean)
plot(daily_infection,type="l")
# points(estimated_inf_mean,col="red")
daily_deaths <- apply(fit$draws("daily_deaths", format = "matrix"),2,mean)
plot(daily_deaths,type="l")
daily_cases <- apply(fit$draws("daily_cases", format = "matrix"),2,mean)
weekly_death <- apply(fit$draws("weekly_deaths",format= "matrix"),2,mean)
# plot(daily_cases)


daily_inf <- data.frame(region1 = daily_infection[1:final_time],
                        region2 = daily_infection[(final_time+1):(2*final_time)],
                        region3 = daily_infection[(2*final_time+1):(3*final_time)],
                        index = 1:final_time)
# 
daily_inf_long <- daily_inf %>% pivot_longer(cols = starts_with("Region"),names_to = "variable", values_to = "value")

ggplot(daily_inf_long, aes(x = index, y = value, color = variable)) +
  geom_line() +
  labs(title = "Simulated daily infection over time",
       x = "Time (day)",
       y = "Daily infection",
       color = "Rt Type") +
  theme_minimal()
# #-----------------------------------------------------------------------------------------------------------
daily_deaths <- data.frame(region1 = daily_deaths[1:final_time],
                           region2 = daily_deaths[(final_time+1):(2*final_time)],
                           region3 = daily_deaths[(2*final_time+1):(3*final_time)],
                           index = 1:final_time)

weekly_deaths <- data.frame(region1 = weekly_death[1:week],
                        region2 = weekly_death[(week+1):(2*week)],
                        region3 = weekly_death[(2*week+1):(3*week)],
                        index = 1:week)
# 
weekly_death_long <- weekly_deaths %>% pivot_longer(cols = starts_with("region"),names_to = "variable", values_to = "value")
daily_death_long <- daily_deaths %>% pivot_longer(cols = starts_with("region"),names_to = "variable", values_to = "value")

ggplot(daily_death_long, aes(x = index, y = value, color = variable)) +
  geom_line() +
  labs(title = "Simulated daily death over time",
       x = "Time (week)",
       y = "daily death",
       color = "death Type") +
  theme_minimal()
# 
daily_cases <- data.frame(region1 = daily_cases[1:final_time],
                        region2 = daily_cases[(final_time+1):(2*final_time)],
                        region3 = daily_cases[(2*final_time+1):(3*final_time)],
                        index = 1:final_time)
# 
# 
save(daily_inf,daily_cases,daily_deaths,stan_data_single_region, file=paste0("data/",'simulated_data_single_region.Rdata'))
