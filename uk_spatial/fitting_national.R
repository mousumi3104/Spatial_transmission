library(rstan)
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

load("data/final_pop_2020_ltla.Rdata")
load("data/england_death_2020.Rdata")       ## weekly data

# load("~/mousumi_codes/uk_spatial/data/final_pop_2020_ltla.Rdata")
# load("~/mousumi_codes/uk_spatial/data/england_death_2020.Rdata")
# load("~/mousumi_codes/uk_spatial/data/uk_ltla_mobility_matrix.Rdata")

pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

death_data <- death_data %>% select(all_of(pop_2020$area_name))  
death_data <- apply(death_data,1,sum)

pop_national <- sum(pop_2020$population)

#--------- start date of the epidemic ----------------------------------------------------------

fitting_start <- which(cumsum(death_data) > 10)[1]     # week when cumulative death exceeds above 10
epidemic_start <- 1#fitting_start - 4    # 30 days before the infection start (4 weeks and 1 week for the delay in data collection)
inf_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", epidemic_start), "-1"))
end_date <- as.Date("31-12-2020",format = "%d-%m-%Y")
week_inf_end <- as.integer(strftime(end_date, format = "%U"))+1

#--------- data arrangement ------------------------------------------------------------------

final_time <- as.numeric(end_date - inf_start_date)+1

death_national <- death_data[epidemic_start : week_inf_end] 
len_data <- length(death_national)

day_week_index <- array(0,final_time)
for (t in 1:final_time){ 
  day_week_index[t] = ceiling(t/7) -(fitting_start-epidemic_start-1)
}
day_week_index[1:((fitting_start-epidemic_start)*7)] = 1
week <- day_week_index[length(day_week_index)]  

initial_seeding_day = 6

pop = pop_national              

#-----------distributions ----------------------------------------------------------------------

si <- rep(0,final_time)
si[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
for (i in 2:final_time){
  si[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
}


mean1 <- 5.1; cv1 <- 0.86; mean2 <-17.8 ; cv2 <- 0.45;
x1 <- rgammaAlt(1e6,mean1,cv1)
x2 <- rgammaAlt(1e6,mean2,cv2)
f <- rep(0,final_time)

f_cached <- ecdf(x1+x2)
ifr <- 0.0103

convolution <- function(u) (ifr*f_cached(u))
f[1] = (convolution(1.5) - convolution(0))
for(i in 2:final_time) {
  f[i] = (convolution(i+.5) - convolution(i-.5)) 
}

stan_data <- list(final_time=final_time,  
                  W = week,       
                  initial_seeding_day=initial_seeding_day,      
                  data_length =len_data,
                  death=death_national,
                  SI=si,   
                  f=f,
                  pop=pop_national,
                  day_week_index = day_week_index,      
                  fitting_start = fitting_start
                  )     

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Example in R using rstan
m <- rstan::stan_model(file="fitting_national.stan")

fit = rstan::sampling(
  object=m,
  data=stan_data,
  iter=1500, 
  warmup=1000, 
  chains=4,
  thin=1, 
  seed=1234,
  control = list(adapt_delta = 0.95, max_treedepth = 12))     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
                                                              # opposite for the higher adapt_delta)
                                                              # default adapt_delta=0.80, max.treedepth = 10)
                                                              # max_treedepth is for time efficiency concern (higher -> more time)  
out <- rstan::extract(fit)

parnames <- names(out)
estimated_infection <-
estimated_weekly_death <- out$weekly_deaths
estimated_rt <- 

#pairs(fit,pars = "weekly_effect_d")
posterior_fit <- as.array(fit)       # for posterior distribution
summary_fit <- summary(fit)
Rhat <- summary_fit$summary[,"Rhat"]
ess_bulk <- summary_fit$summary[,"n_eff"]

save(fit,stan_data, file=paste0("results/",'national_fitting.Rdata'))

#--------------- plot -----------------------------------------------------------------#

complete_dates <- seq.Date(inf_start_date, end_date, by = "day")
mondays <- complete_dates[weekdays(complete_dates) == "Monday"]

estimated_deaths <-  out$weekly_deaths
estimated_deaths_mean <- colMeans(estimated_deaths)
estimated_deaths_li_1 <- colQuantiles(estimated_deaths,prob=0.025)
estimated_deaths_ui_1 <- colQuantiles(estimated_deaths,prob=0.975)
estimated_deaths_li_2 <- colQuantiles(estimated_deaths,prob=0.25)
estimated_deaths_ui_2 <- colQuantiles(estimated_deaths,prob=0.75)

estimated_rt <- out$Rt
estimated_rt_mean <-  colMeans(estimated_rt)
estimated_rt_li_1 <- colQuantiles(estimated_rt,prob=0.025)
estimated_rt_ui_1 <- colQuantiles(estimated_rt,prob=0.975)
estimated_rt_li_2 <- colQuantiles(estimated_rt,prob=0.25)
estimated_rt_ui_2 <- colQuantiles(estimated_rt,prob=0.75)

data_estimated_rt <- data.frame(Rt = estimated_rt_mean,
                        Rt_min_1 = estimated_rt_li_1,
                        Rt_max_1 = estimated_rt_ui_1,
                        Rt_min_2 = estimated_rt_li_2,
                        Rt_max_2 = estimated_rt_ui_2,
                        time = complete_dates)

                        
data_estimated_death <- data.frame(estimated_deaths = estimated_deaths_mean,
                                  death_min_1 = estimated_deaths_li_1,
                                  death_max_1 = estimated_deaths_ui_1,
                                  death_min_2 = estimated_deaths_li_2,
                                  death_max_2 = estimated_deaths_ui_2,
                                  reported_death = death_national,
                                  week = mondays)


data_Rt_95 <- data.frame(time = data_estimated_rt$time, Rt_min = data_estimated_rt$Rt_min_1, 
                            Rt_max = data_estimated_rt$Rt_max_1, key = rep("nintyfive", length(data_estimated_rt$time)))

data_Rt_50 <- data.frame(time = data_estimated_rt$time, Rt_min = data_estimated_rt$Rt_min_2, 
                         Rt_max = data_estimated_rt$Rt_max_2, key = rep("fifty", length(data_estimated_rt$time)))

data_death_95 <- data.frame(time = data_estimated_death$week, death_min = data_estimated_death$death_min_1, 
                            death_max = data_estimated_death$death_max_1, key = rep("nintyfive", length(data_estimated_death$week)))

data_death_50 <- data.frame(time = data_estimated_death$week, death_min = data_estimated_death$death_min_2, 
                            death_max = data_estimated_death$death_max_2, key = rep("fifty", length(data_estimated_death$week)))


data_Rt <- rbind(data_Rt_95, data_Rt_50)
levels(data_Rt$key) <- c("ninetyfive", "fifty")

data_death <- rbind(data_death_95, data_death_50)
levels(data_death$key) <- c("ninetyfive", "fifty")

data_Rt_1 <- data.frame(time = data_estimated_rt$time, Rt = rep(1,length(data_estimated_rt$time)))

p1 <- ggplot(data_estimated_rt)+
  geom_ribbon(data = data_Rt, aes(x=time,ymin = Rt_min, ymax = Rt_max, fill=key))+
  geom_line(data = data_estimated_rt, aes(x=time,y=Rt),color = "black",size=0.8)+
  geom_line(data = data_Rt_1, aes(x=time, y = Rt))+
  xlab("")+
  ylab("Estimated Rt")+
  scale_x_date(date_breaks = "4 week", labels = date_format("%d %b"))+
  scale_fill_manual(name = "", labels = c("50%", "95%"),
                    values = c(alpha("deepskyblue4", 0.35), 
                               alpha("deepskyblue4", 0.25))) + 
  ggtitle("England (National model)")+
  theme_pubr(base_family="sans") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14),
        legend.position = "None")  +
  guides(fill=guide_legend(ncol=1))

plot(p1)

p2 <-ggplot(data_estimated_death)+
  geom_ribbon(data = data_death, aes(x=time, ymin = death_min, ymax = death_max, fill=key))+
  geom_point(data = data_estimated_death, aes(x=week,y=reported_death),color = "coral",size =2)+
  geom_line(data = data_estimated_death, aes(x=week,y=estimated_deaths_mean),color = "black",size=0.8)+
  xlab("")+
  ylab("Estimated Rt")+
  scale_x_date(date_breaks = "4 week", labels = date_format("%d %b"))+
  scale_fill_manual(name = "", labels = c("50%", "95%"),
                    values = c(alpha("deepskyblue4", 0.35), 
                               alpha("deepskyblue4", 0.25))) + 
  scale_shape_manual(values = 16)+
  ggtitle("Weekly reported deaths")+
  theme_pubr(base_family="sans") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14),
        legend.position = "None")  +
  guides(fill=guide_legend(ncol=1))

plot(p2)

p <- plot_grid(p1,p2,ncol=1,nrow=2,rel_heights = c(1,1))
plot(p)

filename = paste0("figures/national_model",".png")
ggsave(filename, plot = p, width = 8, height = 8, units = "in")
