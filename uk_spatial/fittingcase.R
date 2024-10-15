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
load("data/uk_ltla_mobility_matrix.Rdata")

# load("~/mousumi_codes/uk_spatial/data/final_pop_2020_ltla.Rdata")
# load("~/mousumi_codes/uk_spatial/data/england_death_2020.Rdata")
# load("~/mousumi_codes/uk_spatial/data/uk_ltla_mobility_matrix.Rdata")

death_data <- death_data %>% select(all_of(pop_2020$area_name))  
death_data$total_death <- apply(death_data,1, sum)
threshold_week <- which(cumsum(death_data$total_death) > 10 )[1] -1      #the index of that week
inf_start_week <- threshold_week - 4

start_of_year <- as.Date("01-01-2020",format="%d-%m-%y")
start_ind_2020 <- as.integer(strftime(start_of_year,format = "%u"))
first_monday_2020 <- start_of_year - (start_ind_2020 - 1)
inf_start_date <- first_monday_2020 + (inf_start_week*7)

end_of_year <- as.Date("31-12-2020",format = "%d-%m-%y")

final_time <- as.numeric(difftime(end_of_year , inf_start_date, units = "days"))


M_regions <-2#ncol(death_data)-1     # number of region

death_data <- death_data[inf_start_week:(inf_start_week+floor(final_time/7)),1:M_regions]   #(adjust accordingly)
week <- nrow(death_data)

C <- mob_matrix_norm[1:M_regions,1:M_regions]    # mobility matrix 

for (ind in 1:M_regions){
  C[ind,ind] <- 1-(sum(C[,ind])-C[ind,ind])
}

# C <- diag(M_regions)

final_time = 7*week
initial_seeding_day = 6
# initial_seeding = rep(5,M_regions)
pop = pop_2020$population[1:M_regions]

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

convolution <- function(u) (0.0103*f_cached(u))
f[1] = (convolution(1.5) - convolution(0))
for(i in 2:final_time) {
  f[i] = (convolution(i+.5) - convolution(i-.5)) 
}

week_index <- array(0,(final_time-1))
for (t in 1:final_time){
  week_index[t] = floor(t/7)+1      
}
week_index <- c(1,week_index[1:(final_time-1)])

stan_data <- list(M_regions= M_regions,
                  final_time=final_time,
                  W = week,
                  initial_seeding_day=initial_seeding_day,
                  death=death_data,#[seq(1,nrow(death_data),by=2),],
                  SI=si,
                  f=f,
                  pop=pop,C=C,
                  week_index=week_index
                  )     # this is the ifr for uk from the code of swapnil's nature npi

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Example in R using rstan
m <- rstan::stan_model(file="fittingcase_zro_inf.stan")

fit = rstan::sampling(
  object=m,
  data=stan_data,
  iter=1000, 
  warmup=700, 
  chains=1,
  thin=1, 
  seed=12345,
  control = list(adapt_delta = 0.9, max_treedepth = 10))     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
                                                              # opposite for the higher adapt_delta)
                                                              # default adapt_delta=0.80, max.treedepth = 10)
                                                              # max_treedepth is for time efficiency concern (higher -> more time)  
out <- rstan::extract(fit)
parnames <- names(out)
#pairs(fit,pars = "weekly_effect_d")
posterior_fit <- as.array(fit)       # for posterior distribution
summary_fit <- summary(fit)
Rhat <- summary_fit$summary[,"Rhat"]
ess_bulk <- summary_fit$summary[,"n_eff"]

###############################################################################################################

par(mfrow = c(2, 1), mai=c(0.8,1,0.2,0.3))
##.  plotting rt.  ###############################################
rt_samples <- out$Rt
mean_rt <- apply(rt_samples[,,1],2,mean)
ci_rt <- apply(rt_samples,2,quantile,c(0.05,0.95))
plot(seq(length(mean_rt)),mean_rt,type='l',col='red',xlab="",ylab="Rt")
grid(nx = 5, ny =NA ,lty = 2,col = "gray",lwd = 2) 
polygon(c(1:length(mean_rt), rev(1:length(mean_rt))),
        c(ci_rt[1, ], rev(ci_rt[2, ])),
        col = rgb(1, 0, 0, 0.1), border = NA)
abline(a=1,b=0,h=1)


#######.  plotting cases.... the model fitting data ####################################
death_samples <- out$weekly_deaths
mean_deaths1 <- apply(death_samples[,,1],2,mean)
plot(death_data$Hartlepool,xlab="Time(week)",ylab=paste("Weekly deaths","\n", colnames(death_data)[1]))
lines(mean_deaths1,col="red")
ci_deaths <- apply(death_samples,2,quantile,c(0.05,0.95))
grid(nx = 5, ny = NA,lty = 2,col = "gray",lwd = 2) 
polygon(c(1:length(mean_deaths1), rev(1:length(mean_deaths1))),
         c(ci_deaths[1, ], rev(ci_deaths[2, ])),
         col = rgb(1, 0, 0, 0.1), border = NA)


#points(week,incidence_data_week)
par(mfrow = c(1, 1) )
##########################################################
# observed_data_week$time <- week
# ggplot(observed_data_week, aes(x=week,y=Middlesbrough))+
#   geom_point(color="red",shape=19)+
#   geom_point(y=observed_data_week$Middlesbrough, color="blue",shape=19)
#   xlab("time (week)" )+
#   ylab("death_data")


#---------- diagnose ----------------------------------------------------------#
# param <- 10
# no_chain <- 4
# 
# for (ind_par in 1:param){
#   post <- posterior_fit[,,ind_par]
#   
#   for (ind_par in 1:param){ 
#     ggplot(post[])
# }