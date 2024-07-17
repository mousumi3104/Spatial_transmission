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
death_data <- death_data[,1:2]
week <- nrow(death_data)

M_regions <-2#ncol(death_data)     # number of region
C <- mob_matrix_norm[1:2,1:2]    # mobility matrix for two regions
C[1,1] <- 1-C[2,1]
C[2,2] <- 1-C[1,2]

# C <- diag(M_regions)

final_time = 7*week
initial_seeding_day = 1
initial_seeding = rep(5,M_regions)
pop = pop_2020$population[1:2]

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

convolution <- function(u) (f_cached(u))
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
                  initial_seeding=initial_seeding,
                  death=death_data,
                  SI=si,
                  f=f,
                  pop=pop,C=C,
                  week_index=week_index,
                  ifr = 0.010350686).     # this is the ifr for uk from the code of swapnil's nature npi

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Example in R using rstan
m <- rstan::stan_model(file="fittingcase_zro_inf.stan")
fit = rstan::sampling(
  object=m,
  data=stan_data,
  iter=500, 
  warmup=300, 
  chains=2,
  thin=1, 
  seed=12345,
  control = list(adapt_delta = 0.99, max_treedepth = 15))     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
# opposite for the higher adapt_delta)
out <- rstan::extract(fit)
param_names <- names(out)
summary_fit <- summary(fit)
Rhat <- summary_fit$summary[,"Rhat"]
ess_bulk <- summary_fit$summary[,"n_eff"]

###############################################################################################################

par(mfrow = c(2, 1), mai=c(0.4,0.8,0.05,0.3))
##.  plotting rt.  ###############################################
rt_samples <- out[["rt"]]
mean_rt <- apply(rt_samples,2,mean)
ci_rt <- apply(rt_samples,2,quantile,c(0.025,0.975))
plot(seq(length(mean_rt)),mean_rt,type='l',col='red',xaxt="n",ylab="Rt")
grid(nx = 5, ny =NA ,lty = 2,col = "gray",lwd = 2) 
polygon(c(1:length(mean_rt), rev(1:length(mean_rt))),
        c(ci_rt[1, ], rev(ci_rt[2, ])),
        col = rgb(1, 0, 0, 0.1), border = NA)
abline(a=1,b=0,h=1)


#######.  plotting cases.... the model fitting data ####################################
death_samples <- out$weekly_deaths
mean_deaths1 <- apply(death_samples[,,1],2,mean)
mean_deaths2 <- apply(death_samples[,,2],2,mean)
plot(observed_data_week$Hartlepool)
lines(mean_deaths1)


# ci_cases <- apply(cases_samples,2,quantile,c(0.025,0.975))
# plot(seq(length(mean_cases)),mean_cases,type='l',col='red',ylab="Cases")
# grid(nx = 5, ny = NA,lty = 2,col = "gray",lwd = 2) 
# polygon(c(1:length(mean_cases), rev(1:length(mean_cases))),
#         c(ci_cases[1, ], rev(ci_cases[2, ])),
#         col = rgb(1, 0, 0, 0.1), border = NA)
# points(week,incidence_data_week)

##########################################################
# observed_data_week$time <- week
# ggplot(observed_data_week, aes(x=week,y=Hartlepool))+
#   geom_point(color="red",shape=19)+
#   geom_point(y=observed_data_week$Middlesbrough, color="blue",shape=19)
#   xlab("time (week)" )+
#   ylab("death_data")


par(mfrow = c(1, 1) )
