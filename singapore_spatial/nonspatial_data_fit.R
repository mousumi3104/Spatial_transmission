## to run a code check for the total population and whether there extis a 
##spatial component or not

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
library(rstanarm)

inf <- read.csv(file="~/OneDrive - National University of Singapore/Singapore/code1/my_model/data_sim.csv")

M=1
C <- matrix(c(0.91,0.05,0.04,0.05,0.83,0.12,0.10,0.05,0.85),nrow=3,ncol=3)

data_reg1 <- inf$reg1
data_reg2 <- inf$reg2
data_reg3 <- inf$reg1

##### for non spatial fitting (whole population) ##################
data_inf <- data_reg1 + data_reg2 + data_reg3       #total infection
#data_inf <- total_inf[total_inf!=0]
#data_inf <- data_inf[2:length(data_inf)]
# ######################################################

#data_inf <- total_inf[!apply(total_inf, 1, function(row) any(row == 0)), ]

final_time =length(data_inf)
seed_time <- 10
initial_seeding = data_inf[1:seed_time]
pop = 3*80000      #rep(80000,M)  #3*80000#sum(rep(80000,M))#c(50000,10000,20000)

si <- rep(0,final_time)
si[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
for (i in 2:final_time){
  si[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
}

stan_data <- list(final_time=final_time,
                  initial_seeding_day=seed_time,
                  initial_seeding=initial_seeding,
                  incidence=data_inf,
                  SI=si,
                  pop=pop,M=1)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Example in R using rstan
m <- rstan::stan_model(file="~/OneDrive - National University of Singapore/Singapore/code1/my_model/nonspatial_data_fit.stan")
fit = rstan::sampling(
  object=m,
  data=stan_data,
  iter=5000, 
  warmup=2500, 
  chains=1, 
  thin=1, 
  control = list(adapt_delta = 0.99, max_treedepth = 15))

out <- rstan::extract(fit)
param_names <- names(out)

par(mfrow = c(2, 1), mai=c(0.4,0.8,0.05,0.3))

##plotting rt
rt_samples <- out[["Rt"]]
mean_rt <- apply(rt_samples,2,mean)
ci_rt <- apply(rt_samples,2,quantile,c(0.025,0.975))
plot(seq(length(mean_rt)),mean_rt,type='l',col='black',xaxt="n",ylab="Rt",ylim=c(0.7,2.5))
grid(nx = 5, ny =NA ,lty = 2,col = "gray",lwd = 2) 
polygon(c(1:length(mean_rt), rev(1:length(mean_rt))),
        c(ci_rt[1, ], rev(ci_rt[2, ])),
        col = rgb(0, 0, 0, 0.1), border = NA)
abline(h=2,col="blue")
abline(h=0.8,col="red")
abline(h=1.3,col="green")



##plotting cases.... the model fitting data
cases_samples <- out[["prediction"]]
mean_cases <- apply(cases_samples,2,mean)
ci_cases <- apply(cases_samples,2,quantile,c(0.05,0.95))
plot(seq(length(mean_cases)),mean_cases,type='l',col='red',ylab="total_incidence",ylim=c(0,max(data_inf)))
grid(nx = 5, ny = NA,lty = 2,col = "gray",lwd = 2) 
polygon(c(1:length(mean_cases), rev(1:length(mean_cases))),
        c(ci_cases[1, ], rev(ci_cases[2, ])),
        col = rgb(1, 0, 0, 0.1), border = NA)
points(seq(length(data_inf)),data_inf,type="l")
