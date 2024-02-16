# this is the basic code for three regions where the mobility flow is represented with the contact matrix "C"
# if C=Identity matrix (as defined below), we get the case where people are not moving to any other population
# population is taken same for all regions (will change accordingly)
# Rt is same for all regions over the time (will change accordingly)
# There is no infection to observation distribution is considered. Only infection is taken into account (will change accordingly)

library(rstan)
library(dplyr)
library(EnvStats)

M <- 3     #number of region
pop <- rep(80000,M)#c(80000,80000,80000)     #population
final_time <- 600     #total time
seed_time <- 5          # time of days for initial seeding
C <- matrix(c(0.91,0.05,0.04,0.05,0.83,0.12,0.10,0.05,0.85),nrow=3,ncol=3)
#C <-  matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3) #if there is no mobile population      

Rt <- matrix(1.5,final_time,M)   #reproduction number  

init_seed<- rep(3,M)      #initial seeding

#serial interval  
SI <- rep(0,final_time)        
SI[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
for (i in 2:final_time){
  SI[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
}

# infection to case distribution
mean1 <- 5.1; cv1 <- 0.86       
x1 <- rgammaAlt(1e6,mean1,cv1)
f <- rep(0,final_time)
f_cached2 <- ecdf(x1)
convolution <- function(u) (f_cached2(u))
f[1] = (convolution(1.5) - convolution(0))
for(i in 2:final_time) {
  f[i] = (convolution(i+.5) - convolution(i-.5)) 
}

iar <- 1    # not required 

#stan data
stan_data <- list(M=M,
                  final_time=final_time,
                  seed_time= seed_time,
                  init_seed=init_seed,
                  SI=SI,
                  f=f,
                  pop=pop,
                  C=C,
                  Rt=Rt,iar=iar)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

# Example in R using rstan
m <- stan_model(file="~/OneDrive - National University of Singapore/Singapore/code1/my_model/spatial_trans.stan")
simulated_data = sampling(object=m,data=stan_data,
              iter=1,
              chains=1, thin=1,algorithm = "Fixed_param")
  
y_sim <- simulated_data %>% 
  as.data.frame %>% 
  select(contains("infection"))

#y_mean <- apply(y_sim,2,mean)
plot(1:final_time,y_sim[1:final_time],col="blue",type="l",lwd=2,xlab="time",ylab="infection",ylim=c(0,max(y_mean)))
lines(1:final_time,y_sim[(final_time+1):(2*final_time)],lwd=2,col="red")
lines(1:final_time,y_sim[(2*final_time+1):(3*final_time)],lwd=2,col="green")
