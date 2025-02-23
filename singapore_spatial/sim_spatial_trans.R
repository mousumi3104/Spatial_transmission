# this is the basic code for three regions where the mobility flow is represented with the contact matrix "C"
# if C=Identity matrix (as defined below), we get the case where people are not moving to any other population
# population is taken same for all regions (will change accordingly)
# Rt is same for all regions over the time (will change accordingly)
# There is no infection to observation distribution is considered. Only infection is taken into account (will change accordingly)

library(rstan)
library(dplyr)
library(EnvStats)

M <- 3     #number of region
pop <- rep(80000,M)     #c(80000,60000,30000)    # #population
final_time <- 400     #total time
seed_time <- 5          # time of days for initial seeding
C <- matrix(c(0.81,0.05,0.14,0.05,0.83,0.12,0.60,0.25,0.15),nrow=3,ncol=3)
mobility <- 1      # 0 for no mobility
mobility_reduced <- array(mobility,dim=c(3,3))
C <- C * mobility_reduced
diag(C) <- 1-mobility+diag(C)
#C <-  matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3) #if there is no mobile population      

Rt <- matrix(1,final_time,M)   #reproduction number  

Rt[,1] <- 1.5*Rt[,1]
Rt[,2] <- 1.5*Rt[,2]
Rt[,3] <- 1.5*Rt[,3]

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
m <- stan_model(file="~/OneDrive - National University of Singapore/Singapore/code1/my_model/Spatial_transmission/sim_spatial_trans.stan")
simulated_data = sampling(object=m,data=stan_data,
              iter=1,chains=1, thin=1,algorithm = "Fixed_param")
  
y_sim <- simulated_data %>% 
  as.data.frame %>% 
  select(contains("infection"))

total_infection <- sum(y_sim)
print(total_infection)

print(sprintf("population :%d, total_infection region 1 : %f", pop[1], sum(y_sim[1:400])))
print(sprintf("population :%d, total_infection region 1 : %f", pop[2], sum(y_sim[401:800])))
print(sprintf("population :%d, total_infection region 1 : %f", pop[3], sum(y_sim[801:1200])))

#y_mean <- apply(y_sim,2,mean)
plot(1:final_time,y_sim[1,1:final_time],col="blue",type="l",lwd=2,
     xlab="time",ylab="infection",ylim=c(0,900),main=paste("mobility=", mobility*100,"%"))

lines(1:final_time,y_sim[1,(final_time+1):(2*final_time)],lwd=2,col="red")
lines(1:final_time,y_sim[1,(2*final_time+1):(3*final_time)],lwd=2,col="green")
legend("topright", legend=Rt[1,],col=c("blue","red","green"),lty =1,xpd=TRUE, title = "Rt")

#row <- c(350,420,490) 
#column <- c(490,420,350)

# Add text annotations for matrix values
#for (i in 1:3) {
#  for (j in 1:3) {
#    text(row[j], column[i], labels = C[i, j], cex = 1.5)  # Place matrix values at coordinates (i, j)
#  }
#}
