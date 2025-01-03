distributions <- function(final_time){

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

convolution <- function(u) (0.0103*f_cached(u))      # ifr is 0.0103.  # this is the ifr for uk 
f[1] = (convolution(1.5) - convolution(0))
for(i in 2:final_time) {
  f[i] = (convolution(i+.5) - convolution(i-.5)) 
}

f_case <- rep(0,final_time)
f_case[1] = integrate(function(x) dgammaAlt(x,mean=5.1, cv=0.86), lower=0, upper=1.5)$value
for (i in 2:final_time){
  f_case[i] <- integrate(function(x) dgammaAlt(x,mean=5.1, cv=0.86), lower=i-0.5, upper=i+0.5)$value
}

delay_distribution <- data.frame(si=si,
                                 f=f,
                                 f_case=f_case)
return(delay_distribution)

}