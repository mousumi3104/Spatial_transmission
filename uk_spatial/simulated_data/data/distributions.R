distributions <- function(final_time){

# si_mean <- 6.5
# si_cv <- 0.62
# 
# x <- rgammaAlt(1e6,si_mean,si_cv)
# 
# si_cum_df <-ecdf(x)
# convolution <- function(u) (si_cum_df(u))
# si[1] = convolution(1.5) - convolution(0)
# for(i in 2:final_time){
#   si[i] = (convolution(i+.5) - convolution(i-.5))
# }

si_data <- read.csv("data/serial_interval.csv")
si <- c(si_data$fit, rep(0,final_time-nrow(si_data)))

#### -------------------------------------------------------------------------------------------
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

###-----------------------------------------------------------------------------------------

f_case <- rep(0,final_time)

f_case_mean <- 5.1
f_case_cv <- 0.86

y <- rgammaAlt(1e6,f_case_mean, f_case_cv)
f_case_cum_df <- ecdf(y) 
convolution <- function(u)(f_case_cum_df(u))

f_case[1] = (convolution(1.5) - convolution(0))
for(i in 2:final_time) {
  f_case[i] = (convolution(i+.5) - convolution(i-.5)) 
}

delay_distribution <- data.frame(si=si,
                                 f=f,           # infection to death
                                 f_case=f_case)
return(delay_distribution)

}