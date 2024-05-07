library("tidyr")
library("ggplot2")
library("EnvStats")
library("rstan")
library("tidyverse")

# # # #  To get the mobility matrix by simulation # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#source("~/OneDrive - National University of Singapore/Singapore/code1/my_model/agg_mobility_matrix.R", local = TRUE)
#C <- mobility_matrix(14)   # time cutoff 14 = 2pm

# # # # # data arrangements # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

planning_areas <- read.csv("~/OneDrive - National University of Singapore/Singapore/data/datamall/planning_areas/planning_areas.csv")
neg_area <- c("LIM CHU KANG","NORTH-EASTERN ISLANDS","SIMPANG","WESTERN ISLANDS","STRAITS VIEW",
              "CHANGI BAY", "MARINA EAST","MARINA SOUTH")

planning_areas <- subset(planning_areas, !(Name %in% neg_area))
planning_areas$index <- seq(nrow(planning_areas))
poi <- c("BOON LAY","CENTRAL WATER CATCHMENT","CHANGI","PIONEER","SELETAR",
         "SOUTHERN ISLANDS","SUNGEI KADUT","TUAS",
         "WESTERN WATER CATCHMENT","SINGAPORE RIVER","DOWNTOWN CORE",
         "MUSEUM","NEWTON","ORCHARD","TENGAH","MANDAI","PAYA LEBAR")    
poi_index <- which(planning_areas$Name %in% poi)

planning_areas <- rbind(planning_areas[-poi_index,], planning_areas[poi_index,])
planning_areas$index <- seq(nrow(planning_areas))

M_poi <- length(poi)            #number of commercial area
M_res <- nrow(planning_areas) - M_poi       #number of residential area
pop <-  planning_areas$Population[1:M_res]

# # # # serial interval and infection to onset distribution # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
final_time <- 1000
n_simulation <- 20
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

iar = 1

# # # # reproduction number # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Rt_res <- rep(1.5,M_res)
area_index <- planning_areas$index[planning_areas$Name %in% c("JURONG WEST","BEDOK")]
Rt_res[area_index] <- 0.8
Rt_poi <- rep(0.9,M_poi)

#C <- array(0,dim=c(M_res+M_poi,M_res))
#diag(C) <- 1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
seed_time <- 1          # time of days for initial seeding
init_seed <- rep(0,M_res)         #initial seeding
seeding_index <- planning_areas$index[planning_areas$Name == "BEDOK"]
init_seed[seeding_index] <- 1

mobility <- c(1)
plotting_inf_data <- data.frame(matrix(0,nrow=final_time,ncol=length(mobility)))
plotting_case_data <- data.frame(matrix(0,nrow=final_time,ncol=length(mobility)))

# # # # # # data simulation # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
for (ind in seq(length(mobility))){
  
  load("~/OneDrive - National University of Singapore/Singapore/code1/my_model/C_weekdays.RData")  #mobility matrix
  mobility_reduced <- array(mobility[ind],dim=c((M_res+M_poi),M_res))
  
  C <- C_weekdays * mobility_reduced
  diag(C) <- 1-mobility[ind]+diag(C)
  
  stan_data <- list(M_res = M_res,
                    M_poi = M_poi,
                    pop = pop,
                    final_time = final_time,
                    seed_time = seed_time,
                    init_seed = init_seed,
                    Rt_res=Rt_res,
                    Rt_poi=Rt_poi,
                    SI = SI,
                    f = f,
                    C = C,
                    Rt_res_data = Rt_res,
                    Rt_poi_data = Rt_poi,
                    iar=iar,
                    n_simulation=n_simulation)
  
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  
  m <- rstan::stan_model(file="~/OneDrive - National University of Singapore/Singapore/code1/my_model/outbreak.stan")
  
  simulation = sampling(object=m,data=stan_data,
                        iter=1,
                        chains=1, thin=1, algorithm = "Fixed_param")
  
# # # # # extracting data # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  y_sim <- rstan::extract(simulation) 
  
  simulated_Rt <- y_sim$final_Rt[1,,,]
  mean_Rt <- apply(simulated_Rt,c(2,3),mean)[1,]
  
  simulated_data <- y_sim$final_infection[1,,,]
  sim_total_inf <- apply(simulated_data, 1, function(slice) apply(slice,1,sum))
  
  quantile_total <- apply(sim_total_inf, 1, function(row) quantile(row, c(.05, 0.95), na.rm = TRUE))
  lower_total <- quantile_total[1,]
  upper_total <- quantile_total[2,]
  mean_total <- apply(sim_total_inf, 1, mean)

##### plot total_infection and for every individual area  #################################
  
  plot(mean_total, type="l", main="total_infection",ylim = c(min(lower_total), max(upper_total)),
       xlab="Time", ylab="Total infection")
  polygon(c(1:final_time, rev(1:final_time)), 
          c(lower_total, rev(upper_total)), 
          col = rgb(0, 0, 1, alpha = 0.2), border = NA)
  
  # Initialize arrays to store lower and upper bounds of the bands
  lower_band <- array(NA, dim = c(final_time, M_res))
  upper_band <- array(NA, dim = c(final_time, M_res))
  sim_mean <- array(NA, dim= c(final_time, M_res))
  
  for (j in 1:M_res) {  
    for (t in 1:final_time) {
      data_slice <- simulated_data[,t ,j]
      quantiles <- quantile(data_slice, c(0.05, 0.95), na.rm = TRUE)
      lower_band[t, j] <- quantiles[1]
      upper_band[t, j] <- quantiles[2]
      
      sim_mean[t,j] <- mean(simulated_data[,t,j])
    }
  }
  
  for (j in 1:M_res) {  
    plot(sim_mean[,j],type="l", main=planning_areas$Name[j],ylim = c(min(lower_band), max(upper_band)),
         xlab="Time", ylab="Total infection")
    polygon(
      c(1:final_time, rev(1:final_time)), 
      c(lower_band[,j], rev(upper_band[,j])), 
      col = rgb(0, 0, 1, alpha = 0.2), 
      border = NA
    )}
}










# Now lower_band and upper_band contain the lower and upper bounds of the bands for each time point and planning area





# # # # plot # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#data_long <- reshape2::melt(plotting_inf_data, id.vars = "time")
#ggplot(data_long, aes(x = time, y = value, color = variable)) +
#  geom_line(linewidth = 1.5) +  # Adjust line width here
#  scale_color_manual(values = c("deepskyblue2", "coral2", "darkolivegreen4")) +  # Custom colors
#  labs(title = "Effects of mobility ",
#       x = "Time",
#       y = "Total incidence",
#       color = "Mobility Type")+
#  theme(plot.title=element_text(size=24,face="bold", family = "Arial"),
#        axis.text = element_text(size=20,face="bold",vjust=0.5),
#        axis.title = element_text(size = 22,face = "bold"),
#        legend.position = c(0.15, 0.85),
#        legend.text = element_text(size = 16),
#        legend.title= element_text(size=16,face="bold"))




# # # infection data # # # # # # # # # # # # # # # # # # # # # # ## 

#data_inf <- data.frame(matrix(nrow = final_time, ncol = M_res))
#for (ind_res in seq_len(M_res)) {
#  data_inf[, ind_res] <- y_sim$infection[,,ind_res]
#}
#colnames(data_inf) <- paste0(planning_areas$Name[1:M_res])

#data_inf$time <- seq(final_time)
#data_inf <- data_inf %>% select("time",everything())

#data_inf$total_inf <- apply(data_inf[,2:ncol(data_inf)],1,sum)

#col_name <- paste("Mobility","=",mobility[ind])
#colnames(plotting_inf_data)[ind] <- col_name
#plotting_inf_data[,ind] <- data_inf$total_inf
#plotting_inf_data$time <- data_inf$time

# # # case data # # # # # # # # # # # # # # # # # # # # # # ## 
#data_case <- data.frame(matrix(0,nrow = final_time, ncol = M_res))
#for (ind_res in seq_len(M_res)) {
#  data_case[, ind_res] <- y_sim$cases[,,ind_res]
#}
#colnames(data_case) <- paste0(planning_areas$Name[1:M_res])

#data_case$time <- seq(final_time)
#data_case <- data_case %>% select("time",everything())

#data_case$total_case <- apply(data_case[,2:ncol(data_case)],1,sum)

#col_name <- paste("Mobility","=",mobility[ind])
#colnames(plotting_case_data)[ind] <- col_name
#plotting_case_data[,ind] <- data_case$total_case
#plotting_case_data$time <- data_case$time


#}

# # # # plot # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#data_long <- reshape2::melt(plotting_inf_data, id.vars = "time")
#ggplot(data_long, aes(x = time, y = value, color = variable)) +
#  geom_line(linewidth = 1.5) +  # Adjust line width here
#  scale_color_manual(values = c("deepskyblue2", "coral2", "darkolivegreen4")) +  # Custom colors
#  labs(title = "Effects of mobility ",
#       x = "Time",
#       y = "Total incidence",
#       color = "Mobility Type")+
#  theme(plot.title=element_text(size=24,face="bold", family = "Arial"),
#        axis.text = element_text(size=20,face="bold",vjust=0.5),
#        axis.title = element_text(size = 22,face = "bold"),
#        legend.position = c(0.15, 0.85),
#        legend.text = element_text(size = 16),
#        legend.title= element_text(size=16,face="bold"))
