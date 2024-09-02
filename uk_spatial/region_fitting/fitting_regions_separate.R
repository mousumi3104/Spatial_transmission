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
load("data/uk_regions_mobility_matrix.Rdata")

# load("~/mousumi_codes/uk_spatial/data/final_pop_2020_ltla.Rdata")
# load("~/mousumi_codes/uk_spatial/data/england_death_2020.Rdata")
# load("~/mousumi_codes/uk_spatial/data/uk_ltla_mobility_matrix.Rdata")

make_plots <- function(data_estimated_rt, data_estimated_death,region){
  
  print(data_estimated_death)
  
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
  
  Rt_threshold <- data.frame(time = data_estimated_rt$time, Rt = rep(1,length(data_estimated_rt$time)))
  
  
  p1 <- ggplot(data_estimated_rt)+
    geom_ribbon(data = data_Rt, aes(x=time,ymin = Rt_min, ymax = Rt_max, fill=key))+
    geom_line(data = data_estimated_rt, aes(x=time,y=Rt),color = "black",linewidth=0.8)+
    geom_line(data = Rt_threshold, aes(x=time, y = Rt))+
    xlab("")+
    ylab("Estimated Rt")+
    scale_x_date(date_breaks = "4 week", labels = date_format("%d %b"))+
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.35), 
                                 alpha("deepskyblue4", 0.25))) + 
    ggtitle(paste(region,"region","(not spatially connected)"))+
    theme_pubr(base_family="sans") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.y = element_text(size = 14),
          legend.position = "None")  +
    guides(fill=guide_legend(ncol=1))
  
  
  p2 <-ggplot(data_estimated_death)+
    geom_ribbon(data = data_death, aes(x=time, ymin = death_min, ymax = death_max, fill=key))+
    geom_point(data = data_estimated_death, aes(x=week,y=reported_death),color = "coral",size =2)+
    geom_line(data = data_estimated_death, aes(x=week,y=estimated_deaths),color = "black",linewidth=0.8)+
    xlab("")+
    ylab("Weekly reported deaths")+
    scale_x_date(date_breaks = "4 week", labels = date_format("%d %b"))+
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.35), 
                                 alpha("deepskyblue4", 0.25))) + 
    scale_shape_manual(values = 16)+
    ggtitle(paste(region,"region"))+
    theme_pubr(base_family="sans") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.y = element_text(size = 14),
          legend.position = "None")  +
    guides(fill=guide_legend(ncol=1))
  
  p <- plot_grid(p1,p2,ncol=1,nrow=2,rel_heights = c(1,1))
  filename = paste0("figures/",region,"_separate.png")
  ggsave(filename, plot = p, width = 8, height = 8, units = "in")
  
}  
#-------------------------------------------------------------------------------------


pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

death_data <- death_data %>% select(all_of(pop_2020$area_name))  

#------- region_pop ----------------------------------------------------------------- 

pop_region <- data.frame(regions = unique(pop_2020$region))

population_region <- rep(0,nrow(pop_region))
for (i in 1:nrow(pop_region)){
  population_region[i] = sum(pop_2020$population[pop_2020$region == pop_region$regions[i]])
}

pop_region$population <- population_region

#------- regions index -----------------------------------------------------------------
region_name <- c("North east","North west","Yorkshire and the humber","East midlands","West midlands","East","London","South east","South west")

north_east_index <- which(pop_2020$region == "North east")
north_west_index <- which(pop_2020$region == "North west")
yorkshire_index <- which(pop_2020$region == "Yorkshire and the humber")
east_midlands_index <- which(pop_2020$region == "East midlands")
west_midlands_index <- which(pop_2020$region == "West midlands")
east_index <- which(pop_2020$region == "East")
london_index <- which(pop_2020$region == "London")
south_east_index <- which(pop_2020$region == "South east")
south_west_index <- which(pop_2020$region == "South west")

#-------- region wise death --------------------------------------------------------------

death_regions <- data.frame(north_east = apply(death_data[,north_east_index],1,sum),
                            north_west = apply(death_data[,north_west_index],1,sum),
                            yorkshire = apply(death_data[,yorkshire_index],1,sum),
                            east_midlands = apply(death_data[,east_midlands_index],1,sum),
                            west_midlands = apply(death_data[,west_midlands_index],1,sum),
                            east = apply(death_data[,east_index],1,sum),
                            london = apply(death_data[,london_index],1,sum),
                            south_east = apply(death_data[,south_east_index],1,sum),
                            south_west = apply(death_data[,south_west_index],1,sum))

#--------- start date of the epidemic ----------------------------------------------------------

death_regions$total_death <- apply(death_data,1, sum)
fitting_start <- which(cumsum(death_regions$total_death) > 10 )[1]       # week when cumulative death exceeds above 10
epidemic_start <- 1
inf_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", epidemic_start), "-1"))
end_date <- as.Date("31-12-2020",format = "%d-%m-%Y")
week_inf_end <- as.integer(strftime(end_date, format = "%U"))+1

#--------- data arrangement ------------------------------------------------------------------

final_time <- as.numeric(end_date - inf_start_date +1)

M_regions <-ncol(death_regions)-1     # number of region

death_regions <- death_regions[,1:M_regions]   #(adjust accordingly)
len_data <- nrow(death_regions)

day_week_index <- array(0,final_time)
for (t in 1:final_time){ 
  day_week_index[t] = ceiling(t/7) -(fitting_start-epidemic_start-1)
}
day_week_index[1:((fitting_start-epidemic_start)*7)] = 1
week <- day_week_index[length(day_week_index)]  

initial_seeding_day = 6
pop = pop_region$population[1:M_regions]          #pop_2020$population[1:M_regions]

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

m <- rstan::stan_model(file="fitting_regions_separate.stan") 

for (i in 1:M_regions){
 stan_data <- list(final_time=final_time,
                  W = week,
                  initial_seeding_day=initial_seeding_day,
                  data_length =len_data,
                  death=death_regions[,i],
                  SI=si,
                  f=f,
                  pop=pop[i],
                  day_week_index = day_week_index,
                  fitting_start = fitting_start)     # this is the ifr for uk from the code of swapnil's nature npi

 options(mc.cores = parallel::detectCores())
 rstan_options(auto_write = TRUE)

 # Example in R using rstan


 fit = rstan::sampling(
   object=m,
   data=stan_data,
   iter=1000, 
   warmup=800, 
   chains=4,
   thin=1, 
   seed=1234,
   control = list(adapt_delta = 0.99, max_treedepth = 15))     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
 # opposite for tigher adapt_delta)
 # default adapt_delta=0.80, max.treedepth = 10                                                              # max_treedepth is for time efficiency concern (higher -> more time)  
 out <- rstan::extract(fit)
   
   complete_dates <- seq.Date(inf_start_date, end_date, by = "day")
   mondays <- complete_dates[weekdays(complete_dates) == "Monday"]
   day <- length(complete_dates)
   week <- length(mondays)
   region <- region_name[i]
   
   estimated_rt <- out$Rt
   data_estimated_rt <- data.frame(Rt = colMeans(estimated_rt),
                                   Rt_min_1 = colQuantiles(estimated_rt,prob=0.025),
                                   Rt_max_1 = colQuantiles(estimated_rt,prob=0.975),
                                   Rt_min_2 = colQuantiles(estimated_rt,prob=0.25),
                                   Rt_max_2 = colQuantiles(estimated_rt,prob=0.75),
                                   time = complete_dates)
   
   estimated_deaths <-  out$weekly_deaths
   data_estimated_death <- data.frame(estimated_deaths = colMeans(estimated_deaths),
                                      death_min_1 = colQuantiles(estimated_deaths,prob=0.025),
                                      death_max_1 = colQuantiles(estimated_deaths,prob=0.975),
                                      death_min_2 = colQuantiles(estimated_deaths,prob=0.25),
                                      death_max_2 = colQuantiles(estimated_deaths,prob=0.75),
                                      reported_death = death_regions[,i],
                                      week = mondays)
   
  results <- make_plots(data_estimated_rt, data_estimated_death,region) 
   
 }

 
