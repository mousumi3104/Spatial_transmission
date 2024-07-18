library(tidyr)
library(ggplot2)
library(EnvStats)
library(rstan)
library(tidyverse)
library(reshape2)
# library(this.path)
# 
# script_directory <- this.path::this.dir()
# setwd(script_directory)

##################################################################################################################

load("~/OneDrive - National University of Singapore/Singapore/code1/my_model/spatial_transmission/uk_spatial/data/final_pop_2020_ltla.Rdata")
load("~/OneDrive - National University of Singapore/Singapore/code1/my_model/spatial_transmission/uk_spatial/data/uk_ltla_mobility_matrix.Rdata")

sum(colnames(mob_matrix_norm) == pop_2020$area_name)       # to check the ltla's are in order in both the datasets

ltla_names <- colnames(mob_matrix_norm)
pop_2020 <- pop_2020 %>% mutate(area_name = factor(area_name, levels = ltla_names)) %>% arrange(area_name)
pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)), tolower(substring(x,2)))})

england_pop <- sum(pop_2020$population)
########### regions ##############################################################################################

north_east_index <- which(pop_2020$region == "North east")
north_west_index <- which(pop_2020$region == "North west")
yorkshire_index <- which(pop_2020$region == "Yorkshire and the humber")
east_midlands_index <- which(pop_2020$region == "East midlands")
west_midlands_index <- which(pop_2020$region == "West midlands")
east_index <- which(pop_2020$region == "East")
london_index <- which(pop_2020$region == "London")
south_east_index <- which(pop_2020$region == "South east")
south_west_index <- which(pop_2020$region == "South west")

regions <- c("North east","North west","Yorkshire and the humber","East midlands","West midlands","East","London","South east","South west")
pop_percentage <- c(sum(pop_2020$population[pop_2020$region == "North east"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "North west"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "Yorkshire and the humber"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "East midlands"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "West midlands"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "East"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "London"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "South east"])/sum(pop_2020$population),
                sum(pop_2020$population[pop_2020$region == "South west"])/sum(pop_2020$population))
            
pop_region <- data.frame(regions = regions, population = round(pop_percentage*100,1))
# # # # # data arrangements # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

M_ltla <- nrow(pop_2020)            #number of commercial area
pop <-  pop_2020$population

# # # # serial interval and infection to onset distribution # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
final_time <- 600        #total time

#serial interval  (covid-19)

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

# # # # mobility matrix # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Rt_list <- array(1.2, final_time)            #reproduction number 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
seed_time <- 1          # time of days for initial seeding

# # # # # # data simulation # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# for (ind_Rt in seq(length(Rt_list))){
  
  # Rt <- Rt_list[ind_Rt]
  # seed_time <- 1          # time of days for initial seeding
  # init_seed <-  rep(0,M_ltla)
  # init_seed[north_east_index] <- floor(100/length(north_east_index))         #initial seeding

############################
ltla_reaching_time <- data.frame(matrix(0,nrow=M_ltla, ncol=length(regions)))
england_reaching_time <- data.frame(matrix(0,nrow = length(regions), ncol=1))

for (ind_regions in seq(length(regions))){
    Rt <- Rt_list
    init_seed <-  rep(0,M_ltla)

    seed_index <- which(pop_2020$region == regions[ind_regions])   # seeding at the ltlas of that regions
    init_seed[seed_index] <- floor(100/length(seed_index))         # initial seeding
#############################
# ltla_reaching_time <- data.frame(matrix(0,nrow=M_ltla, ncol=length(regions)))
# england_reaching_time <- data.frame(matrix(0,nrow = length(regions), ncol=1))
# 
# for (ind_regions in seq(length(c("Hartlepool")))){   
#   Rt <- Rt_list
#   init_seed <-  rep(0,M_ltla)
#   seed_index <- which(pop_2020$area_name == "Hartlepool")
#   init_seed[seed_index] <- 100
  
  stan_data <- list(M = M_ltla,
                    pop = pop,
                    final_time = final_time,
                    seed_time = seed_time,
                    init_seed = init_seed,
                    Rt_complete = Rt_list,
                    SI = SI,
                    f = f,
                    C = mob_matrix_norm,
                    iar=iar)
  
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  
  m <- rstan::stan_model(file="~/OneDrive - National University of Singapore/Singapore/code1/my_model/spatial_transmission/uk_spatial/time_of_infection/uk_simulation.stan")
  simulated_data = sampling(object=m, data=stan_data,
                            iter=1,
                            chains=1, thin=1, algorithm = "Fixed_param")
  
  y_sim <- rstan::extract(simulated_data)
   
  # # # # calculation of reaching time from the data # # # # # # # # # # # # # # # # # # # # # # ## 
  
  data_cum_inf <- data.frame(y_sim$cumm_sum[1,,])
  colnames(data_cum_inf) <- pop_2020$area_name
  data_cum_inf$time <- 1:final_time
  data_cum_inf$total_inf <- apply(y_sim$cumm_sum[1,,],1,sum)      # total infection for whole England
  
  for (ind in seq_len(M_ltla)) {
    col_name <- paste(pop_2020$area_name[ind])
    individual_inf <- floor(data_cum_inf[[col_name]])

    for (i in 2:final_time) {
      if (individual_inf[i] > pop[ind]*0.01 ) {
        ltla_reaching_time[ind,ind_regions] <- i
        break  # Exit the loop once the first change is found
      }
    }
  }
  for (i in 2:final_time) {
    if (data_cum_inf$total_inf[i] > england_pop*0.01 ) {
      england_reaching_time[ind_regions,] <- i
      break  # Exit the loop once the first change is found
    }
  }
}

colnames(ltla_reaching_time) <- regions        #pop_2020$area_name[100+(1:20)]
colnames(england_reaching_time) <- "time"
ltla_reaching_time$area <- pop_2020$area_name
england_reaching_time$area <- regions

# # # # plot # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

ltla_data_long <- reshape2::melt(ltla_reaching_time, id.vars = "area")
# england_data_long <- reshape2::melt(england_reaching_time, variable.name = "regions",value.name = "value")
min_reaching_time <- min(ltla_data_long$value)
p <- ggplot(ltla_data_long, aes(x =variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  geom_point(data = england_reaching_time, aes(x=regions, y= time ), col="red",cex=5)+
  labs(title = paste0("Arrival time, Rt = ",Rt_list[1]),
       x = NULL,
       y = "Arrival time") +
  theme(plot.title = element_text(size = 20, face = "bold", family = "Arial"),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  ylim(30,NA)

# for(i in 1:nrow(pop_region)) {
#   p <- p + annotate("text", x = pop_region$regions[i], y = min_reaching_time , 
#                     label = pop_region$population[i], size = 5, vjust = 1.5)
# }
p <- p + geom_text(data = pop_region, aes(x = regions, y=35,label = population), size = 5, vjust = 1.5,col="midnightblue")
plot(p)


figure_name <- paste0("reaching_timeRt_",Rt_list[1],".png")
path_to_save <- "~/OneDrive - National University of Singapore/uk_mobility_data/figure/reaching_time"
full_path <- file.path(path_to_save,figure_name)
ggsave(p, filename=full_path, width=8, height=8)

