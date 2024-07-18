library("tidyr")
library("ggplot2")
library("EnvStats")
library("rstan")
library("tidyverse")

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

singapore_pop <- sum(planning_areas_res$Population)

regions <- unique(planning_areas$Region[1:M_res])
pop_percentage <- c(sum(planning_areas_res$Population[planning_areas_res$Region == "EAST REGION"])/singapore_pop,
                    sum(planning_areas_res$Population[planning_areas_res$Region == "WEST REGION"])/singapore_pop,
                    sum(planning_areas_res$Population[planning_areas_res$Region == "CENTRAL REGION"])/singapore_pop,
                    sum(planning_areas_res$Population[planning_areas_res$Region == "NORTH-EAST REGION"])/singapore_pop,
                    sum(planning_areas_res$Population[planning_areas_res$Region == "NORTH REGION"])/singapore_pop)

pop_region <- data.frame(regions = regions, population = round(pop_percentage*100),1)
M_poi <- length(poi)            #number of commercial area
M_res <- nrow(planning_areas) - M_poi       #number of residential area
pop <-  planning_areas$Population[1:M_res]

planning_areas_res <- planning_areas[1:M_res,]


# # # # serial interval and infection to onset distribution # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
final_time <- 600        #total time

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

seed_time <- 1          # time of days for initial seeding

load("~/OneDrive - National University of Singapore/Singapore/code1/my_model/mobility_matrix.RData")  #mobility matrix

Rt_res <- 1.5
Rt_com <- 1.5
Rt_res_standata <- Rt_res*matrix(Rt_res,final_time,M_res)
Rt_com_standata <- Rt_com*matrix(Rt_com,final_time,M_poi)

planning_area_reach_time <- data.frame(matrix(0,nrow=M_res,ncol=length(regions)))
singapore_reach_time <- data.frame(matrix(0,nrow=length(regions),ncol=1))

#------------------ data simulation ---------------------------------------------------------

m <- rstan::stan_model(file="~/OneDrive - National University of Singapore/Singapore/code1/my_model/sing_simulation.stan")

for (ind_region in seq(length(regions))){
  
  init_seed <- rep(0,M_res)         #initial seeding
  
  seeding_index <- planning_areas_res$index[planning_areas_res$Region == regions[ind_region]]
  init_seed[seeding_index] <- floor(50/length(seeding_index))
  
  stan_data <- list(M_res = M_res,
                    M_poi = M_poi,
                    pop = pop,
                    final_time = final_time,
                    seed_time = seed_time,
                    init_seed = init_seed,
                    Rt_res = Rt_res_standata,
                    Rt_poi = Rt_com_standata,
                    SI = SI,
                    f = f,
                    C_day = C_weekdays,
                    C_end = C_weekends,
                    Rt = Rt,
                    iar=iar,
                    n_simulation=1)
  
  
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  
  simulated_data = sampling(object=m,data=stan_data,
                            iter=1,
                            chains=1, thin=1,algorithm = "Fixed_param")
  
  
  y_sim <- rstan::extract(simulated_data) 
  
  #--------------- infection data -------------------------------------------------
  
  data_cum_inf <- data.frame(y_sim$cumm_sum[1,,])
  colnames(data_cum_inf) <- planning_areas_res$Name
  data_cum_inf$time <- seq(final_time)
  data_cum_inf$singapore <- apply(y_sim$cumm_sum[1,,],1,sum)
  
  for (ind in seq_len(M_res)) {
    res_name <- planning_areas_res$Name[ind]
    individual_inf <- floor(data_cum_inf[[res_name]])
    for (t in 2:final_time){
      if (individual_inf[t] > planning_areas_res$Population[ind]*0.01){
        planning_area_reach_time[ind,ind_region] <- t
        break
      }
    }
  }
  for (t in 2: final_time){
    if (data_cum_inf$singapore[t] > singapore_pop * 0.01){
      singapore_reach_time[ind_region] <- t
      break
    }
  }
}

colnames(planning_area_reach_time) <- regions  
colnames(singapore_reach_time) <- "time"
planning_area_reach_time$area <- planning_areas_res$Name
singapore_reach_time$area <- regions

#------------ plot -------------------------------------------------------------------------- # 

planarea_data_long <- reshape2::melt(planning_area_reach_time, id.vars = "area")
# england_data_long <- reshape2::melt(england_reaching_time, variable.name = "regions",value.name = "value")
min_reaching_time <- min(planarea_data_long$value)
p <- ggplot(planarea_data_long, aes(x =variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  geom_point(data = singapore_reach_time, aes(x=regions, y= time ), col="red",cex=5)+
  labs(title = paste0("Arrival time, Rt = ",Rt_res),
       x = NULL,
       y = "Arrival time") +
  theme(plot.title = element_text(size = 20, face = "bold", family = "Arial"),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  ylim(30,NA)

p <- p + geom_text(data = pop_region, aes(x = regions, y=31,label = population), size = 5, vjust = 1.5,col="midnightblue")
plot(p)


figure_name <- paste0("reaching_timeRt_",Rt_res,".png")
path_to_save <- "~/OneDrive - National University of Singapore/Singapore/code1/my_model"
full_path <- file.path(path_to_save,figure_name)
ggsave(p, filename=full_path, width=8, height=8)


  

  
