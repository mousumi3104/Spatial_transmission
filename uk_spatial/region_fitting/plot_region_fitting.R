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
library(cowplot)
library(ISOweek)
library(matrixStats)
library(ggpubr)
library(this.path)

script_directory <- this.path::this.dir()
setwd(script_directory)

#----------------------------------- data ------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
load("results/region_fitting.Rdata")
out <- rstan::extract(fit)
death_regions <- stan_data$death

#----------------------------------- plot function -----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
make_plots <- function(i,data_estimated_rt, data_estimated_death,region){
  
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
    ggtitle(paste(region,"region"))+
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

  return(list(p1,p2))
}  

#---------------------------- data and calling function -----------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

epidemic_start <- 1
inf_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", epidemic_start), "-1"))
end_date <- as.Date("31-12-2020",format = "%d-%m-%Y")
complete_dates <- seq.Date(inf_start_date, end_date, by = "day")
mondays <- complete_dates[weekdays(complete_dates) == "Monday"]

region_name <- c("North east","North west","Yorkshire and the humber","East midlands","West midlands","East","London","South east","South west")

data_regions <- data.frame()

plots <- list()

for (i in 1:5){
  
  day <- length(complete_dates)
  week <- length(mondays)
  region <- region_name[i]
  
  estimated_rt <- out$Rt[,,i]
  data_estimated_rt <- data.frame(Rt = colMeans(estimated_rt),
                                  Rt_min_1 = colQuantiles(estimated_rt,prob=0.025),
                                  Rt_max_1 = colQuantiles(estimated_rt,prob=0.975),
                                  Rt_min_2 = colQuantiles(estimated_rt,prob=0.25),
                                  Rt_max_2 = colQuantiles(estimated_rt,prob=0.75),
                                  time = complete_dates)
  
  estimated_deaths <-  out$weekly_deaths[,,i]
  data_estimated_death <- data.frame(estimated_deaths = colMeans(estimated_deaths),
                                     death_min_1 = colQuantiles(estimated_deaths,prob=0.025),
                                     death_max_1 = colQuantiles(estimated_deaths,prob=0.975),
                                     death_min_2 = colQuantiles(estimated_deaths,prob=0.25),
                                     death_max_2 = colQuantiles(estimated_deaths,prob=0.75),
                                     reported_death = death_regions[,i],
                                     week = mondays)
  
  result <- make_plots(i,data_estimated_rt, data_estimated_death,region) 
  plots <- c(plots,result)
}

combined_plot <- plot_grid(plotlist = plots,ncol=2,nrow=5)
plot(combined_plot)

# ggsave(filename = paste0("figures/region_fitting1_5.png"), plot = combined_plot, width=16, height=16, units="in")


#----------------------- region and national model comparison (Rt) ----------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

load("results/region_fitting.Rdata")    # region fitting data

epidemic_start <- 1
inf_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", epidemic_start), "-1"))
end_date <- as.Date("31-12-2020",format = "%d-%m-%Y")
complete_dates <- seq.Date(inf_start_date, end_date, by = "day")


region_name <- c("North east","North west","Yorkshire and the humber","East midlands","West midlands","East","London","South east","South west")
estimated_Rt_region <- data.frame(matrix(nrow = dim(out$Rt)[2] , ncol=length(region_name)))
colnames(estimated_Rt_region) <- region_name

region_colors <- c("North east" = "#e41a1c", 
                   "North west" = "#377eb8", 
                   "Yorkshire and the humber" = "#4daf4a", 
                   "East midlands" = "#984ea3", 
                   "West midlands" = "#ff7f00", 
                   "East" = "#ffff33", 
                   "London" = "#a65628", 
                   "South east" = "#f781bf", 
                   "South west" = "#999999")

for (i in 1:length(region_name)){
  estimated_Rt_region[[i]] <- colMeans(out$Rt[,,i])
}
estimated_Rt_region$time <- complete_dates

melted_Rt_region <- reshape2::melt(estimated_Rt_region, id.vars = "time", variable.name = "region",value.name = "Rt")   # id.vars = x axis, variable, legend, value = yaxis

load("results/national_fitting.Rdata")      # national fitting data
out_national <- rstan::extract(fit)
estimated_Rt_national <- data.frame(Rt = colMeans(out_national$Rt),
                                    Rt_min_1 = colQuantiles(out_national$Rt,prob=0.025),
                                    Rt_max_1 = colQuantiles(out_national$Rt,prob=0.975),
                                    Rt_min_2 = colQuantiles(out_national$Rt,prob=0.25),
                                    Rt_max_2 = colQuantiles(out_national$Rt,prob=0.75),
                                    time= complete_dates)

data_Rt_95 <- data.frame(time = estimated_Rt_national$time, Rt_min = estimated_Rt_national$Rt_min_1, 
                         Rt_max = estimated_Rt_national$Rt_max_1, key = rep("nintyfive", length(estimated_Rt_national$time)))

data_Rt_50 <- data.frame(time = estimated_Rt_national$time, Rt_min = estimated_Rt_national$Rt_min_2, 
                         Rt_max = estimated_Rt_national$Rt_max_2, key = rep("fifty", length(estimated_Rt_national$time)))


quant_Rt_national <- rbind(data_Rt_95, data_Rt_50)
quant_Rt_national$key <- factor(quant_Rt_national$key, levels = c("ninetyfive", "fifty"))

Rt_threshold <- data.frame(time = data_estimated_rt$time, Rt = rep(1,length(data_estimated_rt$time)))
estimated_Rt_national$country <- "National"

p <- ggplot()+
  geom_ribbon(data = quant_Rt_national, aes(x=time, ymin = Rt_min, ymax = Rt_max, fill=key))+
  geom_line(data = melted_Rt_region,aes(x=time,y=Rt,color = region),linewidth =0.6)+
  geom_line(data = estimated_Rt_national, aes(x=time, y = Rt,linetype = "National"),color="midnightblue", linewidth = 1)+
  geom_line(data = Rt_threshold, aes(x=time, y = Rt),linewidth =0.6)+
  labs(x = "", y= " Estimated Rt")+
  scale_x_date(date_breaks = "4 weeks", labels = date_format("%d %b"))+
  scale_fill_manual(name = "", labels = c("50% CI of Rt_national", "95% CI of Rt_national"),
                    values = c(alpha("deepskyblue4", 0.35), 
                               alpha("deepskyblue3", 0.25))) + 
  scale_color_manual(name = "region", values = region_colors) +
  scale_linetype_manual(name = NULL, values = c("National" = "solid")) +
  ggtitle("Estimated Rt")+
  theme_pubr(base_family="sans") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),      # Increase legend title size
        legend.text = element_text(size = 12),       # Increase legend text size
        legend.key.size = unit(1.5, "cm")) +         # Increase legend key size
  guides(fill=guide_legend(ncol=1))

plot(p)

# ggsave(filename = paste0("figures/national_region_rt.png"), plot = p, width=12, height=8, units="in")


















