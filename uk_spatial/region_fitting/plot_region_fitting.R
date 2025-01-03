# this is the plot for estimated rt based on england COVID dataset

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
# load("results/double_mobility/region_connected_rt_double_mob_xyz.Rdata")
load("results/region_disconnected_rt_including_jan.Rdata")
load("results/region_connected_rt_including_jan.Rdata")

inf_start_date <- plot_required_date$inf_start_date
fitting_start <- plot_required_date$fitting_start_date
end_date <- plot_required_date$end_date
plot_end <- as.Date("2020-12-31", format = "%Y-%m-%d")
final_time <- stan_data_connected$final_time
first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")

second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 
Rt <- fit_connected$draws("Rt",format = "matrix")
load("data/final_pop_2020_ltla.Rdata")

M_regions <- stan_data_connected$M_regions
death_regions <- data.frame(stan_data_disconnected$death)
death_data_length <- stan_data_disconnected$death_data_length
death_regions$time <- seq(from = inf_start_date ,to =  end_date, by = "week")
# death_regions <- death_regions %>%
#   filter(time >= fitting_start)

regions <- c("North East","North West","Yorkshire and the Humber","East Midlands","West Midlands","East","London","South East","South West")

#---------- disconnected model ------------------------------------------------------------------------------------
for (i in 1:M_regions){
  
  # load(paste0("results/disconnected/region_disconnected_xyz",i,".rds"))
  fit <- fit_disconnected
  est_Rt_disc <- fit$draws("Rt",format="matrix")
  data_est_Rt_disc <- data.frame(est_Rt_disc_mean = colMeans(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)]),
                                Rt_disc_min_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                Rt_disc_max_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                Rt_disc_min_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                Rt_disc_max_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                time = seq(from=inf_start_date ,to = end_date, by = "day"))
  data_est_Rt_disc <- data_est_Rt_disc %>% filter(time >= fitting_start & time <= plot_end)
  
  
  # 
  # est_Rt_disc <- fit$draws("Rt",format="matrix")
  # data_est_Rt_disc <- data.frame(est_Rt_disc_mean = colMeans(est_Rt_disc),
  #                                Rt_disc_min_1 = colQuantiles(est_Rt_disc,prob=0.025),
  #                                Rt_disc_max_1 = colQuantiles(est_Rt_disc,prob=0.975),
  #                                Rt_disc_min_2 = colQuantiles(est_Rt_disc,prob=0.25),
  #                                Rt_disc_max_2 = colQuantiles(est_Rt_disc,prob=0.75),
  #                                time = seq(from=inf_start_date ,to =  end_date, by = "day"))
  # data_est_Rt_disc <- data_est_Rt_disc %>%
  #   filter(time >= fitting_start)
  
  data_Rt_disc_95 <- data.frame(time = data_est_Rt_disc$time, Rt_disc_min = data_est_Rt_disc$Rt_disc_min_1,
                                Rt_disc_max = data_est_Rt_disc$Rt_disc_max_1, key = rep("nintyfive", length(data_est_Rt_disc$time)))
  
  data_Rt_disc_50 <- data.frame(time = data_est_Rt_disc$time, Rt_disc_min = data_est_Rt_disc$Rt_disc_min_2,
                                Rt_disc_max = data_est_Rt_disc$Rt_disc_max_2, key = rep("fifty", length(data_est_Rt_disc$time)))
  
  data_Rt_disc <- data_Rt_disc_95
  data_Rt_disc$key1 <- "95% CI of \ndisconnected Rt"
  #------- connected model --------------------------------------------------------------------------------------------
  
  fit <- fit_connected
  est_Rt_con <- fit$draws("Rt",format="matrix")
  data_est_Rt_con <- data.frame(est_Rt_con_mean = colMeans(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)]),
                                Rt_con_min_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                Rt_con_max_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                Rt_con_min_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                Rt_con_max_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                time = seq(from=inf_start_date ,to =  end_date, by = "day"))
  data_est_Rt_con <- data_est_Rt_con %>% filter(time >= fitting_start & time <= plot_end)
  
  data_Rt_con_95 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_1,
                               Rt_con_max = data_est_Rt_con$Rt_con_max_1, key = rep("nintyfive", length(data_est_Rt_con$time)))
  
  data_Rt_con_50 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_2,
                               Rt_con_max = data_est_Rt_con$Rt_con_max_2, key = rep("fifty", length(data_est_Rt_con$time)))
  
  #  data_Rt_con <- rbind(data_Rt_con_95, data_Rt_con_50)
  #  levels(data_Rt_con$key) <- c("ninetyfive", "fifty")
  data_Rt_con <- data_Rt_con_95
  data_Rt_con$key1 <- "95% CI of \nconnected Rt"

  Rt_threshold <- data.frame(time = data_est_Rt_disc$time, Rt = rep(1,length(data_est_Rt_disc$time)))  # for Rt threshold horizontal line
  
  #---------- plot ------------------------------------------------------------------------------------------------
  
  breaks = sort(c(seq(ymd("2020-4-1"),ymd("2020-12-31"),by="months"),seq(ymd("2020-3-15"),ymd("2020-12-31"),by="months"),ymd("2021-1-1")))
  
  labels = unique(date_format("%b")(data_est_Rt_con$time))
  labels = as.vector(rbind(labels,rep("",length(labels))))
  colors_rt <- c("Connected Rt" = "#1a9850", "Disconnected Rt" = "#b2182b")#, "Simulated Rt"="black")

  plot_rt <- ggplot(data_est_Rt_disc)+
    
    geom_ribbon(data = data_Rt_disc, aes(x = time, ymin = Rt_disc_min, ymax = Rt_disc_max, fill=key1))+
    geom_ribbon(data = data_Rt_con, aes(x = time, ymin = Rt_con_min, ymax = Rt_con_max, fill=key1))+
    geom_line(data = data_est_Rt_disc, aes(x = time,y = est_Rt_disc_mean, color = "Disconnected Rt"), linewidth = 1)+
    geom_line(data = data_est_Rt_con, aes(x = time, y = est_Rt_con_mean, color = "Connected Rt"), linewidth = 1)+
    geom_line(data = Rt_threshold, aes(x=time, y = Rt),color = "black")+
    geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 1)+
    xlab("")+
    ylab("")+
    scale_fill_manual(name = "",
                      values = c("95% CI of \nconnected Rt" = alpha("#1a9850", 0.25),
                                 "95% CI of \ndisconnected Rt" = alpha("#b2182b", 0.25))) +
    scale_color_manual(values = colors_rt)+
    scale_x_date(labels = labels, breaks = breaks, limits = c(fitting_start,ymd("2020-12-31"))) + 
    ggtitle(regions[i])+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 15,margin = margin(r=10),color="black"),
          axis.text.y = element_text(size = 17,margin = margin(r=10),color="black"),
          axis.title.y = element_text(size = 20, margin=margin(r=10)),
          axis.title.x = element_text(size = 20, margin=margin(r=10)),
          plot.title = element_text(size=18, margin = margin(l = 15,b=10),hjust = 0.2),
          axis.ticks.x= element_line(colour=c(rep(c(NA,"black"), t=9),NA)),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks = element_line(linewidth =1),
          panel.grid.major.x = element_line(colour=c(rep(c(NA,"grey94"), t=9),NA)),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),      # Increase legend title size
          legend.text = element_text(size = 20),       # Increase legend text size
          legend.key.size = unit(1.2, "cm"))+
    guides(fill=guide_legend(nrow=1))

  # plot(plot_rt)
  assign(paste0("rt",i),plot_rt)
}

legend_rt <- get_legend(plot_rt)
plot_rt_list <-  list(rt1,rt2,rt3,rt4,rt5,rt6,rt7,rt8,rt9)

for (i in c(1,4,7)){
  plot_rt_list <- list(get(paste0("rt",i)),get(paste0("rt",i+1)),get(paste0("rt",i+2)))
  for (m in 1:length(plot_rt_list)){
    plot_rt_list[[m]] <- plot_rt_list[[m]] + theme(legend.position = "none")
  }
  three_rt <- plot_grid(plotlist = plot_rt_list,  nrow = 1, ncol = 3,rel_widths = c(1,1,1), align = "hv",axis = "tblr")
  
  three_rt <- ggdraw() +
    draw_plot(three_rt, x = 0.01, y = 0, width = 0.95, height = 1) +  
    draw_label(expression(R[t]),  x = 0.02, y = 0.6,  angle=90, size = 22)   # y-label
  
  assign(paste0("three_rt", (i+2)/3),three_rt) 
}

final_rt_list <- list(get(paste0("three_rt",1)),get(paste0("three_rt",2)),get(paste0("three_rt",3)))
for (m in 1:length(final_rt_list)){
  final_rt_list[[m]] <- final_rt_list[[m]] + theme(legend.position = "none")
}

final_rt <- plot_grid(plot_grid(plotlist = final_rt_list,  nrow = 3, ncol = 1,rel_heights = c(1,1,1), align = "hv",axis = "tblr"),
                      legend_rt, nrow=2, rel_heights = c(2.5,0.15))
final_rt <- final_rt + theme(plot.background = element_rect(fill = "white", color = NA))

# p_rt <- ggdraw() +
#   draw_label(expression("Daily"~"estimated"~R[t]~"over regions"), fontface = 'bold', x = 0.5, y = 0.98, hjust = 0.5, size = 20) +  # Title
#   draw_plot(p_rt, y = 0, height = 0.95)  # Add the combined plot
print(final_rt)

ggsave(filename = paste0("figures/estimated_rt_region_jan.png"), plot = final_rt, width=17, height=13, units="in")
#---------------

# ggsave(filename = paste0("figures/estimated_rt.png"), plot = p_rt, width=13, height=10, units="in")
# ggsave(filename = paste0("figures/estimated_data.png"), plot = p_inf, width=12, height=8, units="in")


#---------- fitted death data -------------------------------------------------------------------------
for (i in 1:M_regions){
  
  fit <- fit_disconnected
  
  est_death_disc <- fit$draws("weekly_deaths",format="matrix")
  data_est_death_disc <- data.frame(est_death_disc_mean = colMeans(est_death_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)]),
                                 death_disc_min_1 = colQuantiles(est_death_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.025),
                                 death_disc_max_1 = colQuantiles(est_death_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.975),
                                 death_disc_min_2 = colQuantiles(est_death_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.25),
                                 death_disc_max_2 = colQuantiles(est_death_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.75),
                                 time = seq(from=inf_start_date ,to =  end_date, by = "week"))
  
  data_est_death_disc <- data_est_death_disc %>% filter(time >= fitting_start & time <= plot_end)

  data_death_disc_95 <- data.frame(time = data_est_death_disc$time, death_disc_min = data_est_death_disc$death_disc_min_1,
                                 death_disc_max = data_est_death_disc$death_disc_max_1, key = rep("nintyfive", length(data_est_death_disc$time)))

  # data_death_disc_50 <- data.frame(time = data_est_death_disc$time, death_disc_min = data_est_death_disc$death_disc_min_2,
  #                                death_disc_max = data_est_death_disc$death_disc_max_2, key = rep("fifty", length(data_est_death_disc$time)))

  # data_death_disc <- rbind(data_death_disc_95, data_death_disc_50)
  # levels(data_death_disc$key) <- c("ninetyfive", "fifty")
  
  data_death_disc <- data_death_disc_95
  data_death_disc$key1 <- "95% CI of fitted death\n for disconnected model"
  # 
  #------------------------------------
  fit <- fit_connected
  est_death_con <- fit$draws("weekly_deaths",format="matrix")
  data_est_death_con <- data.frame(est_death_con_mean = colMeans(est_death_con[,(((i-1)*death_data_length)+1):(i*death_data_length)]),
                                 death_con_min_1 = colQuantiles(est_death_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.025),
                                 death_con_max_1 = colQuantiles(est_death_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.975),
                                 death_con_min_2 = colQuantiles(est_death_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.25),
                                 death_con_max_2 = colQuantiles(est_death_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.75),
                                 time = seq(from=inf_start_date ,to =  end_date, by = "week"))

  data_est_death_con <- data_est_death_con %>% filter(time >= fitting_start & time <= plot_end)

  data_death_con_95 <- data.frame(time = data_est_death_con$time, death_con_min = data_est_death_con$death_con_min_1,
                                death_con_max = data_est_death_con$death_con_max_1, key = rep("nintyfive", length(data_est_death_con$time)))

  # data_death_con_50 <- data.frame(time = data_est_death_con$time, death_con_min = data_est_death_con$death_con_min_2,
  #                               death_con_max = data_est_death_con$death_con_max_2, key = rep("fifty", length(data_est_death_con$time)))

  #  data_inf_con <- rbind(data_inf_95, data_inf_50)
  #  levels(data_inf_con$key) <- c("ninetyfive", "fifty")
  
  data_death_con <- data_death_con_95
  data_death_con$key1 <- "95% CI of fitted death\n for connected model"
  
  death_regions$column_to_plot <- death_regions[[i]]
  death_regions <- death_regions %>% filter(time >= fitting_start & time <= plot_end)

  colors_death <- c("Fitted deaths for \ndicconected model" = "red4", "Fitted deaths for \nconnected model" = "green4", "Weekly death \ndata"="coral3")
  
  breaks = sort(c(seq(ymd("2020-4-1"),ymd("2020-12-31"),by="months"),seq(ymd("2020-3-15"),ymd("2020-12-31"),by="months"),ymd("2021-1-1")))
  
  labels = unique(date_format("%b")(data_est_Rt_con$time))
  labels = as.vector(rbind(labels,rep("",length(labels))))
  
  fit_death <-ggplot(data_est_death_disc)+
                geom_ribbon(data = data_death_disc, aes(x = time, ymin = death_disc_min, ymax = death_disc_max, fill=key1))+
                geom_ribbon(data = data_death_con, aes(x = time, ymin = death_con_min, ymax = death_con_max, fill=key1))+
                geom_line(data = data_est_death_disc, aes(x = time,y = est_death_disc_mean, color = "Fitted deaths for \ndicconected model"), linewidth = 1)+
                geom_line(data = data_est_death_con, aes(x = time, y = est_death_con_mean, color = "Fitted deaths for \nconnected model"), linewidth = 1)+
                geom_point(data = death_regions, aes(x = time, y = column_to_plot, color = "Weekly death \ndata")) +
                xlab("")+
                ylab("")+
                scale_fill_manual(name = "", values = c("95% CI of fitted death\n for disconnected model" = alpha("red4", 0.25),
                               "95% CI of fitted death\n for connected model" = alpha("seagreen3", 0.25))) +
                scale_color_manual(values = colors_death)+
                scale_shape_manual(values = 16)+
                ggtitle(regions[i])+
                scale_x_date(labels = labels, breaks = breaks, limits = c(fitting_start, ymd("2020-12-31"))) +
                theme_bw()+
                theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 15,margin = margin(r=10),color="black"),
                      axis.text.y = element_text(size = 17,margin = margin(r=10),color="black"),
                      axis.ticks.x= element_line(colour=c(rep(c(NA,"black"), t=9),NA)),
                      axis.ticks.length = unit(0.3,"cm"),
                      axis.ticks = element_line(linewidth =1),
                      panel.grid.major.x = element_line(colour=c(rep(c(NA,"grey94"), t=9),NA)),
                      panel.grid.minor = element_blank(),
                      axis.title.x = element_text(size = 20, margin = margin(t=10)),
                      plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
                      legend.position = "bottom",
                      legend.title = element_blank(),      # Increase legend title size
                      legend.text = element_text(size = 20),       # Increase legend text size
                      legend.key.size = unit(1.2, "cm"))  +
                      guides(fill=guide_legend(nrow=1))
  # plot(fit_death)
  assign(paste0("death",i),fit_death)

}
legend_death <- get_legend(fit_death)
plot_death_list <-  list(death1,death2,death3,death4,death5,death6,death7,death8,death9)

for (i in c(1,4,7)){
  plot_death_list <- list(get(paste0("death",i)),get(paste0("death",i+1)),get(paste0("death",i+2)))
  for (m in 1:length(plot_death_list)){
    plot_death_list[[m]] <- plot_death_list[[m]] + theme(legend.position = "none")
  }
  three_deaths <- plot_grid(plotlist = plot_death_list,  nrow = 1, ncol = 3,rel_widths = c(1,1,1), align = "hv",axis = "tblr")
  
  three_deaths <- ggdraw() +
    draw_plot(three_deaths, x = 0.01, y = 0, width = 0.95, height = 1) +  
    draw_label("Weekly deaths",  x = 0.015, y = 0.6,  angle=90, size = 17)   # y-label
  
  assign(paste0("three_deaths", (i+2)/3),three_deaths) 
}

final_death_list <- list(get(paste0("three_deaths",1)),get(paste0("three_deaths",2)),get(paste0("three_deaths",3)))
for (m in 1:length(final_death_list)){
  final_death_list[[m]] <- final_death_list[[m]] + theme(legend.position = "none")
}

final_death <- plot_grid(plot_grid(plotlist = final_death_list,  nrow = 3, ncol = 1,rel_heights = c(1,1,1), align = "hv",axis = "tblr"),
                      legend_death, nrow=2, rel_heights = c(2.5,0.15))
final_death <- final_death + theme(plot.background = element_rect(fill = "white", color = NA))

print(final_death)

ggsave(filename = paste0("figures/fitted_region_death.png"), plot = final_death, width=17, height=12, units="in")

#------------------------------------------------------------------------------------------------------
#----- plot original weekly death data ----------------------------------------------------------------
#------------------------------------------------------------------------------------------------------


death_data <- read_excel("data/death_20_21.xlsx")  
load("data/final_pop_2020_ltla.Rdata") 

#-------- death data arrangements -------------------------------------------------------
pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

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

death_data <- death_data %>% select(all_of(pop_2020$area_name))  

death_regions <- data.frame('North East' = apply(death_data[,north_east_index],1,sum),
                            'North West' = apply(death_data[,north_west_index],1,sum),
                            "Yorkshire" = apply(death_data[,yorkshire_index],1,sum),
                            'East Midlands' = apply(death_data[,east_midlands_index],1,sum),
                            'West Midlands' = apply(death_data[,west_midlands_index],1,sum),
                            "East" = apply(death_data[,east_index],1,sum),
                            "London" = apply(death_data[,london_index],1,sum),
                            'South East' = apply(death_data[,south_east_index],1,sum),
                            'South West' = apply(death_data[,south_west_index],1,sum),
                            Week = seq(ymd(20200101),ymd(20210131),by="week"))

first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")
second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 
#---- facet plot for all the death data -------------------------------------------------------------
death_region_long <-
  death_regions %>%
  pivot_longer(cols = !Week,names_to = "region", values_to = "deaths") %>%
  mutate('Weekly deaths' = as.integer(deaths)) 

death_region_long$region <- factor(death_region_long$region, 
                                   levels = c("North.East","North.West",
                                              "Yorkshire","East.Midlands",
                                              "West.Midlands","East","London","South.East","South.West"))

breaks = sort(c(seq(ymd("2020-1-1"),ymd("2021-1-1"),by="months"),seq(ymd("2020-1-15"),ymd("2021-1-31"),by="months"),ymd("2021-02-1")))

labels = unique(date_format("%b")(death_regions$Week))
labels = as.vector(rbind(labels,rep("",length(labels))))
labels = c("",labels,"Jan","")

p <- ggplot(data = death_region_long, aes(Week,deaths)) +
  geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 0.8)+
  geom_line(color = "darkblue", linewidth = 1)+
  scale_x_date(labels=labels,breaks=breaks,limits = c(ymd("20200101"),ymd("20210201"))) + 
  facet_grid(rows = vars(region))+facet_wrap(~region)+
  ylab("Weekly deaths")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 15,margin = margin(r=10),color="black"),
        axis.text.y = element_text(size = 15,margin = margin(r=10),color="black"),
        axis.ticks.x= element_line(colour=c(rep(c("black",NA), t=13),"black")),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major.x = element_line(colour=c(rep(c("white",NA), t=13),"white")),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 20, margin=margin(r=10)),
        axis.title.x = element_text(size = 20, margin=margin(r=10)),
        plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
        strip.text = element_text(size = 20, face = "bold"),
        legend.position = "",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 15),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.y = unit(10, "cm"))+
  guides(fill=guide_legend())

plot(p)
ggsave(filename = paste0("figures/death_data.png"), plot = p, width=17, height=13, units="in")








# ------------ mobility matrix------------------------------------------------------------------------

mobility <- stan_data_connected$C_base
mobility_df <- melt(mobility)
colnames(mobility_df) <- c("X", "Y", "Value")
mobility_df$X <- factor(mobility_df$X, labels = regions)
mobility_df$Y <- factor(mobility_df$Y, labels = regions)
ggplot(mobility_df, aes(x = Y, y = X, fill = Value)) +
  geom_tile() +  # Create the heatmap
  geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 6,fontface="bold") +  # Add values inside the tiles
  scale_fill_gradient(low = "#fff", high = "#5ab4ac") +  # Color gradient for the values
  theme_minimal() +
  labs(title = "Mobility matrix", x = "", y = "", fill = "Value") +  # Axis and legend labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 17,color="black"),  # Rotate x-axis labels and set size
        axis.text.y = element_text(size = 17,color="black"),
        plot.title = element_text(size=18, margin = margin(l = 15,b=10),face = "bold"),
        legend.position = "none")  # Set y-axis text size


#--- for national model-------
# est_inf_disc <- fit$draws("weekly_deaths",format="matrix")       #fit_disconnected$draws("infection", format = "matrix")
# data_est_inf_disc <- data.frame(est_inf_disc_mean = colMeans(est_inf_disc),
#                                 inf_disc_min_1 = colQuantiles(est_inf_disc,prob=0.025),
#                                 inf_disc_max_1 = colQuantiles(est_inf_disc,prob=0.975),
#                                 inf_disc_min_2 = colQuantiles(est_inf_disc,prob=0.25),
#                                 inf_disc_max_2 = colQuantiles(est_inf_disc,prob=0.75),
#                                 time = seq(from=inf_start_date ,to =  end_date, by = "week"))

# data_est_inf_disc <- data_est_inf_disc %>% filter(time >= fitting_start)