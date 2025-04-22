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
library(latex2exp)

# library(pathchwork)

script_directory <- this.path::this.dir()
setwd(script_directory)

#----------------------------------- data ------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# load("results/double_mobility/region_connected_rt_double_mob_xyz.Rdata")
load("results/updated_si/gm_region_disconnected_rt_jan.Rdata")
load("results/updated_si/gm_region_connected_rt_jan.Rdata")

inf_start_date <- plot_required_date$inf_start_date
fitting_start <- plot_required_date$fitting_start_date
end_date <- plot_required_date$end_date

#------------------------------------------------------------------
plot_end <- as.Date("2021-12-31", format = "%Y-%m-%d")
plot_start <- fitting_start

final_time <- stan_data_connected$final_time
first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")

second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 
Rt <- fit_connected$draws("Rt",format = "matrix")
load("data/final_pop_2020_ltla.Rdata")

M_regions <- stan_data_connected$M_regions
death_regions <- data.frame(stan_data_connected$death)
death_data_length <- stan_data_connected$death_data_length
death_regions$time <- seq(from = inf_start_date ,to =  end_date, by = "week")
# death_regions <- death_regions %>%
#   filter(time >= fitting_start)

regions <- c("North East","North West","Yorkshire \nand the Humber","East Midlands","West Midlands","East","London","South East","South West")

#---------- disconnected model ------------------------------------------------------------------------------------
for (i in 1:M_regions){
  
  # load(paste0("results/disconnected/region_disconnected_xyz",i,".rds"))
  fit <- fit_disconnected
  est_Rt_disc <- fit$draws("Rt",format="matrix")
  data_est_Rt_disc <- data.frame(est_Rt_disc_mean = colMeans(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)]),
                                Rt_disc_min_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.05),
                                Rt_disc_max_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.95),
                                Rt_disc_min_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.05),
                                Rt_disc_max_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.95),
                                time = seq(from=inf_start_date ,to = end_date, by = "day"))
  data_est_Rt_disc <- data_est_Rt_disc %>% filter(time >= plot_start & time <= plot_end)
  
  
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
                                Rt_con_min_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.05),
                                Rt_con_max_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.95),
                                Rt_con_min_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.05),
                                Rt_con_max_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.95),
                                time = seq(from=inf_start_date ,to =  end_date, by = "day"))
  data_est_Rt_con <- data_est_Rt_con %>% filter(time >= plot_start & time <= plot_end)
  
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
  
  breaks = sort(c(seq(ymd("2020-3-1"),ymd("2021-1-31"),by="months"),seq(ymd("2020-3-15"),ymd("2021-1-31"),by="months"), ymd("2021-1-31")))
  
  labels = unique(date_format("%b")(breaks))
  labels = c("", as.vector(rbind(labels,rep("",length(labels)))))
  colors_rt <- c("Connected Rt" = "#1a9850", "Disconnected Rt" = "#b2182b")#, "Simulated Rt"="black")
  
  

  plot_rt <- ggplot(data_est_Rt_disc)+
    
    geom_ribbon(data = data_Rt_disc, aes(x = time, ymin = Rt_disc_min, ymax = Rt_disc_max, fill=key1),show.legend = FALSE)+
    geom_ribbon(data = data_Rt_con, aes(x = time, ymin = Rt_con_min, ymax = Rt_con_max, fill=key1),show.legend = FALSE)+
    geom_line(data = data_est_Rt_disc, aes(x = time,y = est_Rt_disc_mean, color = "Disconnected Rt"), linewidth = 1.3)+
    geom_line(data = data_est_Rt_con, aes(x = time, y = est_Rt_con_mean, color = "Connected Rt"), linewidth = 1.3)+
    geom_line(data = Rt_threshold, aes(x=time, y = Rt),color = "black", linetype = "dashed", linewidth = 0.8)+
    geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dotted", color = "grey30", linewidth = 1.2)+
    xlab("")+
    ylab("")+
    scale_fill_manual(name = "",
                      values = c("95% CI of \nconnected Rt" = alpha("#1a9850", 0.2),
                                 "95% CI of \ndisconnected Rt" = alpha("#b2182b", 0.2)),
                      labels = c(TeX(r"(\overset{\normalsize{95\% CI of estimated}}{\overset{\normalsize{$connected~R_t$}}})"),
                                 TeX(r"(\overset{\normalsize{95\% CI of estimated}}{\overset{\normalsize{$disconnected~R_t$}}})"))
    )+
    scale_color_manual(values = colors_rt,
                       labels = c(TeX(r"($Connected~R_t$)"),
                                  TeX(r"($Disconnected~R_t$)")))+
    scale_x_date(labels = labels, breaks = breaks, limits = c(fitting_start,ymd("2020-12-31"))) + 
    scale_y_continuous(labels = c(0,"",2,"",4,"",6),breaks = c(0,1,2,3,4,5,6), limits = c(0,6))+
    ggtitle(regions[i])+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust = 0.4, vjust = 0.4,size = 20,margin = margin(r=10),color="black"),
          axis.text.y = element_text(size = 20,margin = margin(r=10),color="black"),
          axis.title.y = element_text(size = 20, margin=margin(r=10)),
          axis.title.x = element_text(size = 20, margin=margin(r=10)),
          plot.title = element_text(size=23, margin = margin(l=15,b=10),hjust = 0.5),
          plot.margin = margin(t=0,r=0,b=0,l=0),
          axis.ticks.x= element_line(colour=c(rep(c("black",NA), t=floor(length(labels)/2)),NA)),
          axis.ticks.y= element_line(colour=c(rep(c("black",NA), t=floor(length(labels)/2)),NA)),
          axis.ticks.length = unit(0.2,"cm"),
          axis.ticks = element_line(linewidth =1),
          panel.grid.major.x = element_line(colour=c(rep(c("grey98",NA), t=floor(length(labels)/2)),NA)),
          # panel.grid.major.y = element_line(colour=c(rep("grey90", t=length(labels)*2))),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth=1.5),
          legend.position = "bottom",
          legend.title = element_blank(),      # Increase legend title size
          legend.text = element_text(size = 25),       # Increase legend text size
          legend.key.size = unit(1.2, "cm"))+
    guides(fill=guide_legend(nrow=1))
  # plot_rt <- ggdraw() +
  #   draw_plot(plot_rt, x = 0, y = 0, width = 1, height = 0.95) +  
  #   draw_label(index[i],  x = 0.15, y = 0.95,  angle=0, size = 23)   # y-label
  # plot(plot_rt)
  assign(paste0("rt",i),plot_rt)
}

index = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)")
legend_rt <- get_legend(plot_rt)
rt_list <-  list(rt1,rt2,rt3,rt4,rt5,rt6,rt7,rt8,rt9)

for (m in 1:length(rt_list)){
  rt_list[[m]] <- rt_list[[m]] + theme(legend.position = "none")
}

# 
# p_mob
# blankplot1
# blankplot2
# 
# plot_list = c(list(p_mob, blankplot1, blankplot1),rt_list)

final_rt <- ggdraw() +
            draw_plot(plot_grid(plot_grid(plotlist = rt_list,  nrow = 3, ncol = 3,rel_heights = c(1,1,1), rel_widths = c(1,1,1), align = "hv",axis = "tblr",
                                labels = index,
                                label_size = 22,
                                label_x = 0.12,
                                label_y = 1.01),
                      legend_rt, nrow=2, rel_heights = c(2.5,0.15)),x = 0.02, y = 0, width = 0.95, height = 1)+
            draw_label(TeX(r"($Estimated~R_t$)"),  x = 0.02, y = 0.6,  angle=90, size = 22)

final_rt

final_rt <- final_rt + theme(plot.background = element_rect(fill = "white", color = NA))
# ggsave(filename = paste0("figures/updated_si/estimated_rt_region_jan.png"), plot = final_rt, width=15, height=13, units="in", dpi=300)


#---------------

# ggsave(filename = paste0("figures/estimated_rt.png"), plot = p_rt, width=13, height=10, units="in")
# ggsave(filename = paste0("figures/estimated_data.png"), plot = p_inf, width=12, height=8, units="in")


#---------- fitted death data -------------------------------------------------------------------------
end_date <- as.Date("2021-1-31", format = "%Y-%m-%d")
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
  data_death_disc$key1 <- "95% CI of fitted death\n with disconnected model"
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
  data_death_con$key1 <- "95% CI of fitted death\n with connected model"
  
  death_regions$column_to_plot <- death_regions[[i]]
  death_regions <- death_regions %>% filter(time >= fitting_start & time <= end_date)

  colors_death <- c("Fitted deaths with \ndisconnected model" = "red4", "Fitted deaths with \nconnected model" = "green4", "Observed death \ndata"="coral3")
  
  breaks = sort(c(seq(ymd("2020-3-1"),ymd("2021-1-31"),by="months"),seq(ymd("2020-3-15"),ymd("2021-1-31"),by="months"),ymd("2021-1-31")))
  
  labels = unique(date_format("%b")(breaks))
  labels = as.vector(c("",rbind(labels,rep("",length(labels)))))
  
  fit_death <-ggplot(data_est_death_disc)+
                geom_ribbon(data = data_death_disc, aes(x = time, ymin = death_disc_min, ymax = death_disc_max, fill=key1))+
                geom_ribbon(data = data_death_con, aes(x = time, ymin = death_con_min, ymax = death_con_max, fill=key1))+
                geom_line(data = data_est_death_disc, aes(x = time,y = est_death_disc_mean, color = "Fitted deaths with \ndisconnected model"), linewidth = 1.3)+
                geom_line(data = data_est_death_con, aes(x = time, y = est_death_con_mean, color = "Fitted deaths with \nconnected model"), linewidth = 1.3)+
                geom_point(data = death_regions, aes(x = time, y = column_to_plot, color = "Observed death \ndata"), size = 2) +
                geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dotted", color = "grey30", linewidth = 1.2)+
                xlab("")+
                ylab("")+
                scale_fill_manual(name = "", values = c("95% CI of fitted death\n with disconnected model" = alpha("red4", 0.25),
                               "95% CI of fitted death\n with connected model" = alpha("seagreen3", 0.25))) +
                scale_color_manual(values = colors_death)+
                
                # scale_shape_manual(values = 16)+
                ggtitle(regions[i])+
                scale_x_date(labels = labels, breaks = breaks, limits = c(fitting_start, ymd("2021-1-31"))) +
                theme_bw()+
                theme(axis.text.x = element_text(angle = 90,hjust = 0.4, vjust = 0.4,size = 15,margin = margin(r=10),color="black"),
                      axis.text.y = element_text(size = 17,margin = margin(r=10),color="black"),
                      axis.ticks.x= element_line(colour=c(rep(c("black",NA), t=floor(length(labels)/2)),NA)),
                      axis.ticks.length = unit(0.3,"cm"),
                      axis.ticks = element_line(linewidth =1),
                      panel.grid.major.x = element_line(colour=c(rep(c("grey93",NA), t=floor(length(labels)/2)),NA)),
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
index = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)")

for (m in 1:length(plot_death_list)){
  plot_death_list[[m]] <- plot_death_list[[m]] + theme(legend.position = "none")
}

final_fitting_death <- ggdraw() +
  draw_plot(plot_grid(plot_grid(plotlist = plot_death_list,  nrow = 3, ncol = 3,rel_heights = c(1,1,1), rel_widths = c(1,1,1), align = "hv",axis = "tblr",
                                labels = index,
                                label_size = 22,
                                label_x = 0.14,
                                label_y = 1.01),
                      legend_death, nrow=2, rel_heights = c(2.5,0.15)),x = 0.02, y = 0, width = 0.95, height = 1)+
  draw_label("Weekly deaths",  x = 0.02, y = 0.6,  angle=90, size = 22)

final_fitting_death <- final_fitting_death + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(filename = paste0("figures/updated_si/fitted_region_death_jan.png"), plot = final_fitting_death, width=16, height=13, units="in")



#-----------------------------------------------------------------------------------------------------
#------ plot uncertainty of the estimated parameters --------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

load("results/region_disconnected_rt_including_jan.Rdata")
load("results/region_connected_rt_including_jan.Rdata")
regions <- c("North East","North West","Yorkshire and the Humber","East Midlands","West Midlands","East","London","South East","South West")

type = c("connected", "disconnected")
for (j in type){
parameters <-  c("initial_seeding","mu","x","y","z")
stan_fitting <- get(paste0("fit_",j))
fit_summary <- stan_fitting$summary()
for (i in parameters){
  parameter_data <- fit_summary %>% filter(str_starts(variable, i))
  parameter_data$region <- regions
  parameter_data$region <- factor(parameter_data$region, levels = parameter_data$region)
  p <-ggplot(data = parameter_data, aes(x = region, y = mean)) +
    geom_pointrange(data = parameter_data, aes(ymin = q5, ymax = q95), size = 0.5,color= if(j == "connected"){"green4"} else {"red4"})+
    labs(
      x ="",
      y = paste("Estimated",i))+
    coord_flip()+
    theme_bw()+
    scale_y_continuous( limits = c(min(parameter_data$q5),max(parameter_data$q95)))+
    theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 10,color="black"),
          axis.text.y = element_text(size = 10,margin = margin(r=10),color="black"),
          axis.title.y = element_text(size = 15, margin=margin(r=-10)),
          axis.title.x = element_text(size = 15, margin=margin(t=12)),
          plot.title = element_text(size=10, margin = margin(l = 15,b=10),hjust = 0.5),
          axis.ticks.length = unit(0.2,"cm"),
          axis.ticks=element_line(linewidth =0.7),
          panel.grid.major.x=element_line(colour = "grey90"),
          panel.grid.minor.x=element_blank(),
          legend.position = "",
          legend.title = element_blank(),      # Increase legend title size
          legend.text = element_blank())      # Increase legend text size)
  
  plot(p)
  ggsave(filename = paste0("figures/",i,"_england_region_",j,".png"), plot = p, width=5, height=5, units="in")
  
}
}

#------------------------------------------------------------------------------------------------------
#----- plot original weekly death data ----------------------------------------------------------------
#------------------------------------------------------------------------------------------------------


death_data <- read_excel("data/death_20_21.xlsx")
load("data/final_pop_2020_ltla.Rdata") 
#------- regions index -----------------------------------------------------------------

pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

pop_2020$area_name <- sapply(pop_2020$area_name, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

colnames(death_data) <- sapply(colnames(death_data), function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

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

death_regions <- data.frame("North East" = apply(death_data[,north_east_index],1,sum),
                            "North West" = apply(death_data[,north_west_index],1,sum),
                            "Yorkshire and the Humber" = apply(death_data[,yorkshire_index],1,sum),
                            "East Midlands" = apply(death_data[,east_midlands_index],1,sum),
                            "West Midlands" = apply(death_data[,west_midlands_index],1,sum),
                            "East" = apply(death_data[,east_index],1,sum),
                            "London" = apply(death_data[,london_index],1,sum),
                            "South East" = apply(death_data[,south_east_index],1,sum),
                            "South West" = apply(death_data[,south_west_index],1,sum),
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

death_region_long$region <- gsub("\\.", " ", death_region_long$region)
death_region_long$region <- factor(death_region_long$region,levels = c("North East", "North West",
                                                                       "Yorkshire and the Humber", "East Midlands",
                                                                       "West Midlands", "East", "London", 
                                                                       "South East", "South West"))

breaks = sort(c(seq(ymd("2020-1-1"),ymd("2021-1-1"),by="months"),seq(ymd("2020-1-15"),ymd("2021-1-31"),by="months"),ymd("2021-02-1")))

labels = unique(date_format("%b")(death_regions$Week))
labels = as.vector(rbind(labels,rep("",length(labels))))
labels = c("",labels,"Jan","")

p <- ggplot(data = death_region_long, aes(Week,deaths)) +
  geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 0.8)+
  geom_bar(stat = "identity", fill = "darkblue", color = "white", linewidth = 0.1)+
  scale_x_date(labels=labels,breaks=breaks,limits = c(ymd("20200305"),ymd("20210201"))) + 
  facet_grid(rows = vars(region))+facet_wrap(~region, nrow=4, ncol=3)+
  ylab("Weekly deaths")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.4, vjust = 0.4,size = 19,margin = margin(r=10),color="black"),
        axis.text.y = element_text(size = 20,margin = margin(r=10),color="black"),
        axis.ticks.x= element_line(colour=c(rep(c("black",NA), t=13),"black")),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major.x = element_line(colour=c(rep(c("#f0f0f0",NA), t=13),"#f0f0f0")),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 25, margin=margin(r=10)),
        axis.title.x = element_text(size = 15, margin=margin(r=10)),
        plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
        strip.text = element_text(size = 25),
        panel.background = element_rect(fill = "white", color = "black", linewidth=1),
        legend.position = "",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 15),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.y = unit(10, "cm"))+
  guides(fill=guide_legend())

plot(p)
ggsave(filename = paste0("figures/death_data_england.png"), plot = p, width=15, height=10, units="in",dpi=300)


death_region_long$region <- factor(death_region_long$region, 
                                   levels = c("North East", "North West", 
                                              "Yorkshire and the Humber", "East Midlands", 
                                              "West Midlands", "East", "London", 
                                              "South East", "South West"))

labels <- format(breaks, "%b")  # Format month labels

p <- ggplot(data = death_region_long, aes(Week, deaths)) +
  geom_vline(xintercept = as.Date(c(first_lockdown_start, first_lockdown_end, 
                                    second_lockdown_start, second_lockdown_end)),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_bar(stat = "identity", fill = "darkblue", color = "white", linewidth = 0.1) +
  scale_x_date(labels = labels, breaks = breaks, limits = c(ymd("20200305"), ymd("20210201"))) +
  facet_wrap(~region, nrow = 3, ncol = 3) +  # Use facet_wrap OR facet_grid
  ylab("Weekly deaths") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")

plot(p)




####--------------------------------------------------------------------------------------------------
# ------------ mobility matrix------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

regions <- c("North East","North West","Yorkshire \nand the Humber","East Midlands","West Midlands","East","London","South East","South West")

mobility <- stan_data_connected$C_base
mobility_df <- melt(mobility)
colnames(mobility_df) <- c("X", "Y", "Value")
mobility_df$X <- factor(mobility_df$X, labels = regions)
mobility_df$Y <- factor(mobility_df$Y, labels = regions)
p_mob <- ggplot(mobility_df, aes(x = Y, y = X, fill = Value)) +
  geom_tile(color = "lightgrey", linewidth = 0.3) +  # Create the heatmap
  geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 10) +  # Add values inside the tiles,fontface="bold"
  scale_fill_gradient(low = "#fff", high = "#5ab4ac") +  # Color gradient for the values
  theme_minimal() +
  labs(title = "", x = "", y = "", fill = "Value") +  # Axis and legend labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5,size = 25,color="black"),  # Rotate x-axis labels and set size
        axis.text.y = element_text(size = 25,color="black"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        #plot.title = element_text(size=18),#, margin = margin(l = 15,b=10)),
        legend.position = "none")  # Set y-axis text size

ggsave(filename = paste0("figures/england_mobility.png"), plot = p_mob, width=5, height=5, units="in", dpi =300)

# joint two mobility matrix 

p_mob_modified <- plot_grid(NULL, p_mob, NULL,  nrow = 1, rel_widths = c(0.05, 1,0.1))
p_joint <- plot_grid(NULL, p_mob_modified, p, ncol = 1, rel_heights = c(0.01,0.78, 1)) +
  draw_label("(a)", x = 0.20, y = 0.992, size = 30, fontface = "bold", hjust = 0) +
  draw_label("(b)", x = 0.20, y = 0.57, size = 30, fontface = "bold", hjust = 0)



ggsave(filename = paste0("figures/mobility_matrix.png"), plot = p_joint, width=19, height=19, units="in", dpi =300)



#-----------------------------------------------------------------------------------------------------
#--------- attaching together ------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

mobility_matrix <- plot_grid(plotlist = list(p_mob,p),  nrow = 1, ncol = 2,rel_widths = c(1,1.25), align = "hv",axis = "tblr",
                    labels = c("(a)","(b)"),
                    label_size = 15,
                    label_x = 0.25,label_y = 1.005)
mobility_matrix <- mobility_matrix + theme(plot.background = element_rect(fill = "white", color = NA))
print(mobility_matrix)
ggsave(filename = paste0("figures/mobility_matrix.png"), plot = mobility_matrix, width=12, height=5.5, units="in", dpi =300)


