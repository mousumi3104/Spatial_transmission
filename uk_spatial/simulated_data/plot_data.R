library(cmdstanr)
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
library(ggpubr)
library(matrixStats)
library(cowplot)
#library(rstanarm)  
library(this.path)

rm(list = ls())
script_directory <- this.path::this.dir()
setwd(script_directory)

#-------- load data ----------------------------------------------------------------------
load("data/simulated_data.Rdata")
load("data/connected_region_fitting.Rdata")
load("data/fitting_national.Rdata")

M_regions <- 3
final_time <- 70 * 7
#-------- arrangement to plot -----------------------------------------------------------------------------

#---------- true data -------------------------------------------------------------------------------------
true_infection <- daily_infection_data
true_infection$index <- 1:nrow(true_infection)
true_Rt <- Rt

  
#---------- disconnected model ------------------------------------------------------------------------------------
for (i in 1:M_regions){
  
  load("data/disconnected_region_fitting.Rdata")
  fit <- fit_disconnected
  est_Rt_disc <- fit$draws("Rt",format="matrix")
  data_est_Rt_disc <- data.frame(est_Rt_disc_mean = colMeans(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)]),
                                 Rt_disc_min_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                 Rt_disc_max_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                 Rt_disc_min_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                 Rt_disc_max_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                 time = 1 : final_time)

   data_Rt_disc_95 <- data.frame(time = data_est_Rt_disc$time, Rt_disc_min = data_est_Rt_disc$Rt_disc_min_1,
                            Rt_disc_max = data_est_Rt_disc$Rt_disc_max_1, key = rep("nintyfive", length(data_est_Rt_disc$time)))

   data_Rt_disc_50 <- data.frame(time = data_est_Rt_disc$time, Rt_disc_min = data_est_Rt_disc$Rt_disc_min_2,
                            Rt_disc_max = data_est_Rt_disc$Rt_disc_max_2, key = rep("fifty", length(data_est_Rt_disc$time)))

  data_Rt_disc <- data_Rt_disc_95
  data_Rt_disc$key1 <- "95% CI of estimated\ndisconnectded Rt"
#---------------------------

  est_inf_disc <- fit$draws("infection",format="matrix")       #fit_disconnected$draws("infection", format = "matrix")
  data_est_inf_disc <- data.frame(est_inf_disc_mean = colMeans(est_inf_disc[,(((i-1)*final_time)+1):(i*final_time)]),
                                   inf_disc_min_1 = colQuantiles(est_inf_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                   inf_disc_max_1 = colQuantiles(est_inf_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                   inf_disc_min_2 = colQuantiles(est_inf_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                   inf_disc_max_2 = colQuantiles(est_inf_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                   time = 1 : final_time)

  data_inf_disc_95 <- data.frame(time = data_est_inf_disc$time, inf_disc_min = data_est_inf_disc$inf_disc_min_1,
                            inf_disc_max = data_est_inf_disc$inf_disc_max_1, key = rep("nintyfive", length(data_est_inf_disc$time)))

  data_inf_disc_50 <- data.frame(time = data_est_inf_disc$time, inf_disc_min = data_est_inf_disc$inf_disc_min_2,
                            inf_disc_max = data_est_inf_disc$inf_disc_max_2, key = rep("fifty", length(data_est_inf_disc$time)))

  # data_inf_disc <- rbind(data_inf_disc_95, data_inf_disc_50)
  # levels(data_inf_disc$key) <- c("ninetyfive", "fifty")
  data_inf_disc <- data_inf_disc_95
  data_inf_disc$key1 <- "95% CI of estimated\ndisconnectded incidence"

#------- connected model --------------------------------------------------------------------------------------------

  fit <- fit_connected
  est_Rt_con <- fit$draws("Rt",format="matrix")
  data_est_Rt_con <- data.frame(est_Rt_con_mean = colMeans(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)]),
                                Rt_con_min_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                Rt_con_max_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                Rt_con_min_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                Rt_con_max_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                time = 1 : final_time)

 data_Rt_con_95 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_1,
                            Rt_con_max = data_est_Rt_con$Rt_con_max_1, key = rep("nintyfive", length(data_est_Rt_con$time)))

 data_Rt_con_50 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_2,
                            Rt_con_max = data_est_Rt_con$Rt_con_max_2, key = rep("fifty", length(data_est_Rt_con$time)))

#  data_Rt_con <- rbind(data_Rt_con_95, data_Rt_con_50)
#  levels(data_Rt_con$key) <- c("ninetyfive", "fifty")
 data_Rt_con <- data_Rt_con_95
 data_Rt_con$key1 <- "95% CI of estimated\nconnectded Rt"
#-------------------------------

  est_inf_con <- fit$draws("infection",format="matrix")
  data_est_inf_con <- data.frame(est_inf_con_mean = colMeans(est_inf_con[,(((i-1)*final_time)+1):(i*final_time)]),
                                 inf_con_min_1 = colQuantiles(est_inf_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                 inf_con_max_1 = colQuantiles(est_inf_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                 inf_con_min_2 = colQuantiles(est_inf_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                 inf_con_max_2 = colQuantiles(est_inf_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                 time = 1:final_time)

  data_inf_con_95 <- data.frame(time = data_est_inf_con$time, inf_con_min = data_est_inf_con$inf_con_min_1,
                            inf_con_max = data_est_inf_con$inf_con_max_1, key = rep("nintyfive", length(data_est_inf_con$time)))

  data_inf_con_50 <- data.frame(time = data_est_inf_con$time, inf_con_min = data_est_inf_con$inf_con_min_2,
                            inf_con_max = data_est_inf_con$inf_con_max_2, key = rep("fifty", length(data_est_inf_con$time)))

#  data_inf_con <- rbind(data_inf_95, data_inf_50)
#  levels(data_inf_con$key) <- c("ninetyfive", "fifty")
 data_inf_con <- data_inf_con_95
 data_inf_con$key1 <- "95% CI of estimated\nconnectded incidence"

 Rt_threshold <- data.frame(time = data_est_Rt_disc$time, Rt = rep(1,length(data_est_Rt_disc$time)))  # for Rt threshold horizontal line

#---------- plot ------------------------------------------------------------------------------------------------

colors_rt <- c("Estimated disconnected Rt" = "red4", "Estimated connected Rt" = "green4", "Simulated Rt"="black")
colors_incidence <- c("Estimated disconnected\nincidence" = "red4", "Estimated connected\nincidence" = "green4", "Simulated\nincidence"="coral3")
 
plot_rt <- ggplot(data_est_Rt_disc)+
  
  geom_ribbon(data = data_Rt_disc, aes(x = time, ymin = Rt_disc_min, ymax = Rt_disc_max, fill=key1))+
  geom_ribbon(data = data_Rt_con, aes(x = time, ymin = Rt_con_min, ymax = Rt_con_max, fill=key1))+
  geom_line(data = data_est_Rt_disc, aes(x = time,y = est_Rt_disc_mean, color = "Estimated disconnected Rt"), linewidth = 1.4)+
  geom_line(data = data_est_Rt_con, aes(x = time, y = est_Rt_con_mean, color = "Estimated connected Rt"), linewidth = 1.5)+
  geom_line(data = true_Rt, aes(x = index, y = !!sym(paste0("Rt_",i)), color = "Simulated Rt"), linewidth = 1 )+
  geom_line(data = Rt_threshold, aes(x=time, y = Rt),color = "black")+
  xlab(NULL)+
  ylab("")+
  scale_fill_manual(name = "",
                    values = c("95% CI of estimated\ndisconnectded Rt" = alpha("red4", 0.25),
                               "95% CI of estimated\nconnectded Rt" = alpha("green4", 0.25))) +
  scale_color_manual(values = colors_rt)+
  ggtitle(paste("Region",i))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 0.5,size = 15),
        axis.text.y = element_text(size = 15,margin = margin(r=10)),
        axis.title.y = element_text(size = 20, margin=margin(r=10)),
        plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 14),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"))+
  guides(fill=guide_legend(ncol=1))
 
 if (i == 1){
   plot_rt <- plot_rt + ylab(expression(R[t]))+
     theme(axis.title.y = element_text(size = 20, margin = margin(r=15)))
 }else {
   plot_rt <- plot_rt + theme(axis.title.y = element_blank())
 }

# plot(plot_rt)
assign(paste0("rt",i),plot_rt)

plot_inf <-ggplot(data_est_inf_disc)+
  
  geom_ribbon(data = data_inf_disc, aes(x = time, ymin = inf_disc_min, ymax = inf_disc_max, fill=key1))+
  geom_ribbon(data = data_inf_con, aes(x = time, ymin = inf_con_min, ymax = inf_con_max, fill=key1))+
  geom_line(data = data_est_inf_disc, aes(x = time,y = est_inf_disc_mean, color = "Estimated disconnected\nincidence"), linewidth = 1.4)+
  geom_line(data = data_est_inf_con, aes(x = time, y = est_inf_con_mean, color = "Estimated connected\nincidence"), linewidth = 1.4)+
  geom_point(data = true_infection, aes(x = index,y = !!sym(paste0("region",i)),color = "Simulated\nincidence"))+

  xlab("Day")+
  ylab("")+
  scale_fill_manual(name = "",
                    values = c("95% CI of estimated\ndisconnectded incidence" = alpha("red4", 0.25),
                               "95% CI of estimated\nconnectded incidence" = alpha("seagreen3", 0.25))) +
  scale_color_manual(values = colors_incidence)+
  scale_shape_manual(values = 16)+
  ggtitle(NULL)+
  theme_bw()+
  theme(axis.text.x = element_text( hjust = 0.5,size =15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(t=10)),
        legend.position = "right",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 14),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"))  +
  
  guides(fill=guide_legend(ncol=1))
if (i == 1){
  plot_inf <- plot_inf + ylab("Incidence")+
    theme(axis.title.y = element_text(size = 20, margin = margin(r=10)))
}else {
  plot_inf <- plot_inf + theme(axis.title.y = element_blank())
}

assign(paste0("inf",i),plot_inf)

}

legend_rt <- get_legend(plot_rt)
legend_inf <- get_legend(plot_inf) 

# Remove legends from individual plots
rt1 <- rt1 + theme(legend.position = "none")
rt2 <- rt2 + theme(legend.position = "none")
rt3 <- rt3 + theme(legend.position = "none")

inf1 <- inf1 + theme(legend.position = "none")
inf2 <- inf2 + theme(legend.position = "none")
inf3 <- inf3 + theme(legend.position = "none")

p_top <- plot_grid(
  NULL, rt1, NULL, rt2, NULL, rt3, legend_rt, 
  nrow = 1, rel_widths = c(0.055, 1.07, 0.06, 0.968, 0.068, 0.96, 0.625),labels = c("","(a)","","(b)","","(c)",""),
  label_size = 25,label_x = c(0.01,0.1,0,0.01,0.01,0.01,0.01),label_fontface = "plain")

p_bottom <- plot_grid(
  inf1, inf2, inf3, legend_inf, 
  nrow = 1, rel_widths = c(1.08, 1, 1, 0.6), labels = c("(d)","(e)","(f)",""),label_fontface = "plain",label_size = 25,
  label_x = c(0.14,0.1,0.1,0.01),label_y = c(1.08,1.08,1.08,1))

p <- plot_grid(p_top, NULL, p_bottom, ncol = 1, rel_heights = c(1, 0.08, 1))
p <- p + theme(plot.background = element_rect(fill = "white", color = NA))
print(p)

ggsave(filename = paste0("figures/simulated_data_est_rt.png"), plot = p, width = 18, height = 10, dpi = 300)

