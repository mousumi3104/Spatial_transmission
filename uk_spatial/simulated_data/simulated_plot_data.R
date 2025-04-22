# this is the plot for simulated data. The estimated connected rt and estimated disconnected rt

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
library(latex2exp)

# rm(list = ls())
script_directory <- this.path::this.dir()
setwd(script_directory)

#-------- load data ----------------------------------------------------------------------
simulated_data <-readRDS("data/updated_si/simulated_infection.rds")
load("results/updated_si/connected_region_fitting.Rdata")
load("results/updated_si/disconnected_region_fitting.Rdata")

M_regions <- 3

#-------- arrangement to plot -----------------------------------------------------------------------------
#---------- true data -------------------------------------------------------------------------------------
true_infection <- data.frame(region1 = simulated_data$infection[1:(length(simulated_data$infection)/3)],
                             region2 = simulated_data$infection[((length(simulated_data$infection)/3)+1):(2*(length(simulated_data$infection)/3))],
                             region3 = simulated_data$infection[((2*length(simulated_data$infection)/3)+1):(3*(length(simulated_data$infection)/3))],
                             time = 1:(length(simulated_data$infection)/3))
true_infection <- true_infection %>% filter(time <= 350)
true_Rt <- simulated_data$Rt %>% mutate(time = 1:(length(simulated_data$infection)/3))
true_Rt <- true_Rt %>% filter(time <= 350)

# plot true Rt
colors_true_rt = c("Region 1"= "salmon", "Region 2" = "darkolivegreen3", "Region 3"= "gold2")

plot_true_rt <- ggplot(true_Rt)+
  geom_line(data = true_Rt, aes(x = time,y = Rt_1, color = "Region 1"), linewidth = 1.2)+
  geom_line(data = true_Rt, aes(x = time, y = Rt_2, color = "Region 2"), linewidth = 1.2)+
  geom_line(data = true_Rt, aes(x = time, y = Rt_3, color = "Region 3"), linewidth = 1.2)+
  geom_hline(yintercept = 1, color = "black", linewidth = 0.8)+
  # geom_vline(xintercept = 350, linetype = "dashed", color = "red")+
  xlab("Day")+
  ylab(expression(R[t]))+
  scale_color_manual(values = colors_true_rt)+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 0.5,size = 15,color="black",margin=margin(t=5)),
        axis.text.y = element_text(size = 15, margin = margin(r=5),color="black"),
        axis.title.y = element_text(size = 15, margin=margin(r=5)),
        axis.title.x = element_text(size = 15, margin=margin(t=0.5)),
        plot.title = element_text(size=15, margin = margin(l = 15,b=10),hjust = 0.5),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        legend.position = "bottom",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 15),       # Increase legend text size
        legend.key.size = unit(0.7, "cm"),
        legend.margin = margin(t = -10, r = 0, b = 0, l = 0),
        legend.justification = c(0, 0))+
  guides(fill=guide_legend(nrow=1))

print(plot_true_rt)

# ggsave(filename = paste0("figures/updated_si/simulated_original_rt.png"), plot = plot_true_rt, width = 4.55, height = 3, dpi = 300)


daily_inf <- true_infection
colors_inf = c("Region 1"= "salmon", "Region 2" = "darkolivegreen3", "Region 3"= "gold2")

plot_inf <- ggplot(daily_inf)+
  geom_point(data = daily_inf, aes(x = time,y = region1, color = "Region 1"), size = 1)+
  geom_point(data = daily_inf, aes(x = time, y = region2, color = "Region 2"), size = 1)+
  geom_point(data = daily_inf, aes(x = time, y = region3, color = "Region 3"), size = 1 )+
  
  xlab("Day")+
  ylab("Daily infection")+
  scale_color_manual(values = colors_inf)+
  ggtitle("")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 0.5,size = 15,color="black",margin=margin(t=5)),
        axis.text.y = element_text(size = 15, margin = margin(r=5),color="black"),
        axis.title.y = element_text(size = 15, margin=margin(r=5)),
        axis.title.x = element_text(size = 15, margin=margin(t=0.5)),
        plot.title = element_text(size=15, margin = margin(l = 15,b=10),hjust = 0.5),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        legend.position = "bottom",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 15),       # Increase legend text size
        # legend.key.size = unit(2, "cm"),
        legend.margin = margin(t = -10, r = 0, b = 0, l = 0),
        legend.justification = c(0, 0))+
  guides(fill=guide_legend(nrow=1))

print(plot_inf)

# ggsave(filename = paste0("figures/updated_si/simulated_infection.png"), plot = plot_inf, width = 4.55, height = 3, dpi = 300)

 
#---------- estimated RT disconnected model ------------------------------------------------------------------------------------
M_regions <- stan_data_disconnected$M_regions

for (i in 1:M_regions){
  final_time <- stan_data_disconnected$final_time
  est_Rt_disc <- fit_disconnected$draws("Rt",format="matrix")
  data_est_Rt_disc <- data.frame(est_Rt_disc_mean = colMeans(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)]),
                                 Rt_disc_min_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                 Rt_disc_max_1 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                 Rt_disc_min_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                 Rt_disc_max_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                 time = 1 : final_time)
  
 data_est_Rt_disc <- data_est_Rt_disc %>%
    filter(time >=stan_data_connected$fitting_start & time <= 350)

   data_Rt_disc_95 <- data.frame(time = data_est_Rt_disc$time, Rt_disc_min = data_est_Rt_disc$Rt_disc_min_1,
                            Rt_disc_max = data_est_Rt_disc$Rt_disc_max_1, key = rep("nintyfive", length(data_est_Rt_disc$time)))

   data_Rt_disc_50 <- data.frame(time = data_est_Rt_disc$time, Rt_disc_min = data_est_Rt_disc$Rt_disc_min_2,
                            Rt_disc_max = data_est_Rt_disc$Rt_disc_max_2, key = rep("fifty", length(data_est_Rt_disc$time)))

  data_Rt_disc <- data_Rt_disc_95
  data_Rt_disc$key1 <- "95% CI of estimated\ndisconnected Rt"

#------- connected model --------------------------------------------------------------------------------------------

  est_Rt_con <- fit_connected$draws("Rt",format="matrix")
  data_est_Rt_con <- data.frame(est_Rt_con_mean = colMeans(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)]),
                                Rt_con_min_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                Rt_con_max_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                Rt_con_min_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                Rt_con_max_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                time = 1 : final_time)
  
  data_est_Rt_con <- data_est_Rt_con %>%
    filter(time >= stan_data_connected$fitting_start & time <= 350)

 data_Rt_con_95 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_1,
                            Rt_con_max = data_est_Rt_con$Rt_con_max_1, key = rep("nintyfive", length(data_est_Rt_con$time)))

 data_Rt_con_50 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_2,
                            Rt_con_max = data_est_Rt_con$Rt_con_max_2, key = rep("fifty", length(data_est_Rt_con$time)))

#  data_Rt_con <- rbind(data_Rt_con_95, data_Rt_con_50)
#  levels(data_Rt_con$key) <- c("ninetyfive", "fifty")
 data_Rt_con <- data_Rt_con_95
 data_Rt_con$key1 <- "95% CI of estimated\nconnectded Rt"

 Rt_threshold <- data.frame(time = data_est_Rt_disc$time, Rt = rep(1,length(data_est_Rt_disc$time)))  # for Rt threshold horizontal line

#---------- plot ------------------------------------------------------------------------------------------------

 colors_rt <- c("Estimated \ndisconnected Rt" = "#e34a33", "Estimated \nconnected Rt" = "chartreuse3", "True Rt" = "blue4")
 true_Rt <- true_Rt %>% filter(time >= stan_data_connected$fitting_start & time <= 350)
 region = c("Region 1", "Region 2", "Region 3")
 
 
 plot_rt <- ggplot(data_est_Rt_disc) +
   geom_ribbon(data = data_Rt_disc, aes(x = time, ymin = Rt_disc_min, ymax = Rt_disc_max, fill = key1),show.legend = FALSE) +
   geom_ribbon(data = data_Rt_con, aes(x = time, ymin = Rt_con_min, ymax = Rt_con_max, fill = key1),show.legend = FALSE) +
   geom_line(data = data_est_Rt_disc, aes(x = time, y = est_Rt_disc_mean, color = "Estimated \ndisconnected Rt"), linewidth = 1.5) +
   geom_line(data = data_est_Rt_con, aes(x = time, y = est_Rt_con_mean, color = "Estimated \nconnected Rt"), linewidth = 1.5) +
   geom_line(data = true_Rt, aes(x = time, y = !!sym(paste0("Rt_",i)), color = "True Rt"), linewidth = 0.8) + # replace Rt_1 with the correct column name
   geom_line(data = Rt_threshold, aes(x = time, y = Rt), color = "black", linetype= "dashed", linewidth = 0.8) +
   xlab("Day") +
   ylab("") +
   scale_color_manual(
     values = colors_rt,
     labels = c(TeX(r"(\overset{\normalsize{Estimated}}{\overset{\normalsize{$connected~R_t$}}})"),
                TeX(r"(\overset{\normalsize{Estimated}}{\overset{\normalsize{$disconnected~R_t$}}})"),
                TeX(r"(\normalsize{True}{\normalsize{$~R_t$}})"))
   ) +
   scale_fill_manual(
     name = "",
     values = c("95% CI of estimated\nconnectded Rt" = alpha("chartreuse3", 0.15),
                "95% CI of estimated\ndisconnected Rt" = alpha("#b2182b", 0.15))
   ) +
   scale_y_continuous(limits = c(0.6,2.5))+
   # ggtitle(region[i])+
   theme_bw() +
   theme(
     axis.text.x = element_text(angle = 0, hjust = 0.4, vjust = 0.4, size = 25, margin = margin(r = 10), color = "black"),
     axis.text.y = element_text(size = 25, margin = margin(r = 10), color = "black"),
     axis.title.y = element_text(size = 25, margin = margin(r = 10)),
     axis.title.x = element_text(size = 25, margin = margin(r = 10)),
     plot.title = element_text(size = 25, margin = margin(l = 15, b = 10), hjust = 0.5),
     axis.ticks.x = element_line(colour = "black"),
     axis.ticks.length = unit(0.3, "cm"),
     axis.ticks = element_line(linewidth = 1),
     panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
     panel.grid.minor = element_blank(),
     panel.border = element_rect(color = "black", linewidth = 2),
     legend.position = "right",
     legend.title = element_blank(),
     legend.text = element_text(size = 25),
     legend.key.size = unit(2, "cm")
   ) +
   guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))  # Ensure the color legend is displayed
 
 # plot(plot_rt) 
 assign(paste0("rt",i),plot_rt)
} 

index = c("(a)","(b)","(c)")
# Remove legends from individual plots

legend_rt <- get_legend(plot_rt)  # Ensure legend is extracted

rt_list <-  list(rt1,rt2,rt3)

for (m in 1:length(rt_list)){
  rt_list[[m]] <- rt_list[[m]] + theme(legend.position = "none")
}

# final_rt <- ggdraw() +
#   draw_plot(plot_grid(plot_grid(plotlist = rt_list,  nrow = 1, ncol = 3, rel_widths = c(1,1,1), align = "hv",axis = "tblr",
#                                 labels = index,
#                                 label_size = 25,
#                                 label_x = 0.2,
#                                 label_y = 1.05),
#                       legend_rt, nrow=2, rel_heights = c(0.75,0.25)),x = 0.02, y = 0, width = 0.95, height = 0.98)+
#   draw_label(TeX(r"($R_t$)"),  x = 0.02, y = 0.6,  angle=90, size = 22)
# 
# final_rt
# 
# final_rt <- final_rt + theme(plot.background = element_rect(fill = "white", color = NA))
# 
# print(final_rt)
  
# ggsave(filename = paste0("figures/updated_si/simulated_data_est_rt.png"), plot = final_rt, width = 13, height = 5.5, dpi = 300)




#########################################################################################
#-------- plot original Rt (for presentation) -------------------------------------------
#########################################################################################

# load("data/simulated_data.Rdata")
# 
# threshold_Rt <- data.frame(threshold = rep(1,nrow(daily_infection_data)),
#                            time = 1:nrow(daily_infection_data))
# 
# p1 <- ggplot(Rt)+
#   geom_line(data = Rt, aes(x = index,y = Rt_1, color = "Region 1"), linewidth = 1.2)+
#   geom_line(data = threshold_Rt, aes(x=time,y=threshold),color = "black")+
#   xlab("")+
#   ylab("")+
#   ggtitle("Original Rt")+
#   theme_bw()+
#   scale_color_manual(values = c("Region 1" = "coral1")) +
#   theme(axis.text.x = element_blank(),#element_text( hjust = 0.5,size =15),
#         axis.text.y = element_text(size = 20,color = "black"),
#         axis.title.x = element_blank(),#element_text(size = 20, margin = margin(t=10)),
#         legend.position = "right",
#         legend.title = element_blank(),      # Increase legend title size
#         legend.text = element_text(size = 18),       # Increase legend text size
#         legend.key.size = unit(1.2, "cm"),
#         panel.grid = element_blank(),
#         plot.title = element_text(size = 24, hjust = 0))  # Increase title size and center it)
# 
# p2 <- ggplot(Rt)+
#   geom_line(data = Rt, aes(x = index,y = Rt_2, color = "Region 2"), linewidth = 1.2)+
#   geom_line(data = threshold_Rt, aes(x=time,y=threshold),color = "black")+
#   xlab("")+
#   ylab("")+
#   ggtitle("")+
#   theme_bw()+
#   scale_color_manual(values = c("Region 2" = "darkolivegreen3")) +
#   theme(axis.text.x = element_blank(),#element_text( hjust = 0.5,size =15),
#         axis.text.y = element_text(size = 20,color="black"),
#         axis.title.x =element_blank(),# element_text(size = 20, margin = margin(t=10)),
#         legend.position = "right",
#         legend.title = element_blank(),      # Increase legend title size
#         legend.text = element_text(size = 18),       # Increase legend text size
#         legend.key.size = unit(1.2, "cm"),
#         panel.grid = element_blank())  
# 
# p3 <- ggplot(Rt)+
#   geom_line(data = Rt, aes(x = index,y = Rt_3, color = "Region 3"), linewidth = 1.2)+
#   geom_line(data = threshold_Rt, aes(x=time,y=threshold),color = "black")+
#   xlab("Day")+
#   ylab("")+
#   ggtitle("")+
#   theme_bw()+
#   scale_color_manual(values = c("Region 3" = "darkgoldenrod2")) +
#   theme(axis.text.x = element_text( hjust = 0.5,size =18, color = "black"),
#         axis.text.y = element_text(size = 20,color = "black"),
#         axis.title.x = element_text(size = 20, margin = margin(t=10)),
#         legend.position = "right",
#         legend.title = element_blank(),      # Increase legend title size
#         legend.text = element_text(size = 18),       # Increase legend text size
#         legend.key.size = unit(1.2, "cm"),
#         panel.grid = element_blank())
# 
# p1 <- p1 + theme(panel.spacing = unit(0, "lines"))
# p2 <- p2 + theme(panel.spacing = unit(0, "lines"))
# 
# 
# p <- plot_grid(p1,NULL,p2,NULL,p3, ncol = 1,rel_heights = c(0.7,0.001,0.67,0.001,0.9))
# 
#   
# p_final <- ggdraw(p) +
#           # draw_label("Your Title Here", x = 0.5, y = 0.95, size = 20, hjust = 0.5, vjust = 1) +  # Add title
#           draw_label(expression(R[t]), x = 0, y = 0.5, angle = 90, vjust = 1, size = 20) # Add common y-axis label
#           # draw_plot(p, x = 0, y = 0, width = 1, height = 0.9)  # Draw the plot below the title
# 
# 
# p_final
# 
########################################################################################
# #------ original incidence ------------------------------------------------------------
#########################################################################################

# colors_incidence <- c("Region 1" = "coral1", "Region 2" = "darkolivegreen3", "Region 3"="darkgoldenrod2")
# 
# plot_inf <-ggplot(daily_infection_data)+
#   
#   geom_point(data = daily_infection_data, aes(x = time,y = region1, color = "Region 1"), size = 2)+
#   geom_point(data = daily_infection_data, aes(x = time, y = region2, color = "Region 2"), size = 2)+
#   geom_point(data = daily_infection_data, aes(x = time,y = region3,color = "Region 3"), size = 2)+
#   
#   xlab("Day")+
#   ylab("Incidence")+
#   scale_color_manual(values = colors_incidence)+
#   scale_shape_manual(values = 16)+
#   ggtitle(NULL)+
#   theme_bw()+
#   theme(axis.text.x = element_text( hjust = 0.5,size =20,color = "black"),
#         axis.text.y = element_text(size = 20,color= "black"),
#         axis.title.x = element_text(size = 20, margin = margin(t=10)),
#         axis.title.y = element_text(size = 20, margin = margin(r=10)),
#         legend.position = "right",
#         legend.title = element_blank(),      # Increase legend title size
#         legend.text = element_text(size = 18),       # Increase legend text size
#         legend.key.size = unit(1.2, "cm"))  +
#   
#   guides(fill=guide_legend(ncol=1))
# print(plot_inf)
# 


#############################################################################
#--- plot mobility ---------------------------------------------------------
#############################################################################


# regions <- c("Region1","Region2","Region3")
# mobility <- matrix(c(0.7, 0.15, 0.15, 0.07, 0.85, 0.08, 0.21, 0.09, 0.7), nrow=3, ncol=3)
# mobility_df <- melt(mobility)
# colnames(mobility_df) <- c("X", "Y", "Value")
# mobility_df$X <- factor(mobility_df$X, labels = regions)
# mobility_df$Y <- factor(mobility_df$Y, labels = regions)
# ggplot(mobility_df, aes(x = X, y = Y, fill = Value)) +
#   geom_tile() +  # Create the heatmap
#   geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 5) +  # Add values inside the tiles
#   scale_fill_gradient(low = "white", high = "steelblue") +  # Color gradient for the values
#   theme_minimal() +
#   labs(x = "", y = "", fill = "Value") +  # Axis and legend labels
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  # Rotate x-axis labels and set size
#         axis.text.y = element_text(size = 12))  # Set y-axis text size
# 
