# the estimated infections from three different sources ( Fig 3 (d-i))
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
library(cowplot)
library(ISOweek)
library(matrixStats)
library(ggpubr)
library(this.path)
library(ggplotify)
library(latex2exp)

# rm(list = ls())
script_directory <- this.path::this.dir()
setwd(script_directory)

# this is only for Rt and infection plot 
simulated_data <- readRDS("data/updated_si/simulated_infection.rds")
load("results/updated_si/connected_region_fitting.Rdata")
pop <- c(20000000,10000000,15000000)

fit <- fit_connected
init_R <- fit$draws("Rt",format="matrix")
Rt_connected <- fit$draws("Rt",format ="matrix")  # need to arrange
inf <- fit$draws("infection", format = "matrix")  # need to arrange

final_time <- stan_data_connected$final_time
M_regions <- stan_data_connected$M_regions 

# ------------- stan data arrange-------------------------------------------------------------
no_sample <- 300
infection <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_in_own <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_in_mob <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_out_mob <- array(data = NA, dim = c(final_time*M_regions, no_sample))

m <- cmdstan_model("simulated_region.stan")  
ind <- sample(1:800,no_sample)

for (k in 1:no_sample){ 
    stan_data <- list(M = stan_data_connected$M_regions,
                      pop = stan_data_connected$pop,
                      final_time = stan_data_connected$final_time,
                      initial_seeding_day = stan_data_connected$initial_seeding_day,
                      init_seed = as.vector(rep(inf[1,1],M_regions)),
                      SI = stan_data_connected$SI,
                      # f1 = stan_data_connected$f1,
                      # f2 = stan_data_connected$f2,
                      Rt = matrix(Rt_connected[ind[k],],final_time,M_regions),
                      C = stan_data_connected$C
                      )
    
    simulated_data <- m$sample(data =stan_data,
                               iter_sampling = 1,
                               chains = 1,
                               thin = 1, 
                               fixed_param = TRUE)  
    
    infection[,k] <- as.matrix(simulated_data$draws("infection"))
    infection_in_own[,k] <- as.matrix(simulated_data$draws("infection_in_own"))
    infection_in_mob[,k] <- as.matrix(simulated_data$draws("infection_in_mob"))
    infection_out_mob[,k] <- as.matrix(simulated_data$draws("infection_out_mob"))
    # weekly_deaths[,k] <- as.matrix(simulated_data$draws("weekly_deaths"))
}


time = 1:final_time

infection_data <- data.frame(inf_mean = rowMeans(infection),
                             inf_min1 = rowQuantiles(infection,prob = 0.05),
                             inf_max1 = rowQuantiles(infection,prob = 0.95))

inf_in_own_data <- data.frame(inf_mean = rowMeans(infection_in_own),
                             inf_min1 = rowQuantiles(infection_in_own,prob = 0.05),
                             inf_max1 = rowQuantiles(infection_in_own,prob = 0.95))

inf_in_mob_data <- data.frame(inf_mean = rowMeans(infection_in_mob),
                             inf_min1 = rowQuantiles(infection_in_mob,prob = 0.05),
                             inf_max1 = rowQuantiles(infection_in_mob,prob = 0.95))

inf_out_mob_data <- data.frame(inf_mean = rowMeans(infection_out_mob),
                             inf_min1 = rowQuantiles(infection_out_mob,prob = 0.05),
                             inf_max1 = rowQuantiles(infection_out_mob,prob = 0.95),time)

# death_data <- data.frame(death_mean = rowMeans(weekly_deaths),
                         # death_min1 = rowQuantiles(weekly_deaths, prob = 0.05),
                         # death_max1 = rowQuantiles(weekly_deaths, prob = 0.95))


plot_inf_list <- list()
plot_prop_list <- list()

for (m in 1:M_regions){

  plot_est_infection <- data.frame(time = stan_data_connected$fitting_start:final_time, 
                               est_infection_mean = colMeans(inf[,(((m-1)*final_time)+stan_data_connected$fitting_start):(m*final_time)]),
                               est_infection_min1 = colQuantiles(inf[,(((m-1)*final_time)+stan_data_connected$fitting_start):(m*final_time)], prob=0.01),
                               est_infection_max1 = colQuantiles(inf[,(((m-1)*final_time)+stan_data_connected$fitting_start):(m*final_time)], prob=0.99))
  
  plot_est_infection <- plot_est_infection %>% filter(time <=350)
  
  plot_data_infection <- data.frame(time = stan_data_connected$fitting_start:final_time, 
                              data_infection = stan_data_connected$data_inf[stan_data_connected$fitting_start:final_time,m])
  
  plot_data_infection <- plot_data_infection %>% filter(time <= 350)
  
  data_inf_95 <- data.frame(time = plot_est_infection$time, inf_min = plot_est_infection$est_infection_min1,
                               inf_max = plot_est_infection$est_infection_max1, key = rep("95% CI of total infection", length(plot_est_infection$time)))
  
  plot_inf_in_own <- inf_in_own_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_own$time <- 1:final_time
  plot_inf_in_own <- plot_inf_in_own %>% filter(time <=350)
  
  data_inf_in_own_95 <- data.frame(time = plot_inf_in_own$time, inf_min = plot_inf_in_own$inf_min1,
                            inf_max = plot_inf_in_own$inf_max1, key = rep("95% CI of infection in own", length(plot_inf_in_own$time)))
  
  plot_inf_in_mob <- inf_in_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_mob$time <- 1:final_time
  plot_inf_in_mob <- plot_inf_in_mob %>% filter(time <=350)
  
  data_inf_in_mob_95 <- data.frame(time = plot_inf_in_mob$time, inf_min = plot_inf_in_mob$inf_min1,
                                   inf_max = plot_inf_in_mob$inf_max1, key = rep("95% CI of infection in mob", length(plot_inf_in_mob$time)))
  
  plot_inf_out_mob <- inf_out_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_out_mob$time <- 1:final_time
  plot_inf_out_mob <- plot_inf_out_mob %>% filter(time <=350)
  
  data_inf_out_mob_95 <- data.frame(time = plot_inf_out_mob$time, inf_min = plot_inf_out_mob$inf_min1,
                               inf_max = plot_inf_out_mob$inf_max1, key = rep("95% CI of infection out mob", length(plot_inf_out_mob$time)))
  
  data_inf <- rbind(data_inf_95, data_inf_in_own_95, data_inf_in_mob_95, data_inf_out_mob_95)
  data_inf$key <- factor(data_inf$key, levels = c("95% CI of total infection", "95% CI of infection in own", "95% CI of infection in mob","95% CI of infection out mob"))
  
data_stack_plot <-  data.frame(in_mob = plot_inf_in_mob$inf_mean, in_own = plot_inf_in_own$inf_mean, out_mob = plot_inf_out_mob$inf_mean, time = plot_inf_in_mob$time)

data_stack_plot1 <- data.frame(
  time = rep(data_stack_plot$time, times = 3),
  types = factor(rep(c("in_own", "in_mob", "out_mob"), each = length(data_stack_plot$time)),levels = c("in_own", "in_mob", "out_mob")),
  infections = c(data_stack_plot$in_own, data_stack_plot$in_mob, data_stack_plot$out_mob)
)

p  <- ggplot(data_stack_plot1, aes(x = time, y = infections, fill = types)) +
   geom_area(position = "stack") +
   geom_ribbon(data = data_inf_out_mob_95, aes(x = time, ymin = inf_min, ymax = inf_max), fill = "#7570b3", alpha = 0.5, inherit.aes = FALSE) +
   geom_ribbon(data = data_inf_in_mob_95, aes(x = time, ymin = plot_inf_out_mob$inf_mean  + inf_min, ymax = plot_inf_out_mob$inf_mean + inf_max), fill = "deepskyblue4", alpha = 0.5, inherit.aes = FALSE) +
   geom_ribbon(data = data_inf_in_own_95, aes(x = time, ymin = plot_inf_out_mob$inf_mean  + plot_inf_in_mob$inf_mean + inf_min, ymax = plot_inf_out_mob$inf_mean  + plot_inf_in_mob$inf_mean + inf_max), fill = "#ff7f00", alpha = 0.5, inherit.aes = FALSE) +
   geom_line(data = plot_inf_out_mob, aes(x = time, y = inf_mean),color = "#7570b3",inherit.aes = FALSE, linewidth = 1)+
   geom_line(data = plot_inf_in_mob, aes(x = time, y = plot_inf_out_mob$inf_mean  + inf_mean),color = "deepskyblue4",inherit.aes = FALSE,linewidth = 1)+
   geom_line(data = plot_inf_in_own, aes(x = time, y = plot_inf_out_mob$inf_mean  + plot_inf_in_mob$inf_mean +inf_mean),color = "#ff7f00",inherit.aes = FALSE, linewidth = 1)+
  scale_fill_manual(values = c("in_own" = alpha("#ff7f00",0.3), "in_mob" = alpha("deepskyblue4",0.3), "out_mob" = alpha("#7570b3",0.3)),
                    labels = c(TeX(r"(\overset{\normalsize{Infections driven by}}{\normalsize{own infected population}})"),
                               TeX(r"(\overset{\normalsize{Mobility induced}}{\normalsize{infections within region}})") ,
                               TeX(r"(\overset{\normalsize{Mobility induced infections}}{\normalsize{outside the region}})"))) +
    labs(
     x = "",
     y = ""
   ) +
   theme_bw() + scale_y_continuous(limits = c(0, 4000))+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 25,margin = margin(r=10),color="black"),
        axis.text.y = element_text(size = 25,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 25, margin=margin(r=10)),
        axis.title.x = element_text(size = 25, margin=margin(r=10)),
        plot.title = element_text(size=25, margin = margin(l = 15,b=10),hjust = 0.2),
        axis.ticks.x= element_line(colour="black"),#c(rep(c(NA,"black"), t=floor(length(labels)/2)),NA)),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 2),
        legend.position = "right",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 25),       # Increase legend text size
        legend.key.size = unit(2, "cm"))+
  guides(fill=guide_legend(nrow=1))

legend_inf <- get_legend(p) 
p<- p+theme(legend.position = "none")
print(p)
plot_inf_list[[m]] <- ggplotGrob(p)


#  -------plot_proportion -----------------------------------------------------------------------------------

data_stack_plot$in_mob_prop = data_stack_plot$in_mob / (data_stack_plot$in_mob + data_stack_plot$in_own + data_stack_plot$out_mob)
data_stack_plot$in_own_prop = data_stack_plot$in_own / (data_stack_plot$in_mob + data_stack_plot$in_own + data_stack_plot$out_mob)
data_stack_plot$out_mob_prop = data_stack_plot$out_mob / (data_stack_plot$in_mob + data_stack_plot$in_own + data_stack_plot$out_mob)

data_stack_plot <- lapply(data_stack_plot, function(df) {
  df[is.nan(df)] <- 0
  return(df)
})

data_stack_plot1 <- 
  data.frame(
    time = rep(data_stack_plot$time, times = 3),
    types = factor(rep(c("in_own_prop", "in_mob_prop", "out_mob_prop"), each = length(data_stack_plot$time)),levels = c("in_own_prop", "in_mob_prop", "out_mob_prop")),
    infections = c(data_stack_plot$in_own_prop, data_stack_plot$in_mob_prop, data_stack_plot$out_mob_prop)
  )

p_prop  <- ggplot(data_stack_plot1, aes(x = time, y = infections, fill = types)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = c("in_own_prop" = alpha("#ff7f00",0.45), "in_mob_prop" = alpha("deepskyblue4",0.45), "out_mob_prop" = alpha("#7570b3",0.45))) +
  labs(
    x = "Day",
    y = ""
  ) +
  theme_bw() + scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 25,margin = margin(r=10),color="black"),
        axis.text.y = element_text(size = 25,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 25, margin=margin(r=10)),
        axis.title.x = element_text(size = 25, margin=margin(r=10)),
        plot.title = element_text(size=25, margin = margin(l = 15,b=10),hjust = 0.2),
        axis.ticks.x= element_line(colour="black"),#c(rep(c(NA,"black"), t=floor(length(labels)/2)),NA)),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 2),
        legend.position = "right",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 25),       # Increase legend text size
        legend.key.size = unit(2, "cm"))+
  guides(fill=guide_legend(nrow=1))
p_prop<- p_prop+theme(legend.position = "none")
print(p_prop)
plot_prop_list[[m]] <- ggplotGrob(p_prop)
}

#--------------------------------------------------------------------------------------------------------------------
source("simulated_plot_data.R")

for (i in 1:3) {plot_inf_list[[i]] <- as.ggplot(plot_inf_list[[i]])}
for (i in 1:3) {plot_prop_list[[i]] <- as.ggplot(plot_prop_list[[i]])}
index = c("","(a)","","(b)","","(c)","","(d)","","(e)","","(f)","(g)","(h)","(i)")
labels_rt <- index[1:6]
labels_inf <- index[(6+1):(6+6)]
labels_prop <- index[(12+1):(12+3)]

# Arrange the first row (Rt plots)
p_rt <- plot_grid(
  plot_grid(plotlist = list(NULL,rt_list[[1]],NULL, rt_list[[2]],NULL,rt_list[[3]]), nrow = 1, ncol = 6, rel_widths = c(0.08,1,0.08,1,0.08,1), align = "hv",
            labels = labels_rt, label_size = 25, label_x = 0.24, label_y = 1.14,label_fontface = "plain"),
  legend_rt, nrow = 2, rel_heights = c(0.6, 0.2)  # Adjust legend height
)

# p_rt <- ggdraw()+
#   draw_plot(p_rt, x = 0.025, y = 0, width = 0.97, height = 0.93)
plot(p_rt)

# Arrange the second row (Infection plots)
p_inf <- plot_grid(plotlist = list(NULL,plot_inf_list[[1]], NULL, plot_inf_list[[2]],NULL, plot_inf_list[[3]]),nrow = 1, ncol = 6, rel_widths = c(0.02,1,0.02,1,0.02,1),  align = "hv",
            labels = labels_inf, label_size = 25, label_x = 0.3, label_y = 1.14,label_fontface = "plain")

p_prop <- plot_grid(plotlist = plot_prop_list,  nrow = 1, ncol = 3, rel_widths = c(1,1,1), align = "hv",
            labels = labels_prop, label_size = 25, label_x = 0.3, label_y = 1.14,label_fontface = "plain")
  

p_bottom <- plot_grid(plot_grid(plotlist = list(p_inf,NULL,p_prop),  nrow = 3, ncol = 1, rel_heights = c(1,0.02,1.05), align = "hv"),
                      legend_inf,nrow = 2, rel_heights = c(0.8,0.2))

p_bottom <- ggdraw()+
  draw_plot(p_bottom, x = 0, y = 0, width = 0.999, height = 0.95)

# Final combined plot
main_plot <- plot_grid(p_rt, p_bottom, nrow = 2, rel_heights = c(1,2)) +
  draw_label("Daily infections", x = 0.015, y = 0.45, angle = 90, size = 25) +
  draw_label(TeX(r"($R_t$)"), x = 0.015, y = 0.85, angle = 90, size = 25)

top_labels <- ggdraw() +
  draw_label("Region 1", x = 0.25, y = 0.5, size = 25) +
  draw_label("Region 2", x = 0.57, y = 0.5, size = 25) +
  draw_label("Region 3", x = 0.9, y = 0.5, size = 25)

final_plot <- plot_grid(top_labels, main_plot, nrow = 2, rel_heights = c(0.05, 1))
final_plot <- plot_grid(final_plot,NULL, nrow = 2,rel_heights = c(1,0.02))
# plot(final_plot)

final_plot <- final_plot + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("figures/updated_si/simulated_data_fitting.png", plot=final_plot, width = 15, height = 15, dpi = 300)



