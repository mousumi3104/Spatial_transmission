# the estimated infection at own region and at different regions due to the mobility
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

# rm(list = ls())
script_directory <- this.path::this.dir()
setwd(script_directory)

# this is only for Rt and infection plot 
load("data/simulated_data.Rdata")
load("results/connected_region_fitting.Rdata")
pop <- c(20000000,10000000,15000000)

fit <- fit_connected
init_R <- fit$draws("Rt",format="matrix")
Rt_connected <- fit$draws("Rt",format ="matrix")  # need to arrange
inf <- fit$draws("infection", format = "matrix")  # need to arrange

final_time <- stan_data_connected$final_time
M_regions <- stan_data_connected$M_regions 

# ------------- stan data arrange-------------------------------------------------------------
no_sample <- 20
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
                             inf_max1 = rowQuantiles(infection_out_mob,prob = 0.95))

# death_data <- data.frame(death_mean = rowMeans(weekly_deaths),
                         # death_min1 = rowQuantiles(weekly_deaths, prob = 0.05),
                         # death_max1 = rowQuantiles(weekly_deaths, prob = 0.95))


for (m in 1:M_regions){

  plot_est_infection <- data.frame(time = stan_data_connected$fitting_start:final_time, 
                               est_infection_mean = colMeans(inf[,(((m-1)*final_time)+stan_data_connected$fitting_start):(m*final_time)]),
                               est_infection_min1 = colQuantiles(inf[,(((m-1)*final_time)+stan_data_connected$fitting_start):(m*final_time)], prob=0.01),
                               est_infection_max1 = colQuantiles(inf[,(((m-1)*final_time)+stan_data_connected$fitting_start):(m*final_time)], prob=0.99))
  
  plot_data_infection <- data.frame(time = stan_data_connected$fitting_start:final_time, 
                              data_infection = stan_data_connected$data_inf[stan_data_connected$fitting_start:final_time,m])
  
  data_inf_95 <- data.frame(time = plot_est_infection$time, inf_min = plot_est_infection$est_infection_min1,
                               inf_max = plot_est_infection$est_infection_max1, key = rep("95% CI of total infection", length(plot_est_infection$time)))
  
  plot_inf_in_own <- inf_in_own_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_own$time <- 1:final_time
  
  data_inf_in_own_95 <- data.frame(time = plot_inf_in_own$time, inf_min = plot_inf_in_own$inf_min1,
                            inf_max = plot_inf_in_own$inf_max1, key = rep("95% CI of infection in own", length(plot_inf_in_own$time)))
  
  plot_inf_in_mob <- inf_in_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_mob$time <- 1:final_time
  
  data_inf_in_mob_95 <- data.frame(time = plot_inf_in_mob$time, inf_min = plot_inf_in_mob$inf_min1,
                                   inf_max = plot_inf_in_mob$inf_max1, key = rep("95% CI of infection in mob", length(plot_inf_in_mob$time)))
  
  plot_inf_out_mob <- inf_out_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_out_mob$time <- 1:final_time
  
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

p <- ggplot(data_stack_plot1, aes(x = time, y = infections, fill = types)) +
   # geom_ribbon(aes(x = time, ymin = inf_min, ymax = inf_max, fill = key), alpha =0.25,show.legend = FALSE)+
   geom_area(position = "stack") +
   geom_point(data = plot_data_infection, aes(x=time,y=data_infection, color = "observed_data"),inherit.aes = FALSE, size =1.5)+
   
   geom_ribbon(data = plot_est_infection, aes(x = time, ymin = est_infection_min1, ymax = est_infection_max1),fill = "red4",alpha = 0.55, show.legend = FALSE,inherit.aes = FALSE)+
  geom_line(data = plot_est_infection, aes(x=time, y=est_infection_mean, color = "fitted_infection"),inherit.aes = FALSE,linewidth=1) +
  scale_fill_manual(values = c("in_own" = "#c7e9b4", "in_mob" = "#fdae61", "out_mob" = "#2b83ba"),
                    labels = c("Infections driven by \nown infections", "Mobility induced \ninfections within region", "Mobility induced infections \noutside the region") ) +
  scale_color_manual(values = c("observed_data" = "#636363", "fitted_infection" = "#d7191c"),
                     labels = c("Observed \nData", "Estimated \nInfection"))+  # Labels for legends)+

# scale_fill_manual(values = c("in_own" = "#fdc086", "in_mob" = "#ffff99", "out_mob" = "#8da0cb"),
  #                   labels = c("Infections driven by \nown infections", "Mobility induced \ninfections within region", "Mobility induced infections \noutside the region") ) +
  #  # coord_fixed(ratio = 1) +
  # scale_color_manual(values = c("observed_data" = "maroon", "fitted_infection" = "red"),  # Custom colors for legends
    # labels = c("Observed \nData", "Estimated \nInfection"))+  # Labels for legends
    labs(
     title = sprintf("Region %d",m),
     x = "Time",
     y = ""
   ) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 40,hjust = 0.4, vjust = 0.4,size = 22,color="black"),
        axis.text.y = element_text(size = 22,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 22, margin=margin(r=10)),
        axis.title.x = element_text(size = 22, margin=margin(t=10)),
        plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
        legend.position = "bottom",
        plot.margin = margin(1,12,1,1),
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 20),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.y = unit(10, "cm"))+
        
        guides(fill=guide_legend(nrow=1))
assign(paste0("region",m),p)
 # print(p)
}

legend_inf <- get_legend(p) 

plot_inf_list <-  list(region1, region2, region3)
for (i in 1:length(plot_inf_list)){
  plot_inf_list[[i]] <- plot_inf_list[[i]] + theme(legend.position = "none")
}
p_inf <- plot_grid(plot_grid(plotlist=plot_inf_list, nrow = 1, ncol = 3, rel_widths = c(1,1,1), align = "hv",axis = "tblr"),
                   legend_inf,nrow=2,rel_heights = c(2.5,0.38))

p_inf <- ggdraw() +
  draw_plot(p_inf, x = 0.01, y = 0, width = 0.95, height = 1) +  
  draw_label("Number of infections",  x = 0.02, y = 0.6,  angle=90, size = 22)   # Title
  
plot(p_inf)

ggsave("figures/infection_in_out.png", plot=p_inf, width = 16, height = 5.5, dpi = 300)

