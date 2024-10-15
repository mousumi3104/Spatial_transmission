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

script_directory <- this.path::this.dir()
setwd(script_directory)

load("results/region_connected_rt_xyz.Rdata")
load("data/final_pop_2020_ltla.Rdata")

inf_start_date <- as.Date("27-01-2020",format = "%d-%m-%Y")
fitting_start <- as.Date("09-03-2020", format = "%d-%m-%Y")
end_date <- as.Date("31-12-2020", format = "%d-%m-%Y")

first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")
second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 

death_regions <- data.frame(stan_data_connected$death)
death_data_length <- stan_data_connected$death_data_length
death_regions$time <- seq(from=inf_start_date, to =  end_date, by = "week")

regions <- c("North East","North West","Yorkshire and the Humber","East Midlands","West Midlands","East","London","South East","South West")

fit <- fit_connected
Rt_connected <- fit$draws("Rt",format ="matrix")  # need to arrange
inf <- fit$draws("infection", format = "matrix")  # need to arrange
initial_seeding <- fit$draws("initial_seeding", format="matrix")
x1 <- fit$draws("x1", format ="matrix")
y1 <- fit$draws("y1", format ="matrix")
z1 <- fit$draws("z1", format ="matrix")
mu <- fit$draws("mu", format ="matrix")
weekly_var <- fit$draws("weekly_var", format = "matrix")
weekly_effect_d <- fit$draws("weekly_effect_d", format = "matrix")  # need to arrange
ifr_noise <- fit$draws("ifr_noise", format = "matrix")

final_time <- stan_data_connected$final_time
M_regions <- stan_data_connected$M_regions 

# ------------- stan data arrange-------------------------------------------------------------
no_sample <- 200
infection <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_in_own <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_in_mob <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_out_mob <- array(data = NA, dim = c(final_time*M_regions, no_sample))
weekly_deaths <- array(data = NA, dim = c(ceiling(final_time/7)*M_regions, no_sample))

m <- cmdstan_model("uk_region_simulation.stan")  
ind <- sample(1:2000,no_sample)

for (k in 1:no_sample){ 
    stan_data <- list(M_regions = stan_data_connected$M_regions,
                      final_time = stan_data_connected$final_time,
                      W = stan_data_connected$W,
                      gmobility = stan_data_connected$gmobility,
                      initial_seeding_day = stan_data_connected$initial_seeding_day,
                      initial_seeding = as.vector(initial_seeding[ind[k],]),
                      SI = stan_data_connected$SI,
                      f = stan_data_connected$f,
                      pop = stan_data_connected$pop,
                      C_base = stan_data_connected$C_base,
                      C_lockdown = stan_data_connected$C_lockdown,
                      Rt = matrix(Rt_connected[ind[k],],final_time,M_regions),
                      I = stan_data_connected$I,
                      ifr_noise = as.vector(ifr_noise[ind[k],]))
    
    simulated_data <- m$sample(data =stan_data,
                               iter_sampling = 1,
                               chains = 1,
                               thin = 1, 
                               fixed_param = TRUE)  
    
    infection[,k] <- as.matrix(simulated_data$draws("infection"))
    infection_in_own[,k] <- as.matrix(simulated_data$draws("infection_in_own"))
    infection_in_mob[,k] <- as.matrix(simulated_data$draws("infection_in_mob"))
    infection_out_mob[,k] <- as.matrix(simulated_data$draws("infection_out_mob"))
    weekly_deaths[,k] <- as.matrix(simulated_data$draws("weekly_deaths"))
}


time = seq(from=inf_start_date ,to =  end_date, by = "day")

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

death_data <- data.frame(death_mean = rowMeans(weekly_deaths),
                         death_min1 = rowQuantiles(weekly_deaths, prob = 0.05),
                         death_max1 = rowQuantiles(weekly_deaths, prob = 0.95))


for (m in 1:M_regions){
  
  plot_infection <- infection_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_infection$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_infection <- plot_infection %>% filter(time <= as.Date("20-12-2020", format = "%d-%m-%Y"))
  
  data_inf_95 <- data.frame(time = plot_infection$time, inf_min = plot_infection$inf_min1,
                               inf_max = plot_infection$inf_max1, key = rep("95% CI of total infection", length(plot_infection$time)))
  
  plot_inf_in_own <- inf_in_own_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_own$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_inf_in_own <- plot_inf_in_own %>% filter(time <= as.Date("20-12-2020", format = "%d-%m-%Y"))
  
  data_inf_in_own_95 <- data.frame(time = plot_inf_in_own$time, inf_min = plot_inf_in_own$inf_min1,
                            inf_max = plot_inf_in_own$inf_max1, key = rep("95% CI of infection in own", length(plot_inf_in_own$time)))
  
  plot_inf_in_mob <- inf_in_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_mob$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_inf_in_mob <- plot_inf_in_mob %>% filter(time <= as.Date("20-12-2020", format = "%d-%m-%Y"))
  
  data_inf_in_mob_95 <- data.frame(time = plot_inf_in_mob$time, inf_min = plot_inf_in_mob$inf_min1,
                                   inf_max = plot_inf_in_mob$inf_max1, key = rep("95% CI of infection in mob", length(plot_inf_in_mob$time)))
  
  plot_inf_out_mob <- inf_out_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_out_mob$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_inf_out_mob <- plot_inf_out_mob %>% filter(time <= as.Date("20-12-2020", format = "%d-%m-%Y"))
  
  data_inf_out_mob_95 <- data.frame(time = plot_inf_out_mob$time, inf_min = plot_inf_out_mob$inf_min1,
                               inf_max = plot_inf_out_mob$inf_max1, key = rep("95% CI of infection out mob", length(plot_inf_out_mob$time)))
  
  data_inf <- rbind(data_inf_95, data_inf_in_own_95, data_inf_in_mob_95, data_inf_out_mob_95)
  data_inf$key <- factor(data_inf$key, levels = c("95% CI of total infection", "95% CI of infection in own", "95% CI of infection in mob","95% CI of infection out mob"))
  
  total_inf_out_mob <- 
  total_inf_in_mob <- 
  total_inf_in_own <-  plot_inf_in_own 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  colors_infection <- c("Total infection" = "red4", "Infections driven by \nlocal infections" = "green4", "Mobility induced \ninfections within region"="coral3", "Mobility induced \ninfections outside region" = "blue4")
  
  plot_inf <- ggplot(data_inf)+
    geom_ribbon(aes(x = time, ymin = inf_min, ymax = inf_max, fill = key), alpha =0.25,show.legend = FALSE)+
    geom_line(data = plot_infection, aes(x = time, y = inf_mean, color = "Total infection"), linewidth = 1.3)+
    geom_line(data = plot_inf_in_own, aes(x = time, y = inf_mean, color = "Infections driven by \nlocal infections"), linewidth = 1.3)+
    geom_line(data = plot_inf_in_mob, aes(x = time, y = inf_mean, color = "Mobility induced \ninfections within region"), linewidth = 1.3)+
    geom_line(data = plot_inf_out_mob, aes(x = time, y = inf_mean, color = "Mobility induced \ninfections outside region"), linewidth = 1.3)+
    geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 1)+
    
    xlab("")+
    ylab("Daily infection")+
    # scale_fill_manual(name = "",
    #                   values = c("95% CI of total infection" = alpha("red4", 0.25),
    #                              "95% CI of infection in own" = alpha("green4", 0.25),
    #                              "95% CI of infection in mob" = alpha("coral3", 0.25),
    #                              "95% CI of infection out mob" = alpha("blue4", 0.25))) +
    scale_color_manual(values = colors_infection)+
    scale_x_date(date_labels = "%b %y", date_breaks = "1 month") + 
    ggtitle(regions[m])+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 50,hjust = 0.4, vjust = 0.4,size = 17,color="black"),
          axis.text.y = element_text(size = 20,margin = margin(r=10),color="black"),
          axis.title.y = element_text(size = 20, margin=margin(r=10)),
          axis.title.x = element_text(size = 20, margin=margin(r=10)),
          plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank(),      # Increase legend title size
          legend.text = element_text(size = 20),       # Increase legend text size
          legend.key.size = unit(1.2, "cm"),
          legend.spacing.y = unit(10, "cm"))+
    guides(fill=guide_legend(ncol=1))
  
  assign(paste0("inf",m),plot_inf)
}

legend_inf <- get_legend(plot_inf) 
plot_inf_list <-  list(inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf8,inf9)
for (i in 1:length(plot_inf_list)){
  plot_inf_list[[i]] <- plot_inf_list[[i]] + theme(legend.position = "none")
}
p_inf <- plot_grid(do.call(plot_grid, c(plot_inf_list, nrow = 3, ncol = 3)), legend_inf, nrow = 1, rel_widths = c(4,1))
p_inf <- p_inf + theme(plot.background = element_rect(fill = "white", color = NA))
print(p_inf)

# ggsave(filename = paste0("figures/inf_due_mob_london.png"), plot = inf7, width=12, height=5, units="in")





