
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

# this is only for Rt and infection plot 

load("results/ltla_connected_ne_jan.Rdata")
load("data/final_pop_2020_ltla.Rdata")
pop_north_east <- pop_2020 %>% filter(region == "North east")

inf_start_date <- plot_required_date$inf_start_date
fitting_start <- plot_required_date$fitting_start_date
end_date <- plot_required_date$end_date

plot_start <- fitting_start
plot_end <- as.Date("2021-1-1", format = "%Y-%m-%d")

first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")
second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 

death_regions <- data.frame(stan_data_connected$death)
death_data_length <- stan_data_connected$death_data_length
death_regions$time <- seq(from=inf_start_date, to =  end_date, by = "week")

regions <- c("Hartlepool","Middlesbrough","Redcar and Cleveland","Stockton-on-Tees","Darlington","County Durham","Northumberland",
             "Gateshead","Newcastle upon Tyne","North Tyneside","South Tyneside","Sunderland")

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
no_sample <- 20
infection_separate <- array(data = NA, dim = c(final_time*M_regions, no_sample))
# weekly_deaths <- array(data = NA, dim = c(ceiling(final_time/7)*M_regions, no_sample))

m_separate <- cmdstan_model("ne_ltla_simulation_separate.stan")  
ind <- sample(1:nrow(Rt_connected),no_sample)

for (k in 1:no_sample){ 
  print(k)
  stan_data_separate <- list(M_regions = stan_data_connected$M_regions,
                    final_time = stan_data_connected$final_time,
                    W = stan_data_connected$W,
                    initial_seeding_day = stan_data_connected$initial_seeding_day,
                    initial_seeding = as.vector(initial_seeding[ind[k],]),
                    SI = stan_data_connected$SI,
                    f = stan_data_connected$f,
                    pop = stan_data_connected$pop,
                    Rt = matrix(Rt_connected[ind[k],],final_time,M_regions),
                    I = stan_data_connected$I,
                    ifr_noise = as.vector(ifr_noise[ind[k],]))
  
  simulated_data_separate <- m_separate$sample(data =stan_data_separate,
                             iter_sampling = 1,
                             chains = 1,
                             thin = 1, 
                             fixed_param = TRUE)  
  
  infection_separate[,k] <- as.matrix(simulated_data_separate$draws("infection"))
}

m_joint <- cmdstan_model("ne_ltla_simulation.stan")  
infection_joint <- array(data = NA, dim = c(final_time*M_regions, no_sample))

for (k in 1:no_sample){ 
  print(k)
  stan_data_joint <- list(M_regions = stan_data_connected$M_regions,
                    final_time = stan_data_connected$final_time,
                    W = stan_data_connected$W,
                    initial_seeding_day = stan_data_connected$initial_seeding_day,
                    initial_seeding = as.vector(initial_seeding[ind[k],]),
                    SI = stan_data_connected$SI,
                    f = stan_data_connected$f,
                    pop = stan_data_connected$pop,
                    C_base = diag(M_regions),#stan_data_connected$C_base,
                    C_lockdown = diag(M_regions),#stan_data_connected$C_lockdown,
                    Rt = matrix(Rt_connected[ind[k],],final_time,M_regions),
                    I = stan_data_connected$I,
                    ifr_noise = as.vector(ifr_noise[ind[k],]))
  
  simulated_data_joint <- m_joint$sample(data =stan_data_joint,
                             iter_sampling = 1,
                             chains = 1,
                             thin = 1, 
                             fixed_param = TRUE)  
  
  infection_joint[,k] <- as.matrix(simulated_data_joint$draws("infection"))
}


time = seq(from=inf_start_date ,to =  end_date, by = "day")

infection_data_separate <- data.frame(inf_mean = rowMeans(infection_separate),
                             inf_min1 = rowQuantiles(infection_separate,prob = 0.05),
                             inf_max1 = rowQuantiles(infection_separate,prob = 0.95))

infection_data_joint <- data.frame(inf_mean = rowMeans(infection_joint),
                                      inf_min1 = rowQuantiles(infection_joint,prob = 0.05),
                                      inf_max1 = rowQuantiles(infection_joint,prob = 0.95))


for (m in 1:M_regions){
  
  plot_infection_separate <- infection_data_separate[(((m-1)*final_time)+1):(m*final_time),]
  plot_infection_separate$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_infection_separate <- plot_infection_separate %>% filter(time >= plot_start & time < plot_end)
  
  data_inf_95_separate <- data.frame(time = plot_infection_separate$time, inf_min = plot_infection_separate$inf_min1,
                            inf_max = plot_infection_separate$inf_max1, key = rep("95% CI of total infection", length(plot_infection_separate$time)))
  
  
  plot_infection_joint <- infection_data_joint[(((m-1)*final_time)+1):(m*final_time),]
  plot_infection_joint$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_infection_joint <- plot_infection_joint %>% filter(time >= plot_start & time < plot_end)
  
  data_inf_95_joint <- data.frame(time = plot_infection_joint$time, inf_min = plot_infection_joint$inf_min1,
                                     inf_max = plot_infection_joint$inf_max1, key = rep("95% CI of total infection", length(plot_infection_joint$time)))
  
  
  data_inf <- rbind(data_inf_95_joint, data_inf_95_separate)
  data_inf$key <- factor(data_inf$key, levels = c("95% CI of total infection joint","95% CI of total infection separate"))
  
  breaks = sort(c(seq(ymd("20200101"),ymd("20210201"),by="months"),seq(ymd("2020-1-15"),ymd("20210131"),by="months")))
  labels = unique(date_format("%b")(sort(c(seq(ymd("20200101"),ymd("20210201"),by="months"),seq(ymd("2020-1-15"),ymd("20210131"),by="months")))))
  labels = c(as.vector(rbind(rep("",length(labels)),labels)),"","Jan","")
  colors_rt <- c("joint Rt" = "blue", "separate Rt" = "orange")#, "Simulated Rt"="black")
  
  
  p <- ggplot(plot_infection_joint) +
    geom_line(data = plot_infection_joint, aes(x = time, y = inf_mean, color = "joint Rt"),inherit.aes = FALSE,linewidth=2) +
    geom_line(data = plot_infection_separate, aes(x = time, y = inf_mean, color = "separate Rt"),inherit.aes = FALSE,linewidth=1) +
    geom_vline(xintercept = as.Date(c(first_lockdown_end)), linetype = "dashed", color = "black", linewidth = 1)+
    
    scale_color_manual(values = c("joint Rt" = "red", "separate Rt" = "blue"),
                       labels = c("joint Rt","separate Rt"))+  # Labels for legends)+
    labs(
      title = regions[m],    #sprintf("Region %d",m),
      x = "",
      y = ""
    ) +
    scale_x_date(labels = labels, breaks=breaks, limits = c(plot_start, plot_end)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 17,color="black"),
          axis.text.y = element_text(size = 22,margin = margin(r=10),color="black"),
          axis.title.y = element_text(size = 22, margin=margin(r=10)),
          axis.title.x = element_text(size = 22, margin=margin(t=10)),
          plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
          axis.ticks.x= element_line(colour=c(rep(c("black",NA), t=floor(length(labels)/2)))),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks=element_line(linewidth =1),
          panel.grid.major.x=element_line(colour=c(rep(c( "grey94",NA), t=floor(length(labels)/2)))),
          panel.grid.minor=element_blank(),
          legend.position = "bottom",
          plot.margin = margin(1,12,1,1),
          legend.title = element_blank(),      # Increase legend title size
          legend.text = element_text(size = 25),       # Increase legend text size
          legend.key.size = unit(1.2, "cm"),
          legend.spacing.y = unit(10, "cm"))+
    
    guides(fill=guide_legend(nrow=1))
  assign(paste0("region",m),p)
  # print(p)
}

legend_inf <- get_legend(p)
plot_inf_list <-  list(region1,region2,region3,region4,region5,region6,region7,region8,region9,region10,region11,region12)

for (i in 1:length(plot_inf_list)){
  plot_inf_list[[i]] <- plot_inf_list[[i]] + theme(legend.position = "none")
}
p_inf <- plot_grid(plot_grid(plotlist =  c(plot_inf_list), nrow = 4, ncol = 3,rel_widths = c(1,1,1), align = "hv",axis = "tblr"), 
                   legend_inf, nrow = 2, rel_heights = c(2.5,0.25))
p_inf <- p_inf + theme(plot.background = element_rect(fill = "white", color = NA))
p_inf <- ggdraw() +
  draw_plot(p_inf, x = 0.01, y = 0, width = 0.95, height = 1) +  
  draw_label("Number of infections",  x = 0.02, y = 0.6,  angle=90, size = 22)   # ylabel
print(p_inf)  
# ggsave(filename = paste0("figures/inf_in_out_ne_may_aug.png"), plot = p_inf, width=18, height=15, units="in")

