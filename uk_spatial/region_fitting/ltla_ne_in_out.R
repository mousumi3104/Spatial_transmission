# the estimated infection from three different sources (Fig 6,7)
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
library(latex2exp)
library(ggplotify)


script_directory <- this.path::this.dir()
setwd(script_directory)

# this is only for Rt and infection plot 

load("results/updated_si/gm_ltla_connected_ne_jan.Rdata")
load("data/final_pop_2020_ltla.Rdata")
pop_north_east <- pop_2020 %>% filter(region == "North east")

inf_start_date <- plot_required_date$inf_start_date
fitting_start <- plot_required_date$fitting_start_date
end_date <- plot_required_date$end_date

plot_start <- fitting_start#ymd("20200301")
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
no_sample <- 400
infection <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_in_own <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_in_mob <- array(data = NA, dim = c(final_time*M_regions, no_sample))
infection_out_mob <- array(data = NA, dim = c(final_time*M_regions, no_sample))
# weekly_deaths <- array(data = NA, dim = c(ceiling(final_time/7)*M_regions, no_sample))

m <- cmdstan_model("ne_ltla_simulation.stan")  
ind <- sample(1:nrow(Rt_connected),no_sample)

for (k in 1:no_sample){ 
  print(k)
    stan_data <- list(M_regions = stan_data_connected$M_regions,
                      final_time = stan_data_connected$final_time,
                      W = stan_data_connected$W,
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
    # weekly_deaths[,k] <- as.matrix(simulated_data$draws("weekly_deaths"))
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

# death_data <- data.frame(death_mean = rowMeans(weekly_deaths),
#                          death_min1 = rowQuantiles(weekly_deaths, prob = 0.05),
#                          death_max1 = rowQuantiles(weekly_deaths, prob = 0.95))

plot_inf_list <- list()
plot_prop_list <- list()

for (m in 1:M_regions){
  
  plot_infection <- infection_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_infection$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_infection <- plot_infection %>% filter(time >= plot_start & time <= plot_end)
  
  data_inf_95 <- data.frame(time = plot_infection$time, inf_min = plot_infection$inf_min1,
                            inf_max = plot_infection$inf_max1, key = rep("95% CI of total infection", length(plot_infection$time)))
  
  plot_inf_in_own <- inf_in_own_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_own$time <- seq(from=inf_start_date ,to = end_date, by = "day")
  
  plot_inf_in_own <- plot_inf_in_own %>% filter(time >= plot_start & time <= plot_end)
  
  data_inf_in_own_95 <- data.frame(time = plot_inf_in_own$time, inf_min = plot_inf_in_own$inf_min1,
                            inf_max = plot_inf_in_own$inf_max1, key = rep("95% CI of infection in own", length(plot_inf_in_own$time)))
  
  plot_inf_in_mob <- inf_in_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_in_mob$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_inf_in_mob <- plot_inf_in_mob %>% filter(time >= plot_start & time <= plot_end)
  
  data_inf_in_mob_95 <- data.frame(time = plot_inf_in_mob$time, inf_min = plot_inf_in_mob$inf_min1,
                                   inf_max = plot_inf_in_mob$inf_max1, key = rep("95% CI of infection in mob", length(plot_inf_in_mob$time)))
  
  plot_inf_out_mob <- inf_out_mob_data[(((m-1)*final_time)+1):(m*final_time),]
  plot_inf_out_mob$time <- seq(from=inf_start_date ,to =  end_date, by = "day")
  
  plot_inf_out_mob <- plot_inf_out_mob %>% filter(time >= plot_start & time <= plot_end)
  
  data_inf_out_mob_95 <- data.frame(time = plot_inf_out_mob$time, inf_min = plot_inf_out_mob$inf_min1,
                               inf_max = plot_inf_out_mob$inf_max1, key = rep("95% CI of infection out mob", length(plot_inf_out_mob$time)))
  
  data_inf <- rbind(data_inf_95, data_inf_in_own_95, data_inf_in_mob_95, data_inf_out_mob_95)
  data_inf$key <- factor(data_inf$key, levels = c("95% CI of total infection", "95% CI of infection in own", "95% CI of infection in mob","95% CI of infection out mob"))
  
data_stack_plot <-  data.frame(in_mob = plot_inf_in_mob$inf_mean, in_own = plot_inf_in_own$inf_mean, out_mob = plot_inf_out_mob$inf_mean, time = plot_inf_in_mob$time)

#-------------------------------------------------------------------------------------------------------------------------
#_______ plot stack plot of three types of infections over time___________________________________________________________
#-------------------------------------------------------------------------------------------------------------------------

data_stack_plot1 <- data.frame(
  time = rep(data_stack_plot$time, times = 3),
  types = factor(rep(c("in_own", "in_mob", "out_mob"), each = length(data_stack_plot$time)),levels = c("in_own", "in_mob", "out_mob")),
  infections = c(data_stack_plot$in_own, data_stack_plot$in_mob, data_stack_plot$out_mob)
)
# xlabel <- if (m %in% c(10,11,12)){ "Time"}else{""}

breaks = sort(c(seq(ymd("20200101"),ymd("20210101"),by="months"),seq(ymd("2020-1-15"),ymd("2020-12-15"),by="months")))
labels = unique(date_format("%b")(sort(c(seq(ymd("20200101"),ymd("20210101"),by="months"),seq(ymd("2020-1-15"),ymd("20210131"),by="months")))))
labels = c(as.vector(rbind(rep("",length(labels)),labels)),"")

p <- ggplot(data_stack_plot1, aes(x = time, y = infections, fill = types)) +
  geom_area(position = "stack") +
    geom_ribbon(data = data_inf_out_mob_95, aes(x = time, ymin = inf_min, ymax = inf_max), fill = "#7570b3", alpha = 0.5, inherit.aes = FALSE) +
    geom_ribbon(data = data_inf_in_mob_95, aes(x = time, ymin = plot_inf_out_mob$inf_mean  + inf_min, ymax = plot_inf_out_mob$inf_mean + inf_max), fill = "deepskyblue4", alpha = 0.5, inherit.aes = FALSE) +
    geom_ribbon(data = data_inf_in_own_95, aes(x = time, ymin = plot_inf_out_mob$inf_mean  + plot_inf_in_mob$inf_mean + inf_min, ymax = plot_inf_out_mob$inf_mean  + plot_inf_in_mob$inf_mean + inf_max), fill = "#ff7f00", alpha = 0.5, inherit.aes = FALSE) +
    geom_line(data = plot_inf_out_mob, aes(x = time, y = inf_mean),color = "#7570b3",inherit.aes = FALSE, linewidth = 1)+
    geom_line(data = plot_inf_in_mob, aes(x = time, y = plot_inf_out_mob$inf_mean  + inf_mean),color = "deepskyblue4",inherit.aes = FALSE,linewidth = 1)+
    geom_line(data = plot_inf_in_own, aes(x = time, y = plot_inf_out_mob$inf_mean  + plot_inf_in_mob$inf_mean +inf_mean),color = "#ff7f00",inherit.aes = FALSE, linewidth = 1)+
    geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end, second_lockdown_start,second_lockdown_end)), linetype = "dotted", color = "grey30", linewidth = 1)+
    scale_fill_manual(values = c("in_own" = alpha("#ff7f00",0.3), "in_mob" = alpha("deepskyblue4",0.3), "out_mob" = alpha("#7570b3",0.3)),
                      labels = c(TeX(r"(\overset{\normalsize{Infections driven by}}{\normalsize{own infected population}})"),
                                 TeX(r"(\overset{\normalsize{Mobility induced}}{\normalsize{infections within region}})") ,
                                 TeX(r"(\overset{\normalsize{Mobility induced infections}}{\normalsize{outside the region}})"))) +
    labs(
    title = regions[m],    #sprintf("Region %d",m),
    x = "",
    y = ""
  ) +
  scale_x_date(labels = labels, breaks=breaks, limits = c(plot_start, plot_end)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 0.4, vjust = 0.4,size = 17,color="black"),
        axis.text.y = element_text(size = 22,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 22, margin=margin(r=10)),
        axis.title.x = element_text(size = 22, margin=margin(t=10)),
        plot.title = element_text(size=21, margin = margin(l = 15,b=10),hjust = 0.5),
        axis.ticks.x= element_line(colour=c(rep(c(NA,"black"), t=floor(length(labels)/2)))),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks=element_line(linewidth =1),
        panel.grid.major.x=element_line(colour=c(rep(c( "grey94",NA), t=floor(length(labels)/2)))),
        panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", linewidth=1.5),
        legend.position = "bottom",
        plot.margin = margin(1,12,1,1),
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 25),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.y = unit(10, "cm"))+
  guides(fill=guide_legend(nrow=1))

 legend_inf <- get_legend(p)
 p<- p+theme(legend.position = "none")
 assign(paste0("region",m),p)
 plot_inf_list[[m]] <- ggplotGrob(p)

#-------------- plot proportion --------------------------------------------------------------------------------------------------
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
  geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end, second_lockdown_start,second_lockdown_end)), linetype = "dotted", color = "grey30", linewidth = 1)+
  labs(
    title = regions[m], 
    x = "",
    y = ""
  ) +
  scale_x_date(labels = labels, breaks=breaks, limits = c(plot_start, plot_end)) +
  theme_bw() + scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.4, vjust = 0.4,size = 17,color="black"),
        axis.text.y = element_text(size = 25,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 25, margin=margin(r=10)),
        axis.title.x = element_text(size = 25, margin=margin(r=10)),
        plot.title = element_text(size=25, margin = margin(l = 15,b=10),hjust = 0.2),
        axis.ticks.x= element_line(colour=c(rep(c(NA,"black"), t=floor(length(labels)/2)))),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major.x=element_line(colour=c(rep(c( "grey94",NA), t=floor(length(labels)/2)))),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        legend.position = "right",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 25),       # Increase legend text size
        legend.key.size = unit(2, "cm"))+
  guides(fill=guide_legend(nrow=1))
p_prop<- p_prop+theme(legend.position = "none")
# print(p_prop)
plot_prop_list[[m]] <- ggplotGrob(p_prop)
}

#---- arrange individual plots --------------------------------------------------------------------
# absolute sources
for (i in 1:M_regions) {plot_inf_list[[i]] <- as.ggplot(plot_inf_list[[i]])}

index = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)")

p_inf <- ggdraw() +
  draw_plot(plot_grid(plot_grid(plotlist = plot_inf_list,  nrow = 4, ncol = 3,rel_heights = c(1,1,1,1), rel_widths = c(1,1,1), align = "hv",axis = "tblr",
                                labels = index,
                                label_size = 20,
                                label_x = 0.2,
                                label_y = 1.01),
                      legend_inf, nrow=2, rel_heights = c(2.5,0.1)),x = 0.02, y = 0, width = 0.95, height = 1)+
  draw_label("Daily number of infections", x = 0.02, y = 0.6, angle = 90, size = 25)

p_final <- ggdraw() +
  draw_plot(p_inf, x = 0, y = 0.03, width = 1, height = 0.968)

p_final <- p_final + theme(plot.background = element_rect(fill = "white", color = NA))

print(p_final)  
ggsave(filename = paste0("figures/updated_si/inf_in_out_ne_jan.png"), plot = p_final, width=15, height=15, units="in")

#proportion-------

for (i in 1:M_regions) {plot_prop_list[[i]] <- as.ggplot(plot_prop_list[[i]])}

index = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)")

p_prop <- ggdraw() +
  draw_plot(plot_grid(plot_grid(plotlist = plot_prop_list,  nrow = 4, ncol = 3,rel_heights = c(1,1,1,1), rel_widths = c(1,1,1), align = "hv",axis = "tblr",
                                labels = index,
                                label_size = 20,
                                label_x = 0.25,
                                label_y = 1.01),
                      legend_inf, nrow=2, rel_heights = c(2.5,0.1)),x = 0.02, y = 0, width = 0.95, height = 1)+
  draw_label("Proportion of daily infections", x = 0.025, y = 0.6, angle = 90, size = 25)

p_final <- ggdraw() +
  draw_plot(p_prop, x = 0, y = 0.03, width = 1, height = 0.968)

p_final <- p_final + theme(plot.background = element_rect(fill = "white", color = NA))

print(p_final)  
ggsave(filename = paste0("figures/updated_si/prop_in_out_ne_jan.png"), plot = p_final, width=15, height=15, units="in")



#-------------------------------------------------------------------------------------------------------------------------
#_______ data for stack plot total infection before after lock down __________________________________________________________
#-------------------------------------------------------------------------------------------------------------------------

for (m in 1:M_regions){
  
  print(m)
  
  time_before_fitting_start <- length(seq(inf_start_date,fitting_start-1,by="day"))
  time_before_first_lockdown <- length(seq(fitting_start, first_lockdown_start-1,by="day"))
  time_during_first_lockdown <- length(seq(first_lockdown_start, first_lockdown_end, by="day"))
  time_after_first_lockdown <- length(seq(first_lockdown_end+1, second_lockdown_start-1, by="day"))
  time_during_second_lockdown <- length(seq(second_lockdown_start, second_lockdown_end, by="day"))
  time_after_second_lockdown <- length(seq(second_lockdown_end+1, plot_end-1, by="day"))
  
  
  total_inf_in_own <- infection_in_own[(((m-1)*final_time)+1):(m*final_time),]
  total_inf_in_mob <- infection_in_mob[(((m-1)*final_time)+1):(m*final_time),]
  total_inf_out_mob <- infection_out_mob[(((m-1)*final_time)+1):(m*final_time),]
  
  #-----------------------------------------------------------------------
  
  #before first lockdown  
  inf_in_own_1 <- colSums(total_inf_in_own[(time_before_fitting_start +1) : (time_before_fitting_start+time_before_first_lockdown),])
  inf_in_own_data1 <- data.frame(mean_inf = mean(inf_in_own_1),
                                 inf_min = quantile(inf_in_own_1, prob = 0.05),
                                 inf_max = quantile(inf_in_own_1, prob = 0.95))
  
  inf_in_mob_1 <- colSums(total_inf_in_mob[(time_before_fitting_start +1) : (time_before_fitting_start+time_before_first_lockdown),])
  inf_in_mob_data1 <- data.frame(mean_inf = mean(inf_in_mob_1),
                                 inf_min = quantile(inf_in_mob_1, prob = 0.05),
                                 inf_max = quantile(inf_in_mob_1, prob = 0.95))
  
  inf_out_mob_1 <- colSums(total_inf_out_mob[(time_before_fitting_start +1) : (time_before_fitting_start+time_before_first_lockdown),])
  inf_out_mob_data1 <- data.frame(mean_inf = mean(inf_out_mob_1),
                                  inf_min = quantile(inf_out_mob_1, prob = 0.05),
                                  inf_max = quantile(inf_out_mob_1, prob = 0.95))
  
  before_first_lkdn <- data.frame(group = c("Infections driven by \nown infected population", "Mobility induced \ninfections within region", "Mobility induced infections \noutside the region"),
                                  mean_inf = c(inf_in_own_data1$mean_inf, inf_in_mob_data1$mean_inf, inf_out_mob_data1$mean_inf),
                                  min_inf = c(inf_in_own_data1$inf_min, inf_in_mob_data1$inf_min, inf_out_mob_data1$inf_min),
                                  max_inf = c(inf_in_own_data1$inf_max, inf_in_mob_data1$inf_max, inf_out_mob_data1$inf_max),
                                  inf_perct = c(inf_in_own_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                inf_in_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                inf_out_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf)))
  
  bar_percentage1 <- ggplot(before_first_lkdn, aes(x = "", y = inf_perct, fill = group)) +
    geom_bar(stat = "identity", position = "stack", color = "black", width = 0.35) +  # Stacked bars
    geom_text(aes(label = scales::percent(inf_perct, accuracy = 1)), 
              position = position_stack(vjust = 0.5), size = 5, color = "black") +  # Labels inside the boxes
    scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0, 1)) +  # Percentage scale
    scale_fill_manual(values = c("Infections driven by \nown infected population" = "skyblue", "Mobility induced \ninfections within region" = "tomato", "Mobility induced infections \noutside the region" = "gold")) +  # Custom colors
    labs(
      title = "Before \nfirst lockdown",
      x = "",
      y = "",
      fill = "Infection Source"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=12),  # Hide x-axis labels for cleaner look
      axis.title = element_text(size=12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 15),
      legend.title = element_text(size = 20),
      legend.text = element_blank(),#text(size = 20),
      legend.position = "right"
    )
  
  bar_percentage1  <- bar_percentage1  + theme(legend.position = "none")
  
  #---------------------------------------------------------
  
  #during first lockdown  
  time_1 <- time_before_fitting_start + time_before_first_lockdown
  time_2 <- time_before_fitting_start + time_before_first_lockdown + time_during_first_lockdown
  
  inf_in_own_1 <- colSums(total_inf_in_own[(time_1 + 1) : time_2,])
  inf_in_own_data1 <- data.frame(mean_inf = mean(inf_in_own_1),
                                 inf_min = quantile(inf_in_own_1, prob = 0.05),
                                 inf_max = quantile(inf_in_own_1, prob = 0.95))
  
  inf_in_mob_1 <- colSums(total_inf_in_mob[(time_1 + 1) : (time_2),])
  inf_in_mob_data1 <- data.frame(mean_inf = mean(inf_in_mob_1),
                                 inf_min = quantile(inf_in_mob_1, prob = 0.05),
                                 inf_max = quantile(inf_in_mob_1, prob = 0.95))
  
  inf_out_mob_1 <- colSums(total_inf_out_mob[(time_1 + 1) : (time_2),])
  inf_out_mob_data1 <- data.frame(mean_inf = mean(inf_out_mob_1),
                                  inf_min = quantile(inf_out_mob_1, prob = 0.05),
                                  inf_max = quantile(inf_out_mob_1, prob = 0.95))
  
  during_first_lkdn <- data.frame(group = c("Infections driven by \nown infected population", "Mobility induced \ninfections within region", "Mobility induced infections \noutside the region"),
                                  mean_inf = c(inf_in_own_data1$mean_inf, inf_in_mob_data1$mean_inf, inf_out_mob_data1$mean_inf),
                                  min_inf = c(inf_in_own_data1$inf_min, inf_in_mob_data1$inf_min, inf_out_mob_data1$inf_min),
                                  max_inf = c(inf_in_own_data1$inf_max, inf_in_mob_data1$inf_max, inf_out_mob_data1$inf_max),
                                  inf_perct = c(inf_in_own_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                inf_in_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                inf_out_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf)))
  
  #---------------------------------------------------------------
  bar_percentage2 <- ggplot(during_first_lkdn, aes(x = "", y = inf_perct, fill = group)) +
    geom_bar(stat = "identity", position = "stack", color = "black", width = 0.35) +  # Stacked bars
    geom_text(aes(label = scales::percent(inf_perct, accuracy = 1)), 
              position = position_stack(vjust = 0.5), size = 5, color = "black") +  # Labels inside the boxes
    scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0, 1)) +  # Percentage scale
    scale_fill_manual(values = c("Infections driven by \nown infected population" = "skyblue", "Mobility induced \ninfections within region" = "tomato", "Mobility induced infections \noutside the region" = "gold")) +  # Custom colors
    labs(
      title = "During \nfirst lockdown",
      x = "",
      y = "",
      fill = "Infection Source"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=12),  # Hide x-axis labels for cleaner look
      axis.title = element_text(size=12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 15),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      legend.position = "right"
    )
  
  bar_percentage2  <- bar_percentage2  + theme(legend.position = "none")
  
  #---------------------------------------------------------
  
  # after first lockdown and before second lockdown
  
  time_1 <- time_before_fitting_start + time_before_first_lockdown + time_during_first_lockdown
  time_2 <- time_before_fitting_start + time_before_first_lockdown + time_during_first_lockdown + time_after_first_lockdown
  
  inf_in_own_1 <- colSums(total_inf_in_own[(time_1 + 1) : time_2,])
  inf_in_own_data1 <- data.frame(mean_inf = mean(inf_in_own_1),
                                 inf_min = quantile(inf_in_own_1, prob = 0.05),
                                 inf_max = quantile(inf_in_own_1, prob = 0.95))
  
  inf_in_mob_1 <- colSums(total_inf_in_mob[(time_1 + 1) : (time_2),])
  inf_in_mob_data1 <- data.frame(mean_inf = mean(inf_in_mob_1),
                                 inf_min = quantile(inf_in_mob_1, prob = 0.05),
                                 inf_max = quantile(inf_in_mob_1, prob = 0.95))
  
  inf_out_mob_1 <- colSums(total_inf_out_mob[(time_1 + 1) : (time_2),])
  inf_out_mob_data1 <- data.frame(mean_inf = mean(inf_out_mob_1),
                                  inf_min = quantile(inf_out_mob_1, prob = 0.05),
                                  inf_max = quantile(inf_out_mob_1, prob = 0.95))
  
  after_first_lkdn <- data.frame(group = c("Infections driven by \nown infected population", "Mobility induced \ninfections within region", "Mobility induced infections \noutside the region"),
                                 mean_inf = c(inf_in_own_data1$mean_inf, inf_in_mob_data1$mean_inf, inf_out_mob_data1$mean_inf),
                                 min_inf = c(inf_in_own_data1$inf_min, inf_in_mob_data1$inf_min, inf_out_mob_data1$inf_min),
                                 max_inf = c(inf_in_own_data1$inf_max, inf_in_mob_data1$inf_max, inf_out_mob_data1$inf_max),
                                 inf_perct = c(inf_in_own_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                               inf_in_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                               inf_out_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf)))
  
  #---------------------------------------------------------------
  bar_percentage3 <- ggplot(after_first_lkdn, aes(x = "", y = inf_perct, fill = group)) +
    geom_bar(stat = "identity", position = "stack", color = "black", width = 0.35) +  # Stacked bars
    geom_text(aes(label = scales::percent(inf_perct, accuracy = 1)), 
              position = position_stack(vjust = 0.5), size = 5, color = "black") +  # Labels inside the boxes
    scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0, 1)) +  # Percentage scale
    scale_fill_manual(values = c("Infections driven by \nown infected population" = "skyblue", "Mobility induced \ninfections within region" = "tomato", "Mobility induced infections \noutside the region" = "gold")) +  # Custom colors
    labs(
      title = "Between first \nand second lockdown",
      x = "",
      y = "",
      fill = "Infection Source"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=12),  # Hide x-axis labels for cleaner look
      axis.title = element_text(size=12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 15),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      legend.position = "right")
  
  bar_percentage3  <- bar_percentage3  + theme(legend.position = "none")
  
  #---------------------------------------------------------
  
  # during second lockdown
  
  time_1 <- time_before_fitting_start + time_before_first_lockdown + time_during_first_lockdown + time_after_first_lockdown
  time_2 <- time_before_fitting_start + time_before_first_lockdown + time_during_first_lockdown + time_after_first_lockdown + time_during_second_lockdown
  
  inf_in_own_1 <- colSums(total_inf_in_own[(time_1 + 1) : time_2,])
  inf_in_own_data1 <- data.frame(mean_inf = mean(inf_in_own_1),
                                 inf_min = quantile(inf_in_own_1, prob = 0.05),
                                 inf_max = quantile(inf_in_own_1, prob = 0.95))
  
  inf_in_mob_1 <- colSums(total_inf_in_mob[(time_1 + 1) : (time_2),])
  inf_in_mob_data1 <- data.frame(mean_inf = mean(inf_in_mob_1),
                                 inf_min = quantile(inf_in_mob_1, prob = 0.05),
                                 inf_max = quantile(inf_in_mob_1, prob = 0.95))
  
  inf_out_mob_1 <- colSums(total_inf_out_mob[(time_1 + 1) : (time_2),])
  inf_out_mob_data1 <- data.frame(mean_inf = mean(inf_out_mob_1),
                                  inf_min = quantile(inf_out_mob_1, prob = 0.05),
                                  inf_max = quantile(inf_out_mob_1, prob = 0.95))
  
  during_second_lkdn <- data.frame(group = c("Infections driven by \nown infected population", "Mobility induced \ninfections within region", "Mobility induced infections \noutside the region"),
                                   mean_inf = c(inf_in_own_data1$mean_inf, inf_in_mob_data1$mean_inf, inf_out_mob_data1$mean_inf),
                                   min_inf = c(inf_in_own_data1$inf_min, inf_in_mob_data1$inf_min, inf_out_mob_data1$inf_min),
                                   max_inf = c(inf_in_own_data1$inf_max, inf_in_mob_data1$inf_max, inf_out_mob_data1$inf_max),
                                   inf_perct = c(inf_in_own_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                 inf_in_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                 inf_out_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf)))
  
  #---------------------------------------------------------------

  bar_percentage4 <- ggplot(during_second_lkdn, aes(x = "", y = inf_perct, fill = group)) +
    geom_bar(stat = "identity", position = "stack", color = "black", width = 0.35) +  # Stacked bars
    geom_text(aes(label = scales::percent(inf_perct, accuracy = 1)), 
              position = position_stack(vjust = 0.5), size = 5, color = "black") +  # Labels inside the boxes
    scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0, 1)) +  # Percentage scale
    scale_fill_manual(values = c("Infections driven by \nown infected population" = "skyblue", "Mobility induced \ninfections within region" = "tomato", "Mobility induced infections \noutside the region" = "gold")) +  # Custom colors
    labs(
      title = "During \nsecond lockdown",
      x = "",
      y = "",
      fill = "Infection Source"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=12),  # Hide x-axis labels for cleaner look
      axis.title = element_text(size=12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 15),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      legend.position = "right")
  
  bar_percentage4  <- bar_percentage4  + theme(legend.position = "none")
  
  #---------------------------------------------------------
  
  # after second lockdown
  
  time_1 <- time_before_fitting_start + time_before_first_lockdown + time_during_first_lockdown + time_after_first_lockdown + time_during_second_lockdown
  time_2 <- time_before_fitting_start + time_before_first_lockdown + time_during_first_lockdown + time_after_first_lockdown + time_during_second_lockdown + time_after_second_lockdown
  
  inf_in_own_1 <- colSums(total_inf_in_own[(time_1 + 1) : time_2,])
  inf_in_own_data1 <- data.frame(mean_inf = mean(inf_in_own_1),
                                 inf_min = quantile(inf_in_own_1, prob = 0.05),
                                 inf_max = quantile(inf_in_own_1, prob = 0.95))
  
  inf_in_mob_1 <- colSums(total_inf_in_mob[(time_1 + 1) : (time_2),])
  inf_in_mob_data1 <- data.frame(mean_inf = mean(inf_in_mob_1),
                                 inf_min = quantile(inf_in_mob_1, prob = 0.05),
                                 inf_max = quantile(inf_in_mob_1, prob = 0.95))
  
  inf_out_mob_1 <- colSums(total_inf_out_mob[(time_1 + 1) : (time_2),])
  inf_out_mob_data1 <- data.frame(mean_inf = mean(inf_out_mob_1),
                                  inf_min = quantile(inf_out_mob_1, prob = 0.05),
                                  inf_max = quantile(inf_out_mob_1, prob = 0.95))
  
  after_second_lkdn <- data.frame(group = c("Infections driven by \nown infected population", "Mobility induced \ninfections within region", "Mobility induced infections \noutside the region"),
                                  mean_inf = c(inf_in_own_data1$mean_inf, inf_in_mob_data1$mean_inf, inf_out_mob_data1$mean_inf),
                                  min_inf = c(inf_in_own_data1$inf_min, inf_in_mob_data1$inf_min, inf_out_mob_data1$inf_min),
                                  max_inf = c(inf_in_own_data1$inf_max, inf_in_mob_data1$inf_max, inf_out_mob_data1$inf_max),
                                  inf_perct = c(inf_in_own_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                inf_in_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf),
                                                inf_out_mob_data1$mean_inf / sum(inf_in_own_data1$mean_inf + inf_in_mob_data1$mean_inf + inf_out_mob_data1$mean_inf)))
  
  #---------------------------------------------------------------
  bar_percentage5 <- ggplot(after_second_lkdn, aes(x = "", y = inf_perct, fill = group)) +
    geom_bar(stat = "identity", position = "stack", color = "black", width = 0.35) +  # Stacked bars
    geom_text(aes(label = scales::percent(inf_perct, accuracy = 1)), 
              position = position_stack(vjust = 0.5), size = 5, color = "black") +  # Labels inside the boxes
    scale_y_continuous(breaks = NULL, labels = NULL, limits = c(0, 1)) +  # Percentage scale
    scale_fill_manual(values = c("Infections driven by \nown infected population" = "skyblue", "Mobility induced \ninfections within region" = "tomato", "Mobility induced infections \noutside the region" = "gold")) +  # Custom colors
    labs(
      title = "After \nsecond lockdown",
      x = "",
      y = "",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=12),  # Hide x-axis labels for cleaner look
      axis.title = element_text(size=12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 15),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.position = "bottom"
    )
  
  legend_inf <- get_legend(bar_percentage5)
  bar_percentage5  <- bar_percentage5  + theme(legend.position = "none")
  
  
  plot_inf_list <-  list(bar_percentage1, bar_percentage2, bar_percentage3, bar_percentage4, bar_percentage5)
  
  p_inf <- plot_grid(plotlist =  c(plot_inf_list), nrow = 1, ncol = 5,rel_widths = c(1,1,1,1,1), align = "hv",axis = "tblr")


  p_inf <- p_inf + theme(plot.background = element_rect(fill = "white", color = NA))

  p_inf <- ggdraw() +
  draw_plot(p_inf, x = 0.01, y = 0, width = 1, height = 0.9) +  
  draw_label(regions[m],  x = 0.5, y = 0.95,  angle = 0, size = 22, fontface = "bold")   
  
  assign(paste0("inf_barplot",m), p_inf)

  # ggsave(filename = paste0("figures/infection_source_",regions[m],".png"), plot = p_inf, width=8, height=4, units="in")
}

plot_inf_list <-  list(inf_barplot1, inf_barplot2, inf_barplot3, inf_barplot4)

p_inf1 <- plot_grid(plot_grid(plotlist =  c(plot_inf_list), nrow = 4, ncol = 1,rel_heights = c(1,1,1,1), align = "hv",axis = "tblr"),
          legend_inf, nrow = 2, ncol = 1, rel_heights = c(0.95,0.05), align = "hv",axis = "tblr")
p_inf1 <- p_inf1 + theme(plot.background = element_rect(fill = "white", color = NA))
plot(p_inf1)

plot_inf_list <-  list(inf_barplot5, inf_barplot6, inf_barplot7, inf_barplot8)
p_inf2 <- plot_grid(plot_grid(plotlist =  c(plot_inf_list), nrow = 4, ncol = 1,rel_heights = c(1,1,1,1), align = "hv",axis = "tblr"),
                    legend_inf, nrow = 2, ncol = 1, rel_heights = c(0.95,0.05), align = "hv",axis = "tblr")
p_inf2 <- p_inf2 + theme(plot.background = element_rect(fill = "white", color = NA))
print(p_inf2)

plot_inf_list <-  list(inf_barplot9, inf_barplot10, inf_barplot11, inf_barplot12)
p_inf3 <- plot_grid(plot_grid(plotlist =  c(plot_inf_list), nrow = 4, ncol = 1,rel_heights = c(1,1,1,1), align = "hv",axis = "tblr"),
                    legend_inf, nrow = 2, ncol = 1, rel_heights = c(0.95,0.05), align = "hv",axis = "tblr")
p_inf3 <- p_inf3 + theme(plot.background = element_rect(fill = "white", color = NA))
print(p_inf3)

# ggsave(filename = paste0("figures/inf_barplot_ne_1.png"), plot = p_inf1, width=13, height=15, units="in")
# ggsave(filename = paste0("figures/inf_barplot_ne_2.png"), plot = p_inf2, width=13, height=15, units="in")
# ggsave(filename = paste0("figures/inf_barplot_ne_3.png"), plot = p_inf3, width=13, height=15, units="in")
