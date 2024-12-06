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
library(readxl)

script_directory <- this.path::this.dir()
setwd(script_directory)

#----------------------------------- data ------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# load("results/double_mobility/region_connected_rt_double_mob_xyz.Rdata")
load("results/region_disconnected_rt.Rdata")

load("results/region_connected_rt.Rdata")

inf_start_date <- as.Date("07-02-2020",format = "%d-%m-%Y")
fitting_start <-as.Date("11-03-2020", format = "%d-%m-%Y")
end_date <- as.Date("31-12-2020", format = "%d-%m-%Y")
# final_time <- stan_data_connected$final_time
first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")

second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 
Rt <- fit_connected$draws("Rt",format = "matrix")
load("data/final_pop_2020_ltla.Rdata")

M_regions <- stan_data_connected$M_regions

death_regions <- data.frame(stan_data_connected$death)
death_data_length <- stan_data_connected$death_data_length
death_regions$time <- seq(from = fitting_start-1 ,to =  end_date, by = "day")
final_time <- length(death_regions$time)
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
                                Rt_disc_max_2 = colQuantiles(est_Rt_disc[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75))
                                time = seq(from = fitting_start-1 ,to = end_date, by = "day"))
  data_est_Rt_disc <- data_est_Rt_disc %>% filter(time >= fitting_start)
  
  
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
  #---------------------------
  
  # est_inf_disc <- fit$draws("weekly_deaths",format="matrix")
  # data_est_inf_disc <- data.frame(est_inf_disc_mean = colMeans(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)]),
  #                                inf_disc_min_1 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.025),
  #                                inf_disc_max_1 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.975),
  #                                inf_disc_min_2 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.25),
  #                                inf_disc_max_2 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.75),
  #                                time = seq(from=inf_start_date ,to =  end_date, by = "week"))
  # 
  # # est_inf_disc <- fit$draws("weekly_deaths",format="matrix")       #fit_disconnected$draws("infection", format = "matrix")
  # # data_est_inf_disc <- data.frame(est_inf_disc_mean = colMeans(est_inf_disc),
  # #                                 inf_disc_min_1 = colQuantiles(est_inf_disc,prob=0.025),
  # #                                 inf_disc_max_1 = colQuantiles(est_inf_disc,prob=0.975),
  # #                                 inf_disc_min_2 = colQuantiles(est_inf_disc,prob=0.25),
  # #                                 inf_disc_max_2 = colQuantiles(est_inf_disc,prob=0.75),
  # #                                 time = seq(from=inf_start_date ,to =  end_date, by = "week"))
  # 
  # # data_est_inf_disc <- data_est_inf_disc %>%
  # #   filter(time >= fitting_start)
  # 
  # data_inf_disc_95 <- data.frame(time = data_est_inf_disc$time, inf_disc_min = data_est_inf_disc$inf_disc_min_1,
  #                                inf_disc_max = data_est_inf_disc$inf_disc_max_1, key = rep("nintyfive", length(data_est_inf_disc$time)))
  # 
  # data_inf_disc_50 <- data.frame(time = data_est_inf_disc$time, inf_disc_min = data_est_inf_disc$inf_disc_min_2,
  #                                inf_disc_max = data_est_inf_disc$inf_disc_max_2, key = rep("fifty", length(data_est_inf_disc$time)))
  # 
  # # data_inf_disc <- rbind(data_inf_disc_95, data_inf_disc_50)
  # # levels(data_inf_disc$key) <- c("ninetyfive", "fifty")
  # data_inf_disc <- data_inf_disc_95
  # data_inf_disc$key1 <- "95% CI of estimated\ndisconnectded incidence"
  # 
  #------- connected model --------------------------------------------------------------------------------------------
  
  fit <- fit_connected
  est_Rt_con <- fit$draws("Rt",format="matrix")
  data_est_Rt_con <- data.frame(est_Rt_con_mean = colMeans(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)]),
                                Rt_con_min_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                Rt_con_max_1 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                Rt_con_min_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                Rt_con_max_2 = colQuantiles(est_Rt_con[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                time = seq(from=fitting_start-1 ,to =  end_date, by = "day"))
  data_est_Rt_con <- data_est_Rt_con %>% filter(time >= fitting_start)
  
  data_Rt_con_95 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_1,
                               Rt_con_max = data_est_Rt_con$Rt_con_max_1, key = rep("nintyfive", length(data_est_Rt_con$time)))
  
  data_Rt_con_50 <- data.frame(time = data_est_Rt_con$time, Rt_con_min = data_est_Rt_con$Rt_con_min_2,
                               Rt_con_max = data_est_Rt_con$Rt_con_max_2, key = rep("fifty", length(data_est_Rt_con$time)))
  
  #  data_Rt_con <- rbind(data_Rt_con_95, data_Rt_con_50)
  #  levels(data_Rt_con$key) <- c("ninetyfive", "fifty")
  data_Rt_con <- data_Rt_con_95
  data_Rt_con$key1 <- "95% CI of \nconnected Rt"
  #-------------------------------
  
  # est_inf_con <- fit$draws("weekly_deaths",format="matrix")
  # data_est_inf_con <- data.frame(est_inf_con_mean = colMeans(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)]),
  #                                inf_con_min_1 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.025),
  #                                inf_con_max_1 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.975),
  #                                inf_con_min_2 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.25),
  #                                inf_con_max_2 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.75),
  #                                time = seq(from=inf_start_date ,to =  end_date, by = "week"))
  # 
  # # data_est_inf_con <- data_est_inf_con %>%
  # #   filter(time >= fitting_start)
  # 
  # data_inf_con_95 <- data.frame(time = data_est_inf_con$time, inf_con_min = data_est_inf_con$inf_con_min_1,
  #                               inf_con_max = data_est_inf_con$inf_con_max_1, key = rep("nintyfive", length(data_est_inf_con$time)))
  # 
  # data_inf_con_50 <- data.frame(time = data_est_inf_con$time, inf_con_min = data_est_inf_con$inf_con_min_2,
  #                               inf_con_max = data_est_inf_con$inf_con_max_2, key = rep("fifty", length(data_est_inf_con$time)))
  # 
  # #  data_inf_con <- rbind(data_inf_95, data_inf_50)
  # #  levels(data_inf_con$key) <- c("ninetyfive", "fifty")
  # data_inf_con <- data_inf_con_95
  # data_inf_con$key1 <- "95% CI of estimated\nconnectded incidence"
  # 
  Rt_threshold <- data.frame(time = data_est_Rt_disc$time, Rt = rep(1,length(data_est_Rt_disc$time)))  # for Rt threshold horizontal line
  
  #---------- plot ------------------------------------------------------------------------------------------------
  
  colors_rt <- c("Connected Rt" = "#1a9850", "Disconnected Rt" = "#b2182b")#, "Simulated Rt"="black")
  colors_incidence <- c("Estimated disconnected\nincidence" = "red4", "Estimated connected\nincidence" = "green4", "Simulated\nincidence"="coral3")
  
  plot_rt <- ggplot(data_est_Rt_disc)+
    
    geom_ribbon(data = data_Rt_disc, aes(x = time, ymin = Rt_disc_min, ymax = Rt_disc_max, fill=key1))+
    geom_ribbon(data = data_Rt_con, aes(x = time, ymin = Rt_con_min, ymax = Rt_con_max, fill=key1))+
    geom_line(data = data_est_Rt_disc, aes(x = time,y = est_Rt_disc_mean, color = "Disconnected Rt"), linewidth = 1)+
    geom_line(data = data_est_Rt_con, aes(x = time, y = est_Rt_con_mean, color = "Connected Rt"), linewidth = 1)+
    geom_line(data = Rt_threshold, aes(x=time, y = Rt),color = "black")+
    geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 1)+
    xlab("")+
    ylab(expression(R[t]))+
    scale_fill_manual(name = "",
                      values = c("95% CI of \nconnected Rt" = alpha("#1a9850", 0.25),
                                 "95% CI of \ndisconnected Rt" = alpha("#b2182b", 0.25))) +
    scale_color_manual(values = colors_rt)+
    scale_x_date(date_labels = "%b %y", date_breaks = "1 month") + 
    ggtitle(regions[i])+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 50,hjust = 0.4, vjust = 0.4,size = 17,color="black"),
          axis.text.y = element_text(size = 20,margin = margin(r=10),color="black"),
          axis.title.y = element_text(size = 20, margin=margin(r=10)),
          axis.title.x = element_text(size = 20, margin=margin(r=10)),
          plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank(),      # Increase legend title size
          legend.text = element_text(size = 20),       # Increase legend text size
          legend.key.size = unit(1.2, "cm"))+
    guides(fill=guide_legend(ncol=1))
  
  # if (i == 1){
    # plot_rt <- plot_rt + ylab(expression(R[t]))+
  #     theme(axis.title.y = element_text(size = 20, margin = margin(r=15)))
  # }else {
  #   plot_rt <- plot_rt + theme(axis.title.y = element_blank())
  # }

  # plot(plot_rt)
  assign(paste0("rt",i),plot_rt)

# death_regions$column_to_plot <- death_regions[[i]]
# plot_inf <-ggplot(data_est_inf_disc)+
# 
#   geom_ribbon(data = data_inf_disc, aes(x = time, ymin = inf_disc_min, ymax = inf_disc_max, fill=key1))+
#   geom_ribbon(data = data_inf_con, aes(x = time, ymin = inf_con_min, ymax = inf_con_max, fill=key1))+
#   geom_line(data = data_est_inf_disc, aes(x = time,y = est_inf_disc_mean, color = "Estimated disconnected\nincidence"), linewidth = 1)+
#   geom_line(data = data_est_inf_con, aes(x = time, y = est_inf_con_mean, color = "Estimated connected\nincidence"), linewidth = 1)+
#   geom_point(data = death_regions, aes(x = time, y = column_to_plot, color = "Simulated\nincidence")) +
# 
#   xlab("")+
#   ylab("Weekly deaths")+
#   scale_fill_manual(name = "",
#                     values = c("95% CI of estimated\ndisconnectded incidence" = alpha("red4", 0.25),
#                                "95% CI of estimated\nconnectded incidence" = alpha("seagreen3", 0.25))) +
#   scale_color_manual(values = colors_incidence)+
#   scale_shape_manual(values = 16)+
#   ggtitle(regions[i])+
#   scale_x_date(date_labels = "%b %y", date_breaks = "1 month") + 
#   theme_bw()+
#   theme(axis.text.x = element_text(angle =50, hjust = 0.5,vjust = 0.4,size =17,color="black"),
#         axis.text.y = element_text(size = 15),
#         axis.title.x = element_text(size = 20, margin = margin(t=10)),
#         plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
#         legend.position = "right",
#         legend.title = element_blank(),      # Increase legend title size
#         legend.text = element_text(size = 17),       # Increase legend text size
#         legend.key.size = unit(1.2, "cm"))  +
# 
#   guides(fill=guide_legend(ncol=1))
# # if (i == 1){
#   # plot_inf <- plot_inf + ylab("Incidence")+
# #     theme(axis.title.y = element_text(size = 20, margin = margin(r=10)))
# # }else {
# #   plot_inf <- plot_inf + theme(axis.title.y = element_blank())
# # }
# 
# # plot(plot_inf)
# assign(paste0("inf",i),plot_inf)
}

legend_rt <- get_legend(plot_rt)
# legend_inf <- get_legend(plot_inf) 

# Remove legends from individual plots
plot_rt_list <-  list(rt1,rt2,rt3,rt4,rt5,rt6,rt7,rt8,rt9)
# plot_inf_list <-  list(inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf8,inf9)


for (i in 1:length(plot_rt_list)){
  plot_rt_list[[i]] <- plot_rt_list[[i]] + theme(legend.position = "none")
}

# for (i in 1:length(plot_inf_list)){
#   plot_inf_list[[i]] <- plot_inf_list[[i]] + theme(legend.position = "none")
# }

p_rt <- plot_grid(do.call(plot_grid, c(plot_rt_list,  nrow = 3, ncol = 3)),legend_rt,nrow=1,rel_widths =c(2.5,0.38))
# p_inf <- do.call(plot_grid, c(plot_inf_list, nrow = 3, ncol = 3))

p_rt <- p_rt + theme(plot.background = element_rect(fill = "white", color = NA))
p_rt <- ggdraw() +
  draw_label("Estimated Rt with double mobility", fontface = 'bold', x = 0.5, y = 0.98, hjust = 0.5, size = 20) +  # Title
  draw_plot(p_rt, y = 0, height = 0.95)  # Add the combined plot
print(p_rt)

# p_inf <- p_inf + theme(plot.background = element_rect(fill = "white", color = NA))
# print(p_inf)
#---------------

# ggsave(filename = paste0("figures/estimated_rt_xyz_double_mob.png"), plot = p_rt, width=20, height=15, units="in")
# ggsave(filename = paste0("figures/estimated_data.png"), plot = p_inf, width=12, height=8, units="in")

#----- plot original weekly death data ----------------------------------------------------------------

death_data <- data.frame(stan_data_connected$death )      # form stan_data_arrangements.R
inf_start_date <- as.Date("27-01-2020",format = "%d-%m-%Y")
fitting_start <- as.Date("09-03-2020", format = "%d-%m-%Y")
end_date <- as.Date("31-12-2020", format = "%d-%m-%Y")
first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")

second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 
regions <- c("North East","North West","Yorkshire and \nthe Humber","East Midlands","West Midlands","East","London","South East","South West")
colnames(death_data) <- regions
death_data$time <- seq(from=inf_start_date ,to = end_date, by = "week")
data_long <- pivot_longer(death_data, cols = -time, names_to = "Column", values_to = "Value")

# Plot each column as a curve using ggplot
p <- ggplot(data_long, aes(x = time, y = Value, color = Column)) +
  geom_line(,linewidth = 1.2) + 
  geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 1)+
  labs(title = "Weekly death data for each region",
       x = "",
       y = "Weekly death data",
       color = "") +
  scale_x_date(date_labels = "%b %y", date_breaks = "1 month") + 
  scale_color_manual(values = c("North East" = "#a6cee3", "North West" = "#1f78b4", "Yorkshire and \nthe Humber" = "#b2df8a", 
                                "East Midlands" = "#33a02c", "West Midlands" = "#fb9a99","East"="#e31a1c","London"="#fdbf6f","South East"="#ff7f00","South West"="#cab2d6")) +  # Set specific colors for each region
  theme_bw()+
  theme(axis.text.x = element_text(angle =50, hjust = 0.5,vjust = 0.4,size =14,color="black"),
        axis.text.y = element_text(size = 16,color="black"),
        axis.title.y = element_text(size = 16, margin = margin(t=10)),
        plot.title = element_text(size=18, margin = margin(l = 15,b=10),hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 16),      # Increase legend title size
        legend.text = element_text(size = 14),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.y = unit(0.05, "cm"))
print(p)

ggsave(filename = paste0("figures/death_data.png"), plot = p, width=8, height=5, units="in")


# ---- facet plot for all the death data -------------------------------------------------------------
# daily_death_long <- 
#   daily_death %>%
#   filter(Date >= ymd("20200303")) %>%
#   pivot_longer(!c(Date,England), names_to = "region", values_to = "deaths") %>%
#   mutate(deaths = as.integer(deaths))
# 
# ggplot(data = daily_death_long, aes(Date,deaths)) +
#   geom_line() +
#   facet_grid(rows = vars(region))
# 

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


#--- plot death data---------------------------------------------------------------------------------

#--- 2nd october is monday --------------------------------------------------------------------------
death_data <- read_excel("data/daily_death_data_region.xlsx",sheet = 12)   

colnames(death_data) <- death_data[3,]
daily_death <- death_data[4:314,c(1,4,8:16)]
region_order <- colnames(daily_death)[3:11]
daily_death <- daily_death %>% mutate_all(~ifelse(is.na(.),0,.))
daily_death$Date <- as.Date(daily_death$Date, format="%d/%m/%Y")
daily_death$Date[278:311] <- seq(as.Date("28-11-2020",format = "%d-%m-%Y"),as.Date("31-12-2020",format = "%d-%m-%Y"),by="day")

daily_death <- daily_death %>% filter(Date >= as.Date("02-03-2020",format = "%d-%m-%Y"))
#-------- death data arrangements -------------------------------------------------------

daily_death <- daily_death %>% mutate(week_group = rep(1:(nrow(daily_death) %/% 7 +1), each =7, length.out = nrow(daily_death)))
daily_death[,3:11] <- lapply(daily_death[,3:11], function(x) as.numeric(x))
weekly_death <- daily_death %>% group_by(week_group) %>% summarise(across(3:11, sum, na.rm = TRUE))

weekly_death$Date <- seq(as.Date("02-03-2020",format = "%d-%m-%Y"), as.Date("31-12-2020",format = "%d-%m-%Y"), by="week")
weekly_death$source <- "Current"

weekly_death_long <-
  weekly_death %>%
  pivot_longer(!c(week_group, Date,source), names_to = "region", values_to = "deaths") %>%
  mutate(deaths = as.integer(deaths))

#----- previous weekly data --------------------------------------------------------------

load("data/final_pop_2020_ltla.Rdata") 
load("data/england_death_2020.Rdata")       
load("data/uk_regions_mobility_matrix.Rdata")
#------- regions index -----------------------------------------------------------------

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

death_regions <- data.frame(north_east = apply(death_data[,north_east_index],1,sum),
                            north_west = apply(death_data[,north_west_index],1,sum),
                            yorkshire = apply(death_data[,yorkshire_index],1,sum),
                            east_midlands = apply(death_data[,east_midlands_index],1,sum),
                            west_midlands = apply(death_data[,west_midlands_index],1,sum),
                            east = apply(death_data[,east_index],1,sum),
                            london = apply(death_data[,london_index],1,sum),
                            south_east = apply(death_data[,south_east_index],1,sum),
                            south_west = apply(death_data[,south_west_index],1,sum))
death_regions$Date <- seq(as.Date("30-12-2019",format = "%d-%m-%Y"), as.Date("31-12-2020",format = "%d-%m-%Y"), by="week")
death_regions <- death_regions %>% filter(Date >= as.Date("02-03-2020",format = "%d-%m-%Y"))
death_regions$source <- "Previous"

prev_weekly_death_long <- death_regions %>% 
  pivot_longer(!c(Date, source), names_to = "region", values_to = "deaths") %>%
  mutate(deaths = as.integer(deaths))

combined_death_long <- bind_rows(weekly_death_long,prev_weekly_death_long)



ggplot(data = combined_death_long, aes(Date, deaths, color = source)) +
  geom_line() +
  facet_grid(rows = vars(region)) +
  labs(title = "Weekly Deaths by Region",
       x = "Date",
       y = "Number of Deaths",
       color = "Data Source") +
  theme_minimal()


#---- region wise mobility matrix ---------------

load("data/absolute_matrix_2011.Rdata")
load("data/final_pop_2020_ltla.Rdata") 
# sum(colnames(mobility_matrix_2011)[2:307] == mobility_matrix_2011$place_of_work)

mobility_matrix_2011 <- mobility_matrix_2011 %>% select(- place_of_work)

pop_2020 <- pop_2020 %>% mutate(area_name = factor(area_name, levels = colnames(mobility_matrix_2011))) %>% arrange(area_name)

colnames(mobility_matrix_2011)<- sapply(colnames(mobility_matrix_2011), function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

regions <- c("north_east","north_west","yorkshire","east_midlands","west_midlands","east","london","south_east","south_west")

north_east_index <- which(pop_2020$region == "North east")
north_west_index <- which(pop_2020$region == "North west")
yorkshire_index <- which(pop_2020$region == "Yorkshire and the humber")
east_midlands_index <- which(pop_2020$region == "East midlands")
west_midlands_index <- which(pop_2020$region == "West midlands")
east_index <- which(pop_2020$region == "East")
london_index <- which(pop_2020$region == "London")
south_east_index <- which(pop_2020$region == "South east")
south_west_index <- which(pop_2020$region == "South west")
  
north_east_matrix <- apply(as.matrix(mobility_matrix_2011[north_east_index, north_east_index]), 2, function(col) col/sum(col)) 
row.names(north_east_matrix) <- pop_2020$area_name[north_east_index]

north_west_matrix <- apply(as.matrix(mobility_matrix_2011[north_west_index, north_west_index]), 2, function(col) col/sum(col))
row.names(north_west_matrix) <- pop_2020$area_name[north_west_index]

yorkshire_matrix <- apply(as.matrix(mobility_matrix_2011[yorkshire_index,yorkshire_index]), 2, function(col) col/sum(col))
row.names(yorkshire_matrix) <- pop_2020$area_name[yorkshire_index]

east_midlands_matrix <- apply(as.matrix(mobility_matrix_2011[east_midlands_index,east_midlands_index]), 2, function(col) col/sum(col))
row.names(east_midlands_matrix) <- pop_2020$area_name[east_midlands_index]

west_midlands_matrix <- apply(as.matrix(mobility_matrix_2011[west_midlands_index,west_midlands_index]), 2,function(col) col/sum(col))
row.names(west_midlands_matrix) <- pop_2020$area_name[west_midlands_index]

east_matrix <- apply(as.matrix(mobility_matrix_2011[east_index, east_index]), 2,function(col) col/sum(col))
row.names(east_matrix) <- pop_2020$area_name[east_index]

london_matrix <- apply(as.matrix(mobility_matrix_2011[london_index,london_index]), 2, function(col) col/sum(col))
row.names(london_matrix) <- pop_2020$area_name[london_index]

south_east_matrix <- apply(as.matrix(mobility_matrix_2011[south_east_index,south_east_index]), 2, function(col) col/sum(col))
row.names(south_east_matrix) <- pop_2020$area_name[south_east_index]

south_west_matrix <- apply(as.matrix(mobility_matrix_2011[south_west_index,south_west_index]), 2, function(col) col/sum(col))
row.names(south_west_matrix) <- pop_2020$area_name[south_west_index]

mobility_matrices<- list(north_east_matrix,north_west_matrix,yorkshire_matrix,east_midlands_matrix,west_midlands_matrix,east_matrix,london_matrix,south_east_matrix,south_west_matrix)

for (i in 1:length(regions)){
  
  mobility_matrix <- get(paste0(regions[i],"_matrix"))

  mobility_df <- reshape2::melt(mobility_matrix)
colnames(mobility_df) <- c("X", "Y", "Value")
mobility_df$X <- factor(mobility_df$X, labels = colnames(mobility_matrix))
mobility_df$Y <- factor(mobility_df$Y, labels = colnames(mobility_matrix))
p <- ggplot(mobility_df, aes(x = Y, y = X, fill = Value)) +
  geom_tile() +  # Create the heatmap
  geom_text(data = mobility_df %>% filter(X == Y), aes(label = sprintf("%.2f", Value)), color = "black", size = 4,fontface="bold") +  # Add values inside the tiles
  scale_fill_gradient(low = "white", high = "#5ab4ac") +  # Color gradient for the values
  theme_minimal() +
  labs(title = regions[i], x = "", y = "", fill = "Value") +  # Axis and legend labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10,color="black"),  # Rotate x-axis labels and set size
        axis.text.y = element_text(size = 10,color="black"),
        plot.title = element_text(size=18, margin = margin(l = 15,b=10),face = "bold"),
        legend.position = "none")  # Set y-axis text size

print(p)
ggsave(filename = paste0("figures/",regions[i],"mobility.png"), plot = p, width=10, height=10, units="in")
}
