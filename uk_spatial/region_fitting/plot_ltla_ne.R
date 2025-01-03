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

load("results/ltla_connected_ne_jan.Rdata")
load("results/ltla_disconnected_ne_jan.Rdata")

inf_start_date <- plot_required_date$inf_start_date   #as.Date("10-02-2020",format = "%d-%m-%Y")
fitting_start <- plot_required_date$fitting_start_date        #as.Date("23-03-2020", format = "%d-%m-%Y")
end_date <- plot_required_date$end_date                                           #as.Date("31-12-2020", format = "%d-%m-%Y")
plot_end <- as.Date("2020-12-31", format = "%Y-%m-%d")
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
final_time <- stan_data_connected$final_time

regions <- c("Hartlepool","Middlesbrough","Redcar and Cleveland","Stockton-on-Tees","Darlington","County Durham","Northumberland",
             "Gateshead","Newcastle upon Tyne","North Tyneside","South Tyneside","Sunderland")
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
                                time = seq(from = inf_start_date ,to = end_date, by = "day"))
  data_est_Rt_disc <- data_est_Rt_disc %>% filter(time >= fitting_start & time <= plot_end )
  
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
  
  #----------------------------------------------------------------------------------------------------------------
  Rt_threshold <- data.frame(time = data_est_Rt_disc$time, Rt = rep(1,length(data_est_Rt_disc$time)))  # for Rt threshold horizontal line
  
  #---------- plot ------------------------------------------------------------------------------------------------
  
  breaks = sort(c(seq(ymd("2020-4-1"),ymd("2020-12-31"),by="months"),seq(ymd("2020-3-15"),ymd("2020-12-31"),by="months"),ymd("2021-1-1")))
    
  labels = unique(date_format("%b")(data_est_Rt_con$time))
  labels = as.vector(rbind(labels,rep("",length(labels))))
  
  colors_rt <- c("Connected Rt" = "#1a9850", "Disconnected Rt" = "#b2182b")#, "Simulated Rt"="black")
  colors_incidence <- c("Estimated disconnected\nincidence" = "red4", "Estimated connected\nincidence" = "green4", "Simulated\nincidence"="coral3")
  
  plot_rt <- ggplot(data_est_Rt_disc)+
    geom_ribbon(data = data_Rt_disc, aes(x = time, ymin = Rt_disc_min, ymax = Rt_disc_max, fill=key1))+
    geom_ribbon(data = data_Rt_con, aes(x = time, ymin = Rt_con_min, ymax = Rt_con_max, fill=key1))+
    geom_line(data = data_est_Rt_disc, aes(x = time,y = est_Rt_disc_mean, color = "Disconnected Rt"), linewidth = 1.2)+
    geom_line(data = data_est_Rt_con, aes(x = time, y = est_Rt_con_mean, color = "Connected Rt"), linewidth = 1.2)+
    geom_line(data = Rt_threshold, aes(x=time, y = Rt),color = "black")+
    geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 1)+
    xlab("")+
    ylab("")+
    scale_fill_manual(name = "",
                      values = c("95% CI of \nconnected Rt" = alpha("#1a9850", 0.25),
                                 "95% CI of \ndisconnected Rt" = alpha("#b2182b", 0.25))) +
    scale_color_manual(values = colors_rt)+
    scale_x_date(labels = labels, breaks=breaks, limits = c(fitting_start, ymd("2021-1-1"))) + 
    scale_y_continuous(limits = c(0,6), breaks = seq(0,6,by=1))+
    ggtitle(regions[i])+
    theme_bw()+
    theme(axis.text.x = element_text(size = 15,margin = margin(r=10),color="black", hjust =0.4),
          axis.text.y = element_text(size = 17,margin = margin(r=10),color="black"),
          axis.title.y = element_text(size = 20, margin=margin(r=10)),
          axis.title.x = element_text(size = 20, margin=margin(r=10)),
          plot.title = element_text(size=20, margin = margin(l = 15,b=10),hjust = 0.5),
          axis.ticks.x= element_line(colour=c(rep(c(NA,"black"), t=9),NA)),
          axis.ticks.length = unit(0.3,"cm"),
          axis.ticks=element_line(linewidth =1),
          panel.grid.major.x=element_line(colour=c(rep(c(NA, "grey94"), t=9),NA)),
          panel.grid.minor=element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),      # Increase legend title size
          legend.text = element_text(size = 20),       # Increase legend text size
          legend.key.size = unit(1.2, "cm"))+
    guides(fill=guide_legend(nrow=1))
  
  # if (i %in% c(1,4,7,10)){
  # plot_rt <- plot_rt + ylab(expression(R[t]))+
  #     theme(axis.title.y = element_text(size = 20, margin = margin(r=15)))
  # }else {
  #   plot_rt <- plot_rt + theme(axis.title.y = element_blank())
  # }

  # plot(plot_rt)
  assign(paste0("rt",i),plot_rt)
}

legend_rt <- get_legend(plot_rt)
plot_rt_list <-  list(rt1,rt2,rt3,rt4,rt5,rt6,rt7,rt8,rt9,rt10,rt11,rt12)


for (i in c(1,4,7,10)){
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
  
final_rt_list <- list(get(paste0("three_rt",1)),get(paste0("three_rt",2)),get(paste0("three_rt",3)),get(paste0("three_rt",4)))
for (m in 1:length(final_rt_list)){
  final_rt_list[[m]] <- final_rt_list[[m]] + theme(legend.position = "none")
}

final_rt <- plot_grid(plot_grid(plotlist = final_rt_list,  nrow = 4, ncol = 1,rel_heights = c(1,1,1), align = "hv",axis = "tblr"),
                      legend_rt, nrow=2, rel_heights = c(2.5,0.15))
final_rt <- final_rt + theme(plot.background = element_rect(fill = "white", color = NA))
# final_rt <- ggdraw() +
#   draw_plot(final_rt, x = 0.01, y = 0, width = 0.97, height = 0.95) +
#   draw_label(expression("Estimated"~R[t]~ "over LTLAs of North East region"), fontface = 'bold', x = 0.5, y = 0.98, hjust = 0.5, size = 20)

plot(final_rt)
ggsave(filename = paste0("figures/estimated_rt_ne_jan.png"), plot = final_rt, width=18, height=15, units="in")




#----------------------------------------------------------------------------------------------------
#------ fitted death data ---------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
for (i in 1:M_regions){
  
  fit <- fit_disconnected
  est_inf_disc <- fit$draws("weekly_deaths",format="matrix")
  data_est_inf_disc <- data.frame(est_inf_disc_mean = colMeans(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)]),
                               inf_disc_min_1 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.025),
                               inf_disc_max_1 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.975),
                               inf_disc_min_2 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.25),
                               inf_disc_max_2 = colQuantiles(est_inf_disc[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.75),
                               time = seq(from=inf_start_date ,to =  end_date, by = "week"))

data_est_inf_disc <- data_est_inf_disc %>% filter(time >= fitting_start & time <= plot_end)
# 
data_inf_disc_95 <- data.frame(time = data_est_inf_disc$time, inf_disc_min = data_est_inf_disc$inf_disc_min_1,
                               inf_disc_max = data_est_inf_disc$inf_disc_max_1, key = rep("nintyfive", length(data_est_inf_disc$time)))

# data_inf_disc_50 <- data.frame(time = data_est_inf_disc$time, inf_disc_min = data_est_inf_disc$inf_disc_min_2,
#                                inf_disc_max = data_est_inf_disc$inf_disc_max_2, key = rep("fifty", length(data_est_inf_disc$time)))
# 
# # data_inf_disc <- rbind(data_inf_disc_95, data_inf_disc_50)
# # levels(data_inf_disc$key) <- c("ninetyfive", "fifty")
data_inf_disc <- data_inf_disc_95
data_inf_disc$key1 <- "95% CI of fitted death\n for disconnected model"


fit <- fit_connected
est_inf_con <- fit$draws("weekly_deaths",format="matrix")
data_est_inf_con <- data.frame(est_inf_con_mean = colMeans(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)]),
                               inf_con_min_1 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.025),
                               inf_con_max_1 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.975),
                               inf_con_min_2 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.25),
                               inf_con_max_2 = colQuantiles(est_inf_con[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.75),
                               time = seq(from=inf_start_date ,to =  end_date, by = "week"))

data_est_inf_con <- data_est_inf_con %>% filter(time >= fitting_start & time <= plot_end)
# 
data_inf_con_95 <- data.frame(time = data_est_inf_con$time, inf_con_min = data_est_inf_con$inf_con_min_1,
                              inf_con_max = data_est_inf_con$inf_con_max_1, key = rep("nintyfive", length(data_est_inf_con$time)))
# 
# data_inf_con_50 <- data.frame(time = data_est_inf_con$time, inf_con_min = data_est_inf_con$inf_con_min_2,
#                               inf_con_max = data_est_inf_con$inf_con_max_2, key = rep("fifty", length(data_est_inf_con$time)))
# 
 # data_inf_con <- rbind(data_inf_95, data_inf_50)
# #  levels(data_inf_con$key) <- c("ninetyfive", "fifty")
data_inf_con <- data_inf_con_95
data_inf_con$key1 <- "95% CI of fitted death\n for connected model"

death_regions$column_to_plot <- death_regions[[i]]
death_regions <- death_regions %>% filter(time >= fitting_start & time <= plot_end)

colors_death <- c("Fitted deaths for \ndicconected model" = "red4", "Fitted deaths for \nconnected model" = "green4", "Weekly deaths \ndata"="coral3")

breaks = sort(c(seq(ymd("2020-4-1"),ymd("2020-12-31"),by="months"),seq(ymd("2020-3-15"),ymd("2020-12-31"),by="months"),ymd("2021-1-1")))

labels = unique(date_format("%b")(data_est_Rt_con$time))
labels = as.vector(rbind(labels,rep("",length(labels))))

fit_death <-ggplot(data_est_inf_disc)+
#
  geom_ribbon(data = data_inf_disc, aes(x = time, ymin = inf_disc_min, ymax = inf_disc_max, fill=key1))+
  geom_ribbon(data = data_inf_con, aes(x = time, ymin = inf_con_min, ymax = inf_con_max, fill=key1))+
  geom_line(data = data_est_inf_disc, aes(x = time,y = est_inf_disc_mean, color = "Fitted deaths for \ndicconected model"), linewidth = 1)+
  geom_line(data = data_est_inf_con, aes(x = time, y = est_inf_con_mean, color = "Fitted deaths for \nconnected model"), linewidth = 1)+
  geom_point(data = death_regions, aes(x = time, y = column_to_plot, color = "Weekly deaths \ndata")) +

  xlab("")+
  ylab("")+
  scale_fill_manual(name = "",
                    values = c("95% CI of fitted death\n for disconnected model" = alpha("red4", 0.25),
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
  assign(paste0("death",i),fit_death)
}

legend_death <- get_legend(fit_death)
plot_death_list <-  list(death1,death2,death3,death4,death5,death6,death7,death8,death9,death10,death11,death12)


for (i in c(1,4,7,10)){
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

final_death_list <- list(get(paste0("three_deaths",1)),get(paste0("three_deaths",2)),get(paste0("three_deaths",3)), get(paste0("three_deaths",4)))
for (m in 1:length(final_death_list)){
  final_death_list[[m]] <- final_death_list[[m]] + theme(legend.position = "none")
}

final_death <- plot_grid(plot_grid(plotlist = final_death_list,  nrow = 4, ncol = 1,rel_heights = c(1,1,1), align = "hv",axis = "tblr"),
                        legend_death, nrow=2, rel_heights = c(2.5,0.15))
final_death <- final_death + theme(plot.background = element_rect(fill = "white", color = NA))

print(final_death)
ggsave(filename = paste0("figures/estimated_death_ne.png"), plot = final_death, width=18, height=15, units="in")

#################################################################################################################################
#----- plot original weekly death data for North East region ---------------------------------------------------------------

death_data <- read_excel("data/death_20_21.xlsx")
load("data/final_pop_2020_ltla.Rdata") 
#------- regions index -----------------------------------------------------------------

pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
.ltla.ne
pop_2020$area_name <- sapply(pop_2020$area_name, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

colnames(death_data) <- sapply(colnames(death_data), function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

death_north_east <- death_data %>% select(all_of(pop_2020$area_name[pop_2020$region == "North east"]))  
regions <- colnames(death_north_east)
death_north_east$week <- seq(ymd(20200101),ymd(20210131),by="week")

first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")
second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 
#---- facet plot for all the death data -------------------------------------------------------------
death_ne_long <-
  death_north_east %>%
  pivot_longer(cols = !week,names_to = "region", values_to = "deaths") %>%
  mutate('Weekly deaths' = as.integer(deaths)) 

death_ne_long$region <- factor(death_ne_long$region, 
                                   levels = regions)

breaks = sort(c(seq(ymd("2020-1-1"),ymd("2021-1-1"),by="months"),seq(ymd("2020-1-15"),ymd("2021-1-31"),by="months"),ymd("2021-02-1")))

labels = unique(date_format("%b")(death_north_east$week))
labels = as.vector(rbind(labels,rep("",length(labels))))
labels = c("",labels,"Jan","")

p <- ggplot(data = death_ne_long, aes(week,deaths)) +
  geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 0.8)+
  geom_line(color = "darkblue", linewidth = 1)+
  scale_x_date(labels=labels,breaks=breaks,limits = c(ymd("20200101"),ymd("20210201"))) + 
  facet_grid(rows = vars(region))+facet_wrap(~region, nrow=4, ncol=3)+
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
ggsave(filename = paste0("figures/death_data_ne.png"), plot = p, width=17, height=13, units="in")


##################################################################################################################################
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



##### ###############################################################################################################################
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
