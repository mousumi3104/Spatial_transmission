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
library(ISOweek)
#library(rstanarm)
library(this.path)
library(readxl)
#
script_directory <- this.path::this.dir()
setwd(script_directory)

source("data/data_arrangements_ne.R")
stan_data_separate <- stan_data_arrangements(death_threshold = 4, script_directory)

plot_required_date <- list(inf_start_date = stan_data_separate$inf_start_date,
                           fitting_start_date = stan_data_separate$fitting_start_date,
                           end_date = stan_data_separate$end_date)
plot_end <- as.Date("2020-12-31", format = "%Y-%m-%d")

stan_data_separate$inf_start_date <- NULL
stan_data_separate$fitting_start_date <- NULL
stan_data_separate$end_date <- NULL

M_regions <- stan_data_separate$M_regions

m <- cmdstan_model("fitting_regions_separate.stan")

stan_data <- list(M_regions = stan_data_separate$M_regions,
                   final_time= stan_data_separate$final_time,
                   W = stan_data_separate$W,
                   gmobility = stan_data_separate$gm_non_res,
                   initial_seeding_day = stan_data_separate$initial_seeding_day,
                   death_data_length = stan_data_separate$death_data_length,
                   death= stan_data_separate$death,
                   SI=stan_data_separate$SI,
                   f=stan_data_separate$f,
                   pop=stan_data_separate$pop,
                   day_week_index = stan_data_separate$day_week_index,
                   I = as.matrix(stan_data_separate$I),
                   fitting_start = stan_data_separate$fitting_death_start
                   )

 fit_separated = m$sample(
   data=stan_data,
   iter_sampling = 50,
   iter_warmup =100, 
   parallel_chains = 4,
   # threads_per_chain = 2,
   chains=4,
   thin=1, 
   seed=1234,
   refresh = 40,
   adapt_delta = 0.9, 
   max_treedepth = 12,init = \() list(mu = rep(3.28,M_regions), 
                                      initial_seeding = rep(5,M_regions),
                                      tau = 0.01,
                                      phi = 20,
                                      gamma = 0.5))       # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution

   Rt <- fit_separated$draws("Rt",format="matrix")
   inf <- fit_separated$draws("infection",format="matrix")
   death <- fit_separated$draws("weekly_deaths",format="matrix")
   summary_fit_separated <- fit_separated$summary()
 
 save(fit_separated ,plot_required_date,file = paste0("results/ltla_complete_disconnected_ne_jan.Rdata"))


####----------------------------------------------------------------------------------------------------------------------------
#---- plot joint model and disjoint model --------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------


load("results/ltla_disconnected_ne_jan.Rdata")
load("results/ltla_complete_disconnected_ne_jan.Rdata")

inf_start_date <- plot_required_date$inf_start_date   #as.Date("10-02-2020",format = "%d-%m-%Y")
fitting_start <- plot_required_date$fitting_start_date        #as.Date("23-03-2020", format = "%d-%m-%Y")
end_date <- plot_required_date$end_date                                           #as.Date("31-12-2020", format = "%d-%m-%Y")
plot_end <- as.Date("2020-12-31", format = "%Y-%m-%d")
first_lockdown_start <- as.Date("2020-03-23", format = "%Y-%m-%d")
first_lockdown_end <- as.Date("2020-05-10", format = "%Y-%m-%d")

second_lockdown_start <- as.Date("2020-11-05", format = "%Y-%m-%d")  
second_lockdown_end <- as.Date("2020-12-02", format = "%Y-%m-%d") 
Rt <- fit_connected$draws("Rt",format = "matrix")

M_regions <- stan_data_disconnected$M_regions

death_regions <- data.frame(stan_data_disconnected$death)
death_data_length <- stan_data_disconnected$death_data_length
death_regions$time <- seq(from = inf_start_date ,to =  end_date, by = "week")
final_time <- stan_data_disconnected$final_time

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
  data_Rt_disc$key1 <- "95% CI of \njoint model Rt"
  
  #------- separated model --------------------------------------------------------------------------------------------
  
  fit <- fit_separated
  est_Rt_sep <- fit$draws("Rt",format="matrix")
  data_est_Rt_sep <- data.frame(est_Rt_sep_mean = colMeans(est_Rt_sep[,(((i-1)*final_time)+1):(i*final_time)]),
                                Rt_sep_min_1 = colQuantiles(est_Rt_sep[,(((i-1)*final_time)+1):(i*final_time)],prob=0.025),
                                Rt_sep_max_1 = colQuantiles(est_Rt_sep[,(((i-1)*final_time)+1):(i*final_time)],prob=0.975),
                                Rt_sep_min_2 = colQuantiles(est_Rt_sep[,(((i-1)*final_time)+1):(i*final_time)],prob=0.25),
                                Rt_sep_max_2 = colQuantiles(est_Rt_sep[,(((i-1)*final_time)+1):(i*final_time)],prob=0.75),
                                time = seq(from=inf_start_date ,to =  end_date, by = "day"))
  data_est_Rt_sep <- data_est_Rt_sep %>% filter(time >= fitting_start & time <= plot_end)
  
  data_Rt_sep_95 <- data.frame(time = data_est_Rt_sep$time, Rt_sep_min = data_est_Rt_sep$Rt_sep_min_1,
                               Rt_sep_max = data_est_Rt_sep$Rt_sep_max_1, key = rep("nintyfive", length(data_est_Rt_sep$time)))
  
  data_Rt_sep_50 <- data.frame(time = data_est_Rt_sep$time, Rt_sep_min = data_est_Rt_sep$Rt_sep_min_2,
                               Rt_sep_max = data_est_Rt_sep$Rt_sep_max_2, key = rep("fifty", length(data_est_Rt_sep$time)))
  
  #  data_Rt_sep <- rbind(data_Rt_sep_95, data_Rt_sep_50)
  #  levels(data_Rt_sep$key) <- c("ninetyfive", "fifty")
  data_Rt_sep <- data_Rt_sep_95
  data_Rt_sep$key1 <- "95% CI of \nseparated model Rt"
  
  #----------------------------------------------------------------------------------------------------------------
  Rt_threshold <- data.frame(time = data_est_Rt_disc$time, Rt = rep(1,length(data_est_Rt_disc$time)))  # for Rt threshold horizontal line
  
  #---------- plot ------------------------------------------------------------------------------------------------
  
  breaks = sort(c(seq(ymd("2020-4-1"),ymd("2020-12-31"),by="months"),seq(ymd("2020-3-15"),ymd("2020-12-31"),by="months"),ymd("2021-1-1")))
  
  labels = unique(date_format("%b")(data_est_Rt_sep$time))
  labels = as.vector(rbind(labels,rep("",length(labels))))
  
  colors_rt <- c("Separated model Rt" = "skyblue3", "Joint model Rt" = "#b2182b")#, "Simulated Rt"="black")
  
  plot_rt <- ggplot(data_est_Rt_disc)+
    geom_ribbon(data = data_Rt_disc, aes(x = time, ymin = Rt_disc_min, ymax = Rt_disc_max, fill=key1))+
    geom_ribbon(data = data_Rt_sep, aes(x = time, ymin = Rt_sep_min, ymax = Rt_sep_max, fill=key1))+
    geom_line(data = data_est_Rt_disc, aes(x = time,y = est_Rt_disc_mean, color = "Joint model Rt"), linewidth = 1.2)+
    geom_line(data = data_est_Rt_sep, aes(x = time, y = est_Rt_sep_mean, color = "Separated model Rt"), linewidth = 1.2)+
    geom_line(data = Rt_threshold, aes(x=time, y = Rt),color = "black")+
    geom_vline(xintercept = as.Date(c(first_lockdown_start,first_lockdown_end,second_lockdown_start,second_lockdown_end)), linetype = "dashed", color = "black", linewidth = 1)+
    xlab("")+
    ylab("")+
    scale_fill_manual(name = "",
                      values = c("95% CI of \nseparated model Rt" = alpha("skyblue3", 0.25),
                                 "95% CI of \njoint model Rt" = alpha("#b2182b", 0.25))) +
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
  assign(paste0("rt",i),plot_rt)
  plot(plot_rt)
  
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

final_rt <- ggdraw() +
  draw_plot(final_rt, x=0, y = -0.002, height = 0.96,width=1)+
  draw_label("Estimated Rt from joint and separated model", fontface = "bold", size = 20, hjust = 0.5,color = "black",x = 0.5, y=0.98,angle = 0) 
  

final_rt <- final_rt + theme(plot.background = element_rect(fill = "white", color = NA))
plot(final_rt)
# final_rt <- ggdraw() +
#   draw_plot(final_rt, x = 0.01, y = 0, width = 0.97, height = 0.95) +
#   draw_label(expression("Estimated"~R[t]~ "over LTLAs of North East region"), fontface = 'bold', x = 0.5, y = 0.98, hjust = 0.5, size = 20)


ggsave(filename = paste0("figures/separtaed_model_estimated_rt_ne_jan.png"), plot = final_rt, width=18, height=15, units="in")

#------------------------------------------------------------------------------------------------------------------------------
#-------------plot fitted data ------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------


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
  
  
  fit <- fit_separated
  est_inf_sep <- fit_separated$draws("weekly_deaths",format="matrix")
  data_est_inf_sep <- data.frame(est_inf_sep_mean = colMeans(est_inf_sep[,(((i-1)*death_data_length)+1):(i*death_data_length)]),
                                 inf_sep_min_1 = colQuantiles(est_inf_sep[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.025),
                                 inf_sep_max_1 = colQuantiles(est_inf_sep[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.975),
                                 inf_sep_min_2 = colQuantiles(est_inf_sep[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.25),
                                 inf_sep_max_2 = colQuantiles(est_inf_sep[,(((i-1)*death_data_length)+1):(i*death_data_length)],prob=0.75),
                                 time = seq(from=inf_start_date ,to =  end_date, by = "week"))
  
  data_est_inf_sep <- data_est_inf_sep %>% filter(time >= fitting_start & time <= plot_end)
  # 
  data_inf_sep_95 <- data.frame(time = data_est_inf_sep$time, inf_sep_min = data_est_inf_sep$inf_sep_min_1,
                                inf_sep_max = data_est_inf_sep$inf_sep_max_1, key = rep("nintyfive", length(data_est_inf_sep$time)))
  # 
  # data_inf_sep_50 <- data.frame(time = data_est_inf_sep$time, inf_sep_min = data_est_inf_sep$inf_sep_min_2,
  #                               inf_sep_max = data_est_inf_sep$inf_sep_max_2, key = rep("fifty", length(data_est_inf_sep$time)))
  # 
  # data_inf_sep <- rbind(data_inf_95, data_inf_50)
  # #  levels(data_inf_sep$key) <- c("ninetyfive", "fifty")
  data_inf_sep <- data_inf_sep_95
  data_inf_sep$key1 <- "95% CI of fitted death\n for separated model"
  
  death_regions$column_to_plot <- death_regions[[i]]
  death_regions <- death_regions %>% filter(time >= fitting_start & time <= plot_end)
  
  colors_death <- c("Fitted deaths for \ndicconected model" = "red4", "Fitted deaths for \nseparated model" = "skyblue3", "Weekly deaths \ndata"="coral3")
  
  breaks = sort(c(seq(ymd("2020-4-1"),ymd("2020-12-31"),by="months"),seq(ymd("2020-3-15"),ymd("2020-12-31"),by="months"),ymd("2021-1-1")))
  
  labels = unique(date_format("%b")(data_est_Rt_sep$time))
  labels = as.vector(rbind(labels,rep("",length(labels))))
  
  fit_death <-ggplot(data_est_inf_disc)+
    #
    geom_ribbon(data = data_inf_disc, aes(x = time, ymin = inf_disc_min, ymax = inf_disc_max, fill=key1))+
    geom_ribbon(data = data_inf_sep, aes(x = time, ymin = inf_sep_min, ymax = inf_sep_max, fill=key1))+
    geom_line(data = data_est_inf_disc, aes(x = time,y = est_inf_disc_mean, color = "Fitted deaths for \ndicconected model"), linewidth = 1)+
    geom_line(data = data_est_inf_sep, aes(x = time, y = est_inf_sep_mean, color = "Fitted deaths for \nseparated model"), linewidth = 1)+
    geom_point(data = death_regions, aes(x = time, y = column_to_plot, color = "Weekly deaths \ndata")) +
    
    xlab("")+
    ylab("")+
    scale_fill_manual(name = "",
                      values = c("95% CI of fitted death\n for disconnected model" = alpha("red4", 0.25),
                                 "95% CI of fitted death\n for separated model" = alpha("skyblue3", 0.25))) +
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















