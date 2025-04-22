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
#library(rstanarm)
library(this.path)
library(gridExtra)
library(grid)
# 
script_directory <- this.path::this.dir()
setwd(script_directory)


#--------------------------------------------------------------------------------------------------------
M <- 3
week <- 40
final_time <- 300
pop <- rep(1000000, M)#c(2000000,1000000,1500000)
seed_time <- 1
initial_seeding_day <-6

source("data/distributions.R")
dist <- distributions(final_time)
si <- dist$si
f_death <- dist$f
f_case <- dist$f_case

############################################################################################################
# Rt for generation of the simulated data in the main manuscript (fig 2)

# a <- read.csv("data/GBR-estimates.csv")
# Rt_1 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
# Rt_1 <- Rt_1 %>% arrange(date)
# Rt_1 <- Rt_1 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2022-01-01")) %>% distinct(date, .keep_all = TRUE) 
# # choose only the first option
# Rt_1 <- Rt_1 %>% group_by(date) %>% slice_min(row_number())
# 
# a <- read.csv("data/IND-estimates.csv")
# Rt_2 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
# Rt_2 <- Rt_2 %>% arrange(date)
# Rt_2 <- Rt_2 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2022-01-01")) %>% distinct(date, .keep_all = TRUE)
# Rt_2 <- Rt_2 %>% group_by(date) %>% slice_min(row_number())
# 
# a <- read.csv("data/ITA-estimates.csv")
# Rt_3 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
# Rt_3 <- Rt_3 %>% arrange(date)
# Rt_3 <- Rt_3 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2022-01-01")) %>% distinct(date, .keep_all = TRUE)
# Rt_3 <- Rt_3 %>% group_by(date) %>% slice_min(row_number())

# Rt <- cbind(Rt_1 = 0.1 + Rt_1$Rt,Rt_2 = Rt_2$Rt ,Rt_3=Rt_3$Rt)
# Rt <- data.frame(Rt[1:final_time,])

##################################################################################################################

rt <- data.frame(row1 = c(0.5,1.5,1.5), row2 = c(1.5,1.5,1.5), row3 = c(1.5,1.5,1.5)) #for supplementary Fig A
C1 <- matrix(c(0.9,0.05,0.05,0.07,0.85,0.08,0.21,0.09,0.7),nrow=3,ncol=3)   # mobility matrix 
C2 <- diag(M)     # for disconnected region

m <- cmdstan_model("simulated_region.stan")   # compiling the model

for (row in 1:3){
  
    rt_row <- rt[[paste0("row",row)]]
    Rt <- data.frame(Rt_1 = rep(rt_row[1], final_time), Rt_2 = rep(rt_row[2],final_time), Rt_3 = rep(rt_row[3], final_time))
    C <- if (row %in% c(1,2)) C1 else C2

    stan_data <- list(M = M,
                      pop=pop,
                      final_time=final_time,
                      initial_seeding_day = initial_seeding_day,
                      init_seed= c(10,10,10),
                      SI=si,
                      Rt=as.matrix(Rt),
                      C=C)

    fit <- m$sample(
      data =stan_data, #_single_region,# stan_data,
      iter_sampling = 1,
      chains = 1,
      thin = 1, 
      fixed_param =TRUE) 
    
      # m_single_region <- cmdstan_model("simulated_data_single_region.stan")
      # out <- fit$draws(format = "matrix")
      daily_infection <- apply(fit$draws("infection",format = "matrix"),2,mean)

      simulated_data <- list(No_of_regions=M,
                       population = pop,
                       final_time = final_time,
                       # week = week,
                       mobility = C,
                       infection = daily_infection,
                       Rt = Rt,
                       init_seed = init_seed,
                       initial_seeding_day = initial_seeding_day)



#save data for simulated infection
# saveRDS(simulated_data, file = "data/updated_si/simulated_infection.rds")
# ---------------------------------------------------------------------------------------------------------------

     daily_inf_in_own <- apply(fit$draws("infection_in_own",format = "matrix"),2,mean)
     daily_inf_in_mob <- apply(fit$draws("infection_in_mob",format = "matrix"),2,mean)
     daily_inf_out_mob <- apply(fit$draws("infection_out_mob",format = "matrix"),2,mean)

     daily_inf <- data.frame(Region1 = daily_infection[1:final_time],
                        Region2 = daily_infection[(final_time+1):(2*final_time)],
                        Region3 = daily_infection[(2*final_time+1):(3*final_time)],
                        index = 1:final_time)

     inf_in_own <- data.frame(region1 = daily_inf_in_own[1:final_time],
                        region2 = daily_inf_in_own[(final_time+1):(2*final_time)],
                        region3 = daily_inf_in_own[(2*final_time+1):(3*final_time)],
                        index = 1:final_time)

     inf_in_mob <- data.frame(region1 = daily_inf_in_mob[1:final_time],
                         region2 = daily_inf_in_mob[(final_time+1):(2*final_time)],
                         region3 = daily_inf_in_mob[(2*final_time+1):(3*final_time)],
                         index = 1:final_time)

     inf_out_mob <- data.frame(region1 = daily_inf_out_mob[1:final_time],
                          region2 = daily_inf_out_mob[(final_time+1):(2*final_time)],
                          region3 = daily_inf_out_mob[(2*final_time+1):(3*final_time)],
                          index = 1:final_time)

     daily_inf_long <- daily_inf %>% pivot_longer(cols = starts_with("Region"),names_to = "variable", values_to = "value")
     inf_in_own_long <- inf_in_own %>% pivot_longer(cols = starts_with("Region"),names_to = "variable", values_to = "value")
     inf_in_mob_long <- inf_in_mob %>% pivot_longer(cols = starts_with("Region"),names_to = "variable", values_to = "value")
     inf_out_mob_long <- inf_out_mob %>% pivot_longer(cols = starts_with("Region"),names_to = "variable", values_to = "value")

     print(sum(daily_infection))
 
     data_stack_plot1 <-  data.frame(in_mob = inf_in_mob$region1, in_own = inf_in_own$region1, out_mob = inf_out_mob$region1, time = inf_out_mob$index)

     data_stack_plot_1_final <- data.frame(time = rep(data_stack_plot1$time, times = 3),
                                           types = factor(rep(c("in_own", "in_mob", "out_mob"), each = length(data_stack_plot1$time)),levels = c("in_own", "in_mob", "out_mob")),
                                          infections = c(data_stack_plot1$in_own, data_stack_plot1$in_mob, data_stack_plot1$out_mob))
     
     data_stack_plot2 <-  data.frame(in_mob = inf_in_mob$region2, in_own = inf_in_own$region2, out_mob = inf_out_mob$region2, time = inf_out_mob$index)

     data_stack_plot_2_final <- data.frame(time = rep(data_stack_plot2$time, times = 3),
                                           types = factor(rep(c("in_own", "in_mob", "out_mob"), each = length(data_stack_plot2$time)),levels = c("in_own", "in_mob", "out_mob")),
                                           infections = c(data_stack_plot2$in_own, data_stack_plot2$in_mob, data_stack_plot2$out_mob))
     
      data_stack_plot3 <-  data.frame(in_mob = inf_in_mob$region3, in_own = inf_in_own$region3, out_mob = inf_out_mob$region3, time = inf_out_mob$index)

      data_stack_plot_3_final <- data.frame(time = rep(data_stack_plot3$time, times = 3),
                                            types = factor(rep(c("in_own", "in_mob", "out_mob"), each = length(data_stack_plot3$time)),levels = c("in_own", "in_mob", "out_mob")),
                                            infections = c(data_stack_plot3$in_own, data_stack_plot3$in_mob, data_stack_plot3$out_mob))

 ######### infection #########

      p <- ggplot(daily_inf_long, aes(x = index, y = value, color = variable)) +
           geom_point(size = 3) +
           labs(title = "",
           x = "",
           y = "",
           color = "") +
           theme_bw() +
           ylim(0, ifelse(row==1, 5000, 10000)) +
           theme(axis.text.x = element_text(angle = 0, hjust = 0.4, vjust = 0.4, size = 25, margin = margin(r = 10), color = "black"),
           axis.text.y = element_text(size = 25, margin = margin(r = 10), color = "black"),
           axis.title.y = element_text(size = 25, margin = margin(r = 10)),
           axis.title.x = element_text(size = 25, margin = margin(r = 10)),
           axis.ticks.x = element_line(colour = "black"),
           axis.ticks.length = unit(0.3, "cm"),
           axis.ticks = element_line(linewidth = 1),
           panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
           panel.grid.minor = element_blank(),
           panel.border = element_rect(color = "black", linewidth = 1.5),
           plot.margin = margin(0, 14, 0, 0),
           legend.position = "bottom",
           legend.title = element_blank(),
           legend.text = element_text(size = 25),
           legend.key.size = unit(1.2, "cm")) +
           scale_color_manual(values = c("Region1" = "#8c510a", "Region2" = "gold2", "Region3" = "#01665e")) +
           guides(fill = guide_legend(ncol = 1)) 
          # Add the annotation using annotate.  This is now inside the ggplot call.
           
        legend_top <- get_legend(p)
        p <- p + theme(legend.position = "none")
        assign(paste0("row1",row),p)

  ###### infection differnt sources region 1

p1 <- ggplot(data_stack_plot_1_final, aes(x = time, y = infections, fill = types)) +
  # geom_ribbon(aes(x = time, ymin = inf_min, ymax = inf_max, fill = key), alpha =0.25,show.legend = FALSE)+
  geom_area(position = "stack") +
  # geom_point(data = daily_inf, aes(x=index,y=Region1),color = "#8c510a",inherit.aes = FALSE, size =1.5)+
  # geom_ribbon(data = plot_est_infection, aes(x = time, ymin = est_infection_min1, ymax = est_infection_max1),fill = "red4",alpha = 0.45, show.legend = FALSE,inherit.aes = FALSE)+
  # geom_line(data = plot_est_infection, aes(x=time, y=est_infection_mean, color = "fitted_infection"),inherit.aes = FALSE,linewidth=1) +
  scale_fill_manual(values = c("in_own" = alpha("#ff7f00",0.45), "in_mob" = alpha("deepskyblue4",0.45), "out_mob" = alpha("#7570b3",0.45)),
                    labels = c(TeX(r"(\overset{\normalsize{Infections driven by}}{\normalsize{own infected population}})"),
                               TeX(r"(\overset{\normalsize{Mobility induced}}{\normalsize{infections within region}})") ,
                               TeX(r"(\overset{\normalsize{Mobility induced infections}}{\normalsize{outside the region}})"))) +
  labs(
    x = "",
    y = "",
    title= ""
  ) +
  theme_bw() +
  ylim(0, ifelse(row == 1, 5000, 10000))+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 25,margin = margin(r=10),color="black"),
        axis.text.y = element_text(size = 25,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 25, margin=margin(r=10)),
        axis.title.x = element_text(size = 25, margin=margin(r=10)),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        plot.margin = margin(0, 14, 0, 0),
        legend.position = "right",
        legend.title = element_blank(),      # Increase legend title size
        legend.text = element_text(size = 25),       # Increase legend text size
        legend.key.size = unit(1.2, "cm"))+
        guides(fill=guide_legend(nrow=1))
  legend_bottom <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  # print(p1)
  assign(paste0("row2",row),p1)

p2 <- ggplot(data_stack_plot_2_final, aes(x = time, y = infections, fill = types)) +
  # geom_ribbon(aes(x = time, ymin = inf_min, ymax = inf_max, fill = key), alpha =0.25,show.legend = FALSE)+
  geom_area(position = "stack") +
  # geom_point(data = daily_inf, aes(x=index,y=Region2),color = "#d8b365",inherit.aes = FALSE, size =1.5)+
  # geom_ribbon(data = plot_est_infection, aes(x = time, ymin = est_infection_min1, ymax = est_infection_max1),fill = "red4",alpha = 0.45, show.legend = FALSE,inherit.aes = FALSE)+
  # geom_line(data = plot_est_infection, aes(x=time, y=est_infection_mean, color = "fitted_infection"),inherit.aes = FALSE,linewidth=1) +
  scale_fill_manual(values = c("in_own" = alpha("#ff7f00",0.45), "in_mob" = alpha("deepskyblue4",0.45), "out_mob" = alpha("#7570b3",0.45)),                    labels = c(TeX(r"(\overset{\normalsize{Infections driven by}}{\normalsize{own infected population}})"),
                               TeX(r"(\overset{\normalsize{Mobility induced}}{\normalsize{infections within region}})") ,
                               TeX(r"(\overset{\normalsize{Mobility induced infections}}{\normalsize{outside the region}})"))) +
  labs(
    x = "",
    y = "",
    title= ""
  ) +
  theme_bw() +
  ylim(0, ifelse(row == 1, 5000, 10000))+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 25,margin = margin(r=10),color="black"),
        axis.text.y = element_text(size = 25,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 25, margin=margin(r=10)),
        axis.title.x = element_text(size = 25, margin=margin(r=10)),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        plot.margin = margin(0, 14, 0, 0),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        legend.position = "none")
 # print(p2)
 assign(paste0("row3",row),p2)

p3 <- ggplot(data_stack_plot_3_final, aes(x = time, y = infections, fill = types)) +
  # geom_ribbon(aes(x = time, ymin = inf_min, ymax = inf_max, fill = key), alpha =0.25,show.legend = FALSE)+
  geom_area(position = "stack") +
  # geom_point(data = daily_inf, aes(x=index,y=Region3),color = "#01665e",inherit.aes = FALSE, size =1.5)+
  # geom_ribbon(data = plot_est_infection, aes(x = time, ymin = est_infection_min1, ymax = est_infection_max1),fill = "red4",alpha = 0.45, show.legend = FALSE,inherit.aes = FALSE)+
  # geom_line(data = plot_est_infection, aes(x=time, y=est_infection_mean, color = "fitted_infection"),inherit.aes = FALSE,linewidth=1) +
  scale_fill_manual(values = c("in_own" = alpha("#ff7f00",0.45), "in_mob" = alpha("deepskyblue4",0.45), "out_mob" = alpha("#7570b3",0.45)),                    labels = c(TeX(r"(\overset{\normalsize{Infections driven by}}{\normalsize{own infected population}})"),
                               TeX(r"(\overset{\normalsize{Mobility induced}}{\normalsize{infections within region}})") ,
                               TeX(r"(\overset{\normalsize{Mobility induced infections}}{\normalsize{outside the region}})"))) +
  labs(
    x = "Day",
    y = "",
    title= ""
  ) +
  theme_bw() +
  ylim(0, ifelse(row == 1, 5000, 10000))+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.4, vjust = 0.4,size = 25,margin = margin(r=10),color="black"),
        axis.text.y = element_text(size = 25,margin = margin(r=10),color="black"),
        axis.title.y = element_text(size = 25, margin=margin(r=10)),
        axis.title.x = element_text(size = 25, margin=margin(r=10)),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        plot.margin = margin(0, 14, 0, 0),
        axis.ticks.length = unit(0.3,"cm"),
        axis.ticks = element_line(linewidth =1),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        legend.position = "none")

  # print(p3)
  assign(paste0("row4",row),p3)
}

plot_list1 <- list(row11,row12,row13)
plot_list2 <- list(row21,row22,row23,row31,row32,row33,row41,row42,row43)

index_list1 <- c("(a)","(b)","(c)")
index_list2 <- c("(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)")

p_final1 <- plot_grid(plot_grid(plotlist=plot_list1, nrow = 1, ncol = 3, rel_widths =  c(1,1,1), align = "hv",axis = "tblr",
                                labels = index_list1,
                                label_size = 25,
                                label_x = 0.35,
                                label_y = 0.92),
                               legend_top, nrow =2, rel_heights = c(0.9,0.15), align = "hv",axis = "tblr")

p_final2 <- plot_grid(plot_grid(plotlist = plot_list2, nrow = 3, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,1,1), align = "hv",axis = "tblr",
                                labels = index_list2,
                                label_size = 25,
                                label_x = 0.35,
                                label_y = 0.92),
                      legend_bottom, nrow =2, rel_heights = c(0.9,0.1), align = "hv",axis = "tblr")

joint12 <- plot_grid(p_final1,p_final2, nrow=2, rel_heights = c(1,3),align = "hv",axis = "tblr")+
  draw_label("Three sources of infections", x = 0.015, y = 0.5, angle = 90, size = 27) +
  draw_label("Daily infections", x = 0.015, y = 0.91, angle = 90, size = 27) 
  

label_box <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           fill = "white", color = "black", linewidth = 1) +
  theme_void()

p_final <- ggdraw() +
  draw_plot(joint12, x = 0.01, y = 0, width = 0.94, height = 0.98)  +
  draw_plot(label_box, x = 0.938, y = 0.12, width = 0.043, height = 0.165) +
  draw_label("Region 3", x = 0.96, y = 0.20, angle = 270, size = 26,fontface = "bold") +
  draw_plot(label_box, x = 0.938, y = 0.34, width = 0.043, height = 0.165) +
  draw_label("Region 2", x = 0.96, y = 0.42, angle = 270, size = 26,fontface = "bold") +
  draw_plot(label_box, x = 0.938, y = 0.56, width = 0.043, height = 0.165) +
  draw_label("Region 1", x = 0.96, y = 0.638, angle = 270, size = 26,fontface = "bold") +
  draw_plot(label_box, x = 0.737, y = 0.964, width = 0.22, height = 0.035) +
  draw_label(expression(paste(R[t]^1,"=1.5",", ",R[t]^2,"=1.5",", ",R[t]^3,"=",1.5)),
             x = 0.846, y = 0.981, angle = 0, size = 20, fontface = "bold") +
  draw_plot(label_box, x = 0.424, y = 0.964, width = 0.22, height = 0.035) +
  draw_label(expression(paste(R[t]^1,"=1.5",", ",R[t]^2,"=1.5",", ",R[t]^3,"=",1.5)),
             x = 0.533, y = 0.981, angle = 0, size = 20, fontface = "bold") +
  draw_plot(label_box, x = 0.11, y = 0.964, width = 0.22, height = 0.035) +
  draw_label(expression(paste(R[t]^1,"=0.5",", ",R[t]^2,"=1.5",", ",R[t]^3,"=",1.5)),
             x = 0.219, y = 0.981, angle = 0, size = 20, fontface = "bold")

p_final <- ggdraw() +
  draw_plot(p_final, x = 0, y = 0, width = 1, height = 0.965)  +
  draw_label("With mobility", x = 0.38, y = 0.99, angle = 0, size = 26,fontface = "bold") +
  draw_label("Without mobility", x = 0.85, y = 0.99, angle = 0, size = 26,fontface = "bold") +
  draw_line(x = c(0.11, 0.11), y = c(0.963, 0.973), size = 1.2) +   # vertical bar
  draw_line(x = c(0.64, 0.64), y = c(0.963, 0.973), size = 1.2) +   # vertical bar
  draw_line(x = c(0.11, 0.64), y = c(0.973, 0.973), size = 1.2)  
  # draw_plot(mat_plot1, x = 0.3, y = 0.57, width = 0.4, height = 0.4)+
  # draw_plot(mat_plot, x = 0.55, y = 0.57, width = 0.4, height = 0.4)
  
p_final <- p_final +  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("figures/updated_si/supp_fig8.png", plot=p_final, width = 15, height = 15, dpi = 300)
browseURL("figures/updated_si/supp_fig8.png")

