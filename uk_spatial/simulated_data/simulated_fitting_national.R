library(cmdstanr)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(stringr)
library(abind)
library(bayesplot)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(cowplot)
library(this.path)
# 
script_directory <- this.path::this.dir()
setwd(script_directory)

#--------- data arrangement ------------------------------------------------------------------
week <- 7
final_time <- 50 * week
pop <- c(20000000, 10000000, 15000000)
initial_seeding_day = 6

a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/GBR-estimates.csv")
Rt_1 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_1 <- Rt_1 %>% arrange(date)
Rt_1 <- Rt_1 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2022-01-01")) %>% distinct(date, .keep_all = TRUE) 
# choose only the first option
Rt_1 <- Rt_1 %>% group_by(date) %>% slice_min(row_number())

a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/IND-estimates.csv")
Rt_2 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_2 <- Rt_2 %>% arrange(date)
Rt_2 <- Rt_2 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2022-01-01")) %>% distinct(date, .keep_all = TRUE)
Rt_2 <- Rt_2 %>% group_by(date) %>% slice_min(row_number())

a <- read.csv("~/OneDrive - National University of Singapore/uk_mobility_data/ITA-estimates.csv")
Rt_3 <- data.frame(date = as.Date(a$date), Rt = a$median_R_mean)
Rt_3 <- Rt_3 %>% arrange(date)
Rt_3 <- Rt_3 %>% filter(date >= as.Date("2020-03-01") & date <= as.Date("2022-01-01")) %>% distinct(date, .keep_all = TRUE)
Rt_3 <- Rt_3 %>% group_by(date) %>% slice_min(row_number())

Rt <- cbind(Rt_1 = 0.1 + Rt_1$Rt,Rt_2 = Rt_2$Rt ,Rt_3=Rt_3$Rt)
Rt <- data.frame(Rt[1:N,])

# day_week_index <- c(rep(1,(fitting_case_start -1)*7), rep(2:((final_time/7)),each=7))
# day_week_index <- day_week_index[1:final_time]
# week <- day_week_index[length(day_week_index)]  

#-----------distributions ----------------------------------------------------------------------
si <- rep(0,N)
si[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
for (i in 2:N){
  si[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
}

mean1 <- 5.1; cv1 <- 0.86; mean2 <- 17.8 ; cv2 <- 0.45;
x1 <- rgammaAlt(1e6, mean1, cv1)
x2 <- rgammaAlt(1e6, mean2, cv2)

f1 <- rep(0, N)
f1_cached <- ecdf(x1 + x2)
ifr <- 1
iar <- 1

convolution <- function(u) (ifr * f1_cached(u))
f1[1] <- (convolution(1.5) - convolution(0))
for (i in 2:N) {
  f1[i] <- (convolution(i + .5) - convolution(i - .5))
}

f2 <- rep(0, N)
f2_cached <- ecdf(x1)
convolution <- function(u) (f2_cached(u))
f2[1] <- (convolution(1.5) - convolution(0))
for(i in 2:N) {
  f2[i] <- (convolution(i + .5) - convolution(i - .5))
}


#---- Run the model ----------------------------------------------------------------------------------
for (i in 1:M_regions){
stan_data_single_region <- list(pop=pop[i],
                                final_time = final_time,
                                initial_seeding_day = initial_seeding_day,
                                init_seed = 30,
                                SI=si,
                                f1=f1,
                                f2=f2,
                                Rt=Rt[,i],
                                iar=iar)

m_single_region <- cmdstan_model("simulated_data_single_region.stan")

fit <- m_single_region$sample(
  data =stan_data_single_region,
  iter_sampling = 1,
  chains = 1,
  thin = 1, 
  fixed_param =TRUE)  

daily_infection <- fit$draws("infection",format = "matrix")
daily_deaths <- fit$draws("daily_deaths", format = "matrix")
daily_cases <- fit$draws("daily_cases", format = "matrix")
weekly_death <- fit$draws("weekly_deaths",format= "matrix")

plot(colMeans(daily_infection))
save(daily_infection, file = paste0("infection",i,".rds"))
}

#----- fitting with the model ------------------------------------------------------------------------------------------------

data_fit <- floor(as.numeric(daily_deaths))
data_length <- length(data_fit)
fitting_start <- 30   #which(cumsum(data_fit) >= 10)    #which(data_national>1)[1]

stan_data <- list(pop=pop,
                  final_time=final_time,
                  initial_seeding_day=initial_seeding_day,
                  data_length = data_length,
                  data_fit = data_fit,
                  SI=si,
                  f1=f1,
                  f2=f2,
                  fitting_start = fitting_start,
                  iar=iar)

m <- cmdstan_model("simulated_fitting_national.stan")

fit <-  m$sample(
  data=stan_data,
  iter_sampling = 200, 
  iter_warmup = 500, 
  parallel_chains = 4,
  chains=2,
  thin=1,
  seed=12345,
  refresh = 40,
  adapt_delta = 0.95,
  max_treedepth = 10)

summary_fit <- fit$summary()
saveRDS(fit, )
# save(fit,stan_data, file=paste0("results/",'national_fitting_cmd.Rdata'))

#--------------- plot -----------------------------------------------------------------#

estimated_infection <- fit$draws("infection",format="matrix")
estimated_inf_mean <- colMeans(estimated_infection)

estimated_deaths <- fit$draws("deaths",format="matrix")
estimated_deaths_mean <- colMeans(estimated_deaths)

estimated_cases <- fit$draws("cases",format = "matrix")
estimated_cases_mean <- colMeans(estimated_cases)

estimated_rt <- fit$draws("Rt",format="matrix")
estimated_rt_mean <- colMeans(estimated_rt)

plot(estimated_rt_mean,type="l")
points(stan_data_single_region$Rt,col="red")

plot(estimated_deaths_mean,type="l")
points(data_fit,"red")



# estimated_rt <- fit$draws("Rt",format="matrix")
# 
# estimated_rt_mean <-  colMeans(estimated_rt)
# estimated_rt_li_1 <- colQuantiles(estimated_rt,prob=0.025)
# estimated_rt_ui_1 <- colQuantiles(estimated_rt,prob=0.975)
# estimated_rt_li_2 <- colQuantiles(estimated_rt,prob=0.25)
# estimated_rt_ui_2 <- colQuantiles(estimated_rt,prob=0.75)
# 
# # estimated_initial_seed <- fit$draws("initial_seeding",format="matrix")
# # estimated_initial_seed_mean <- colMeans(estimated_initial_seed)
# 
# 
# save(estimated_rt_mean,file = paste0("data/estimated_rt_mean.Rdata"))

# data_estimated_rt <- data.frame(Rt = estimated_rt_mean,
#                                 Rt_min_1 = estimated_rt_li_1,
#                                 Rt_max_1 = estimated_rt_ui_1,
#                                 Rt_min_2 = estimated_rt_li_2,
#                                 Rt_max_2 = estimated_rt_ui_2,
#                                 time = 1:final_time)
# 
# 
# data_estimated_death <- data.frame(estimated_deaths = estimated_deaths_mean,
#                                    death_min_1 = estimated_deaths_li_1,
#                                    death_max_1 = estimated_deaths_ui_1,
#                                    death_min_2 = estimated_deaths_li_2,
#                                    death_max_2 = estimated_deaths_ui_2,
#                                    reported_death = death_national,
#                                    week = 1:(final_time/7))
# 
# 
# data_Rt_95 <- data.frame(time = data_estimated_rt$time, Rt_min = data_estimated_rt$Rt_min_1, 
#                          Rt_max = data_estimated_rt$Rt_max_1, key = rep("nintyfive", length(data_estimated_rt$time)))
# 
# data_Rt_50 <- data.frame(time = data_estimated_rt$time, Rt_min = data_estimated_rt$Rt_min_2, 
#                          Rt_max = data_estimated_rt$Rt_max_2, key = rep("fifty", length(data_estimated_rt$time)))
# 
# data_death_95 <- data.frame(time = data_estimated_death$week, death_min = data_estimated_death$death_min_1, 
#                             death_max = data_estimated_death$death_max_1, key = rep("nintyfive", length(data_estimated_death$week)))
# 
# data_death_50 <- data.frame(time = data_estimated_death$week, death_min = data_estimated_death$death_min_2, 
#                             death_max = data_estimated_death$death_max_2, key = rep("fifty", length(data_estimated_death$week)))
# 
# 
# data_Rt <- rbind(data_Rt_95, data_Rt_50)
# levels(data_Rt$key) <- c("ninetyfive", "fifty")
# 
# data_death <- rbind(data_death_95, data_death_50)
# levels(data_death$key) <- c("ninetyfive", "fifty")
# 
# data_Rt_1 <- data.frame(time = data_estimated_rt$time, Rt = rep(1,length(data_estimated_rt$time)))
# 
# p1 <- ggplot(data_estimated_rt)+
#   geom_ribbon(data = data_Rt, aes(x=time,ymin = Rt_min, ymax = Rt_max, fill=key))+
#   geom_line(data = data_estimated_rt, aes(x=time,y=Rt),color = "black",linewidth=0.8)+
#   # geom_line(data = data_Rt_1, aes(x=time, y = Rt))+
#   xlab("")+
#   ylab("Estimated Rt")+
#   # scale_x_date(date_breaks = "4 week", labels = date_format("%d %b"))+
#   scale_fill_manual(name = "", labels = c("50%", "95%"),
#                     values = c(alpha("deepskyblue4", 0.35), 
#                                alpha("deepskyblue4", 0.25))) + 
#   ggtitle("Weekly National model")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         axis.title.y = element_text(size = 14),
#         legend.position = "None")  +
#   guides(fill=guide_legend(ncol=1))
# 
# 
# 
# p2 <-ggplot(data_estimated_death)+
#   geom_ribbon(data = data_death, aes(x=time, ymin = death_min, ymax = death_max, fill=key))+
#   geom_point(data = data_estimated_death, aes(x=week,y=reported_death),color = "coral",size =2)+
#   geom_line(data = data_estimated_death, aes(x=week,y=estimated_deaths_mean),color = "black",size=0.8)+
#   xlab("")+
#   ylab("Weekly reported deaths")+
#   # scale_x_date(date_breaks = "4 week", labels = date_format("%d %b"))+
#   scale_fill_manual(name = "", labels = c("50%", "95%"),
#                     values = c(alpha("deepskyblue4", 0.35), 
#                                alpha("deepskyblue4", 0.25))) + 
#   scale_shape_manual(values = 16)+
#   ggtitle("")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         axis.title.y = element_text(size = 14),
#         legend.position = "None")  +
#   guides(fill=guide_legend(ncol=1))
# 
# plot(p2)
# 
# p <- plot_grid(p1,p2,ncol=1,nrow=2,rel_heights = c(1,1))
# plot(p)
# 
# filename = paste0("figures/national_model",".png")
# # ggsave(filename, plot = p, width = 8, height = 8, units = "in")
