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
library(this.path)
library(ISOweek)
# 
script_directory <- this.path::this.dir()
setwd(script_directory)
#-------------- loading data -------------------------------------------------------------------------
# data_google_mobility <- read.csv("data/2020_GB_Region_Mobility_Report.csv")
# 
# mobility_change <-  data_google_mobility %>% filter(sub_region_1 %in% c("Bath and North East Somerset","Bedford","Blackburn with Darwen","Blackpool","Borough of Halton","Bracknell Forest","Brighton and Hove",
# "Bristol City","Buckinghamshire","Cambridgeshire", "Central Bedfordshire","Cornwall","Derby","Gloucestershire","Hartlepool","Leicester","Medway","Middlesbrough","North Yorkshire","Plymouth","Shropshire","Tyne and Wear",
#  "Wiltshire","York","Windsor and Maidenhead","Swindon","Stockton-on-Tees","South Yorkshire","Slough","Portsmouth","Oxfordshire","Northamptonshire","North East Lincolnshire","Merseyside","Leicestershire","Isle of Wight",
#  "Herefordshire","Greater London","Greater Manchester","Hertfordshire","Kent","Lincolnshire","Milton Keynes","Northumberland","Rutland","Somerset","Southampton","Stoke-on-Trent","Thurrock","Warrington",
#   "Wokingham","Worcestershire","West Sussex","Warwickshire","Torbay","Suffolk","Southend-on-Sea","Reading","Nottingham","North Lincolnshire","Kingston upon Hull","Dorset","Essex","Hampshire","Lancashire","Luton","Norfolk",
#   "North Somerset","Nottinghamshire" ,"Peterborough","Redcar and Cleveland","South Gloucestershire","Staffordshire","Surrey","West Berkshire","West Yorkshire"))
# 
# colnames(mobility_change)[colnames(mobility_change) == "sub_region_1"] <- "area_name" 
# mobility_change <- mobility_change %>% select(-c("country_region_code","country_region","sub_region_2","metro_area","iso_3166_2_code","census_fips_code","place_id"))
# 
# mobility_change <- merge(mobility_change, pop_2020[,c("area_name","region")], by = "area_name", all.x = TRUE)
# mobility_change <- mobility_change %>%
#   left_join(a, by = "area_name", suffix = c("_1", "_2")) %>%
#   mutate(
#     region = coalesce(region_1, region_2)
#   ) %>%
#   select(-region_1,-region_2)
# 
# mobility_change$region <- str_to_sentence(tolower(mobility_change$region))


load("data/final_pop_2020_ltla.Rdata")
load("data/england_death_2020.Rdata")       ## weekly data
load("data/uk_regions_mobility_matrix.Rdata")

# load("~/mousumi_codes/uk_spatial/data/final_pop_2020_ltla.Rdata")
# load("~/mousumi_codes/uk_spatial/data/england_death_2020.Rdata")
# load("~/mousumi_codes/uk_spatial/data/uk_ltla_mobility_matrix.Rdata")

#--------------------------------------------------------------------------------------

pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

death_data <- death_data %>% select(all_of(pop_2020$area_name))  

#------- region_pop ----------------------------------------------------------------- 

pop_region <- data.frame(regions = unique(pop_2020$region))

population_region <- rep(0,nrow(pop_region))
for (i in 1:nrow(pop_region)){
    population_region[i] = sum(pop_2020$population[pop_2020$region == pop_region$regions[i]])
}

pop_region$population <- population_region

#------- regions index -----------------------------------------------------------------

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

death_regions <- data.frame(north_east = apply(death_data[,north_east_index],1,sum),
                            north_west = apply(death_data[,north_west_index],1,sum),
                            yorkshire = apply(death_data[,yorkshire_index],1,sum),
                            east_midlands = apply(death_data[,east_midlands_index],1,sum),
                            west_midlands = apply(death_data[,west_midlands_index],1,sum),
                            east = apply(death_data[,east_index],1,sum),
                            london = apply(death_data[,london_index],1,sum),
                            south_east = apply(death_data[,south_east_index],1,sum),
                            south_west = apply(death_data[,south_west_index],1,sum))

#--------- start date of the epidemic ----------------------------------------------------------

death_regions$total_death <- apply(death_data,1, sum)
fitting_death_start <- which(cumsum(death_regions$total_death) > 10 )[1]       # week when cumulative death exceeds above 10
epidemic_start <- fitting_death_start-4
inf_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", epidemic_start), "-1"))
fitting_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", fitting_death_start), "-1"))
inf_end_date <- as.Date("31-12-2020",format = "%d-%m-%Y")
week_dates <- seq.Date(inf_start_date, inf_end_date, by = "week")
week_numbers <- length(week_dates)

#--------- data arrangement ------------------------------------------------------------------

final_time <- as.numeric(inf_end_date - inf_start_date) +1
M_regions <-ncol(death_regions)-1     # number of region

C_base <- mobmatrix_region_norm   # mobility matrix 
C_lockdown <- diag(M_regions)

initial_seeding_day = 6
pop = pop_region$population[1:M_regions]          #pop_2020$population[1:M_regions]

death_regions <- death_regions[epidemic_start:nrow(death_regions),1:M_regions]   #(adjust accordingly)
death_len_data <- nrow(death_regions)


#--------------- google mobility -----------------------------------------------------------------------------
mobility_change <- readRDS("data/mobility_change.rds")
mobility_change <- mobility_change %>%
  group_by(area_name, date, region) %>%
  summarize(avg_retail_and_recreation_percent_change_from_baseline = mean(retail_and_recreation_percent_change_from_baseline, na.rm = TRUE),
            ave_grocery_and_pharmacy_percent_change_from_baseline = mean(grocery_and_pharmacy_percent_change_from_baseline, na.rm = TRUE),
            ave_parks_percent_change_from_baseline = mean(parks_percent_change_from_baseline, na.rm = TRUE),
            ave_transit_stations_percent_change_from_baseline = mean(transit_stations_percent_change_from_baseline, na.rm = TRUE),
            ave_workplaces_percent_change_from_baseline = mean(workplaces_percent_change_from_baseline, na.rm = TRUE))

mobility_change_region <- mobility_change %>%
  group_by(region, date) %>%
  summarize(avg_retail_and_recreation_percent_change_from_baseline = mean(avg_retail_and_recreation_percent_change_from_baseline, na.rm = TRUE),
            ave_grocery_and_pharmacy_percent_change_from_baseline = mean(ave_grocery_and_pharmacy_percent_change_from_baseline, na.rm = TRUE),
            ave_parks_percent_change_from_baseline = mean(ave_parks_percent_change_from_baseline, na.rm = TRUE),
            ave_transit_stations_percent_change_from_baseline = mean(ave_transit_stations_percent_change_from_baseline, na.rm = TRUE),
            ave_workplaces_percent_change_from_baseline = mean(ave_workplaces_percent_change_from_baseline, na.rm = TRUE),.groups = "keep")

mobility_change_region$date <- as.Date(mobility_change_region$date, format = "%Y-%m-%d")

if (min(mobility_change_region$date) > inf_start_date){
  new_dates <- seq.Date(from = inf_start_date, to = min(mobility_change_region$date), by = "day")
  new_rows <- expand.grid(region = unique(mobility_change_region$region), date = new_dates)
  new_rows <- new_rows %>% mutate(avg_retail_and_recreation_percent_change_from_baseline = 0, 
                                  ave_grocery_and_pharmacy_percent_change_from_baseline = 0, 
                                  ave_parks_percent_change_from_baseline = 0, 
                                  ave_transit_stations_percent_change_from_baseline = 0,
                                  ave_workplaces_percent_change_from_baseline = 0)  
}
mobility_change_region <- rbind(new_rows, mobility_change_region)
desired_order_region <- pop_region$regions
mobility_change_region$region <- factor(mobility_change_region$region, levels = desired_order_region)
mobility_change_region <- mobility_change_region %>% arrange(region)

gmobility <- array(NA, dim = c(final_time+1,5,M_regions))
for (i in 1:M_regions){
  gmobility[,,i] <- as.matrix(mobility_change_region[(((i-1)*(final_time+1))+1):(i*(final_time+1)),3:7])
}

for (m in 1:M_regions) {
  for (i in 1:4){
    for (t in 1:(final_time+1)){
      if (is.nan(gmobility[t,i,m])){
       
        prev_val <- if (t > 1) gmobility[max(which(!is.nan(gmobility[1:(t-1),i,m] ))),i,m]
        next_val <- if (t < final_time) gmobility[(min(which(!is.nan(gmobility[(t+1):final_time,i,m]))) + t),i,m]
        gmobility[t,i,m] <- (prev_val +next_val)/2
      }
    }
  }
}


#------------- lockdown effect --------------------------------------------------------------------------------

lockdown1_started <- as.Date("2020-03-23", format = "%Y-%m-%d")
lockdown1_lifted <- as.Date("2020-05-10", format = "%Y-%m-%d")

lockdown2_started <- as.Date("2020-11-05", format = "%Y-%m-%d")  
lockdown2_lifted <- as.Date("2020-12-02", format = "%Y-%m-%d") 

lockdown_index <-  data.frame(date = seq.Date(from = inf_start_date,as.Date("2020-12-31",format = "%Y-%m-%d"),by="day"))
lockdown_index$r_mobility_index <- rep(0,nrow(lockdown_index))
lockdown_index$r_mobility_index[lockdown_index$date >= lockdown1_started & lockdown_index$date <= lockdown1_lifted] <- 1
lockdown_index$r_mobility_index[lockdown_index$date >= lockdown2_started & lockdown_index$date <= lockdown2_lifted] <- 1

lockdown_index$g_mobility_index <-  rep(3,nrow(lockdown_index))
lockdown_index$g_mobility_index[lockdown_index$date <= lockdown1_lifted] <- 1
lockdown_index$g_mobility_index[lockdown_index$date >= lockdown1_lifted & lockdown_index$date <= lockdown2_lifted] <- 2

lockdown_index$L1 <- ifelse(lockdown_index$g_mobility_index == 1,1,0)
lockdown_index$L2 <- ifelse(lockdown_index$g_mobility_index == 2,1,0)
lockdown_index$L3 <- ifelse(lockdown_index$g_mobility_index == 3,1,0)

lockdown_index <- lockdown_index %>% select(- g_mobility_index, -date)

#-----------------------------------------------------------------------------------------------
day_week_index <- array(0,final_time)
for (t in 1:final_time){ 
  day_week_index[t] = ceiling(t/7) -(fitting_death_start - epidemic_start - 1)
}
day_week_index[1:((fitting_death_start-epidemic_start)*7)] = 1
week <- day_week_index[length(day_week_index)]  

#-----------------------------------------------------------------------------------------------
si <- rep(0,final_time)
si[1] = integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=0, upper=1.5)$value
for (i in 2:final_time){
  si[i] <- integrate(function(x) dgamma(x,shape=6.5, rate=0.62), lower=i-0.5, upper=i+0.5)$value
}

mean1 <- 5.1; cv1 <- 0.86; mean2 <-17.8 ; cv2 <- 0.45;
x1 <- rgammaAlt(1e6,mean1,cv1)
x2 <- rgammaAlt(1e6,mean2,cv2)
f <- rep(0,final_time)

f_cached <- ecdf(x1+x2)

convolution <- function(u) (0.0103*f_cached(u))      # ifr is 0.0103.  # this is the ifr for uk from the code of swapnil's NPI nature
f[1] = (convolution(1.5) - convolution(0))
for(i in 2:final_time) {
  f[i] = (convolution(i+.5) - convolution(i-.5)) 
}

f_case <- rep(0,final_time)
f_case[1] = integrate(function(x) dgammaAlt(x,mean=5.1, cv=0.86), lower=0, upper=1.5)$value
for (i in 2:final_time){
  f_case[i] <- integrate(function(x) dgammaAlt(x,mean=5.1, cv=0.86), lower=i-0.5, upper=i+0.5)$value
}

stan_data <- list(M_regions= M_regions,
                  final_time=final_time,
                  W = week,
                  # gmobility = gmobility,
                  initial_seeding_day = initial_seeding_day,
                  death_data_length = death_len_data,
                  death= as.matrix(death_regions,death_len_data,M_regions),
                  SI=si,
                  f=f,
                  pop=pop,
                  C_base=C_base,
                  C_lockdown = C_lockdown,
                  day_week_index = day_week_index,
                  I = lockdown_index,
                  fitting_death_start = fitting_death_start
)

# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

# Example in R using rstan
m <- cmdstan_model("fitting_regions.stan")

fit <- m$sample(
  data=stan_data,
  iter_sampling = 10, 
  iter_warmup = 10, 
  parallel_chains = 4,
  chains=4,
  thin=1, 
  seed=12345,
  refresh = 20,
  adapt_delta = 0.8, 
  max_treedepth = 10)     # adapt_delta controls acceptance probability (lower -> larger step size, higher acceptance rate, less time, less explored posterior distribution
                                                              # opposite for tigher adapt_delta)
                                                              # default adapt_delta=0.80, max.treedepth = 10                                                              # max_treedepth is for time efficiency concern (higher -> more time)  
out <- fit$draws(format = "matrix")
summary_fit <- fit$summary()
# save(fit,stan_data,file=paste0("results/",'region_fitting_g_mob_cmd.Rdata'))

# bayesplot::mcmc_dens(fit$draws(c("mu[1]","mu[2]","mu[3]")))
# bayesplot::mcmc_scatter(fit$draws(c("mu[1]","mu[2]")))

weekly_death <- apply(fit$draws("weekly_deaths",format="matrix"),2,mean)
plot(weekly_death[((2*47)+1):(3*47)])
points(stan_data$death[,3],col="red")

Rt <- apply(fit$draws("Rt",format ="matrix"),2,mean)
plot(Rt[1:final_time])
points(Rt[(final_time+1):(2*final_time)],col="red")
points(Rt[((7*final_time)+1):(8*final_time)],col="red4")



# missing_case_regions <- data.frame(matrix(0,nrow=length(seq.Date(inf_start_date,as.Date(fitting_case_start)-1, by="day")),ncol = length(colnames(case_regions))))
# colnames(missing_case_regions) <- colnames(case_regions)
# missing_case_regions$date <- seq.Date(inf_start_date,as.Date(fitting_case_start)-1, by="day")
# case_regions <- rbind(missing_case_regions,case_regions)
# missing_case_october <- data.frame(matrix(0,nrow=length(seq.Date(as.Date("2020-10-01"),as.Date("2020-10-31"), by="day")),ncol = length(colnames(case_regions))))
# colnames(missing_case_october) <- colnames(case_regions)
# missing_case_october$date <- seq.Date(as.Date("2020-10-01"),as.Date("2020-10-31"), by="day")
# case_regions <- rbind(case_regions[case_regions$date <"2020-10-01",], missing_case_october,case_regions[case_regions$date>"2020-10-31",])
# case_regions$index <- c(rep(0,length(seq.Date(inf_start_date,as.Date(fitting_case_start)-1, by="day"))), 
#                        rep(1,length(seq.Date(as.Date(fitting_case_start),as.Date("2020-09-30"), by="day"))),
#                        rep(0,length(seq.Date(as.Date("2020-10-01"),as.Date("2020-10-31"), by="day"))),
#                        rep(1,length(seq.Date(as.Date("2020-11-01"),inf_end_date,by="day"))))
# 
# case_regions <- case_regions[!colnames(case_regions) == "date"]
# fitting_case_start <- which(rowSums(case_regions[,1:9]) != 0)[1]

# case_regions <- data.frame(north_east = apply(cases_2020[,north_east_index],1,sum),
#                            north_west = apply(cases_2020[,north_west_index],1,sum),
#                            yorkshire = apply(cases_2020[,yorkshire_index],1,sum),
#                            east_midlands = apply(cases_2020[,east_midlands_index],1,sum),
#                            west_midlands = apply(cases_2020[,west_midlands_index],1,sum),
#                            east = apply(cases_2020[,east_index],1,sum),
#                            london = apply(cases_2020[,london_index],1,sum),
#                            south_east = apply(cases_2020[,south_east_index],1,sum),
#                            south_west = apply(cases_2020[,south_west_index],1,sum))
# 
# case_regions$date <- case_date

# load("~/OneDrive - National University of Singapore/uk_mobility_data/data/cases/cases_2020.Rdata")       ## Daily data
# case_date <- cases_2020$date
# fitting_case_start <- min(case_date)
# cases_2020 <- cases_2020[,!names(cases_2020) == "date"]