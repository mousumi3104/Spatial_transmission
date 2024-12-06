stan_data_arrangements <- function(death_threshold, script_directory){
  
  setwd(script_directory) 
  load("data/final_pop_2020_ltla.Rdata") 
  load("data/uk_regions_mobility_matrix.Rdata")
  death_data <- read_excel("data/daily_death_data_region.xlsx",sheet = 12)       
  
  #-------- death data arrangements -------------------------------------------------------
  colnames(death_data) <- death_data[3,]
  daily_death <- death_data[4:314,c(1,4,8:16)]
  region_order <- colnames(daily_death)[3:11]
  daily_death <- daily_death %>% mutate_all(~ifelse(is.na(.),0,.))
  daily_death$Date <- as.Date(daily_death$Date, format="%d/%m/%Y")
  daily_death$Date[278:311] <- seq(as.Date("28-11-2020",format = "%d-%m-%Y"),as.Date("31-12-2020",format = "%d-%m-%Y"),by="day")
  
  #--------- start date of the epidemic ----------------------------------------------------------
  fitting_death_start <- as.Date(daily_death$Date[which(cumsum(daily_death$England) > death_threshold )[1]],format = "%d/%m/%Y")       # day when cumulative death exceeds above 10
  infection_gen_time <- 45
  epidemic_start <- fitting_death_start - infection_gen_time
  epidemic_end <- as.Date("31-07-2020",format = "%d-%m-%Y")
  final_time <- as.numeric(epidemic_end - epidemic_start) +1
  
  daily_death <- daily_death[as.Date(daily_death$Date) >= as.Date(fitting_death_start),]
  daily_death <- daily_death[as.Date(daily_death$Date) <= as.Date(epidemic_end),]
  daily_death <- daily_death[!is.na(daily_death$Date), ]
  
  #------- region_pop ----------------------------------------------------------------- 
  
  region_order <- c("North east", "North west","Yorkshire and the Humber","East midlands","West midlands","East","London","South east","South west")
  population_region <- rep(0,length(region_order))
  for (i in 1:length(region_order)){
    population_region[i] = sum(pop_2020$population[pop_2020$region == region_order[i]])
  }
  
  #---------- check for missing dates within death data------------------------------------------
  # complete_date <- seq(epidemic_start, epidemic_end,by = "day")
  # as.Date(setdiff(complete_date, daily_death$Date),format = "%d-%m-%Y")
  
  #-----------------------------------------------------------------------------------------------
  day_week_index <- array(0,final_time)
  for (t in 1:final_time){ 
    day_week_index[t] = ceiling(t/7) - ceiling(infection_gen_time/7) + 1
  }
  
  day_week_index[1:infection_gen_time] = 1
  week <- day_week_index[length(day_week_index)]  
  
  #--------- data arrangement ------------------------------------------------------------------
  
  M_regions <- length(region_order)
  
  C_base <- mobmatrix_region_norm   # mobility matrix 
  C_lockdown <- matrix(0.0001,M_regions,M_regions)
  diag(C_lockdown) <- 0.9992
  
  # C_lockdown <- diag(M_regions)
  
  initial_seeding_day = 14
  # pop = pop_region$population[1:M_regions]          #pop_2020$population[1:M_regions]
  
  death_len_data <- nrow(daily_death)
  
  #--------------- google mobility -----------------------------------------------------------------------------
  
  mobility_change <- readRDS("data/mobility_change.rds") 
  
  # mobility_2020 <- read.csv("data/2020_GB_Region_Mobility_Report.csv")
  # gmob_area_name <- data.frame(area = unique(mobility_change$area_name), index = 1:length(unique(mobility_change$area_name)))
  # pop_area_name <- data.frame(area=pop_2020$area_name)
  # 
  # area_distance_matrix <- stringdistmatrix(gmob_area_name$area,pop_area_name$area, method = "jw")
  # matches <- expand.grid(Name1 = gmob_area_name$area, Name2 = pop_area_name$area)
  # matches$Distance <- as.vector(area_distance_matrix)
  # 
  # similar_matches <- matches[matches$Distance < 0.12, ]
  # similar_matches <- similar_matches %>% left_join(pop_2020, by = c("Name2" = "area_name"))
  # 
  # gmob_area_name$pop <- NA
  # gmob_area_name <- gmob_area_name %>% left_join(similar_matches, by = c("area" = "Name1")) %>%
  #                                     mutate(pop = population) %>% select(-population)
  
  mobility_change_area <- mobility_change %>%
    group_by(area_name, date, region) %>%
    summarize(avg_retail_and_recreation_percent_change_from_baseline = mean(retail_and_recreation_percent_change_from_baseline, na.rm = TRUE),
              ave_grocery_and_pharmacy_percent_change_from_baseline = mean(grocery_and_pharmacy_percent_change_from_baseline, na.rm = TRUE),
              ave_parks_percent_change_from_baseline = mean(parks_percent_change_from_baseline, na.rm = TRUE),
              ave_transit_stations_percent_change_from_baseline = mean(transit_stations_percent_change_from_baseline, na.rm = TRUE),
              ave_workplaces_percent_change_from_baseline = mean(workplaces_percent_change_from_baseline, na.rm = TRUE))
  
  mobility_change_region <- mobility_change_area %>%
    group_by(region, date) %>%
    summarize(avg_retail_and_recreation_percent_change_from_baseline = mean(avg_retail_and_recreation_percent_change_from_baseline, na.rm = TRUE),
              ave_grocery_and_pharmacy_percent_change_from_baseline = mean(ave_grocery_and_pharmacy_percent_change_from_baseline, na.rm = TRUE),
              ave_parks_percent_change_from_baseline = mean(ave_parks_percent_change_from_baseline, na.rm = TRUE),
              ave_transit_stations_percent_change_from_baseline = mean(ave_transit_stations_percent_change_from_baseline, na.rm = TRUE),
              ave_workplaces_percent_change_from_baseline = mean(ave_workplaces_percent_change_from_baseline, na.rm = TRUE),.groups = "drop")
  
  mobility_change_region$date <- as.Date(mobility_change_region$date, format = "%Y-%m-%d")
  
  if (min(mobility_change_region$date) > epidemic_start){
    new_dates <- seq.Date(from = epidemic_start, to = min(mobility_change_region$date) -1, by = "day")
    new_rows <- expand.grid(region = unique(mobility_change_region$region), date = new_dates)
    new_rows <- new_rows %>% mutate(avg_retail_and_recreation_percent_change_from_baseline = 0,
                                    ave_grocery_and_pharmacy_percent_change_from_baseline = 0,
                                    ave_parks_percent_change_from_baseline = 0,
                                    ave_transit_stations_percent_change_from_baseline = 0,
                                    ave_workplaces_percent_change_from_baseline = 0)
  }
  mobility_change_region <- rbind(new_rows, mobility_change_region)
  
  mobility_change_region$region <- factor(mobility_change_region$region, levels = region_order)
  mobility_change_region <- mobility_change_region %>% arrange(region)
  
  gmobility <- array(NA, dim = c(final_time,5,M_regions))
  for (i in 1:M_regions){
    gmobility[,,i] <- as.matrix(mobility_change_region[(((i-1)*(final_time))+1):(i*(final_time)),3:7])
  }
  
  for (m in 1:M_regions) {
    for (i in 1:5){
      for (t in 1:(final_time)){
        if (is.nan(gmobility[t,i,m])){
          prev_val <- if (t > 1) gmobility[max(which(!is.nan(gmobility[1:(t-1),i,m] ))),i,m]
          next_val <- if (t < final_time) gmobility[(min(which(!is.nan(gmobility[(t+1):final_time,i,m]))) + t),i,m]
          gmobility[t,i,m] <- (prev_val + next_val)/2
        }
      }
    }
  }
  
  gm <- matrix(NA, nrow = final_time, ncol = M_regions)
  for (m in 1:M_regions){
    gm[,m] = (gmobility[,1,m] + gmobility[,2,m] + gmobility[,3,m] + gmobility[,4,m] + gmobility[,5,m])/5;   
  }
  
  #------------- lockdown effect --------------------------------------------------------------------------------
  
  lockdown1_started <- as.Date("2020-03-23", format = "%Y-%m-%d")
  lockdown1_lifted <- as.Date("2020-05-10", format = "%Y-%m-%d")
  
  lockdown2_started <- as.Date("2020-11-05", format = "%Y-%m-%d")  
  lockdown2_lifted <- as.Date("2020-12-02", format = "%Y-%m-%d") 
  
  lockdown_index <-  data.frame(date = seq.Date(from = epidemic_start,as.Date("2020-12-31",format = "%Y-%m-%d"),by="day"))
  lockdown_index$r_mobility_index <- rep(0,nrow(lockdown_index))
  lockdown_index$r_mobility_index[lockdown_index$date >= lockdown1_started & lockdown_index$date < lockdown1_lifted] <- 1
  lockdown_index$r_mobility_index[lockdown_index$date >= lockdown2_started & lockdown_index$date < lockdown2_lifted] <- 1
  
  lockdown_index$g_mobility_index <-  rep(3,nrow(lockdown_index))
  lockdown_index$g_mobility_index[lockdown_index$date < lockdown1_lifted] <- 1
  lockdown_index$g_mobility_index[lockdown_index$date >= lockdown1_lifted & lockdown_index$date < lockdown2_lifted] <- 2
  
  first_lockdown_end <- which(lockdown_index$g_mobility_index == 2)[1]
  second_lockdown_end <- which(lockdown_index$g_mobility_index == 3)[1]
  
  lockdown_index$L1 <- ifelse(lockdown_index$g_mobility_index == 1,1,0)
  lockdown_index$L2 <- ifelse(lockdown_index$g_mobility_index == 2,1,0)
  lockdown_index$L3 <- ifelse(lockdown_index$g_mobility_index == 3,1,0)
  
  lockdown_index <- lockdown_index %>% select(- g_mobility_index, -date)
  lockdown_index <- lockdown_index[1:final_time,]
  
  
  
  source("data/distributions.R")
  dist <- distributions(final_time)
  
  stan_data <- list(M_regions= M_regions,
                    final_time=final_time,
                    W = week,
                    gmobility = gm,
                    initial_seeding_day = initial_seeding_day,
                    death_data_length = death_len_data,
                    death= as.matrix(daily_death[,3:11]),
                    SI=dist$si,
                    f=dist$f,
                    pop=population_region,
                    C_base=C_base,
                    C_lockdown = C_lockdown,
                    day_week_index = day_week_index,
                    I = as.matrix(lockdown_index),
                    first_lockdown_end = first_lockdown_end,
                    second_lockdown_end = second_lockdown_end,
                    infection_gen_time = infection_gen_time)
  
  return(stan_data)
  
}
# google mobility data arrangements ---------------------------------------------------------------------------

# data_google_mobility <- read.csv(here("data/2020_GB_Region_Mobility_Report.csv"))
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
