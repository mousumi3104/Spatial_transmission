stan_data_arrangements <- function(death_threshold, script_directory){
  
  setwd(script_directory) 
  load("data/final_pop_2020_ltla.Rdata") 
  load("data/england_death_2020.Rdata")       
  load("data/absolute_matrix_2011.Rdata")
  mobility_change <- readRDS("data/mobility_change.rds") 
  
  #------- regions index -----------------------------------------------------------------
  
  pop_2020$region <- sapply(pop_2020$region, function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  pop_2020$area_name <- sapply(pop_2020$area_name, function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  colnames(death_data) <- sapply(colnames(death_data), function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  death_north_east <- death_data %>% select(all_of(pop_2020$area_name[pop_2020$region == "North east"]))  
  
  #--------- start date of the epidemic ----------------------------------------------------------
  fitting_death_start <- which(cumsum(death_north_east$north_east) > death_threshold )[1]       # week when cumulative death exceeds above 10
  infection_gen_time <- 6 
  epidemic_start <- fitting_death_start - infection_gen_time
  inf_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", epidemic_start), "-1"))
  fitting_start_date <- ISOweek2date(paste0(2020, "-W", sprintf("%02d", fitting_death_start), "-1"))
  inf_end_date <- as.Date("31-12-2020",format = "%d-%m-%Y")
  # week_dates <- seq.Date(inf_start_date, inf_end_date, by = "week")
  final_time <- as.numeric(inf_end_date - inf_start_date) +1
  
  #-----------------------------------------------------------------------------------------------
  day_week_index <- array(0,final_time)
  for (t in 1:final_time){ 
    day_week_index[t] = ceiling(t/7) - (fitting_death_start - epidemic_start - 1)
  }
  
  day_week_index[1:((fitting_death_start-epidemic_start)*7)] = 1
  week <- day_week_index[length(day_week_index)]  
  
  #------- region_pop ----------------------------------------------------------------- 
  pop_north_east <- pop_2020 %>% filter(region == "North east")
  
  #--------- data arrangement ------------------------------------------------------------------
  
  M_regions <-nrow(pop_north_east)     # number of region
  mobility_matrix_2011 <- mobility_matrix_2011 %>% select(- place_of_work)
  
  colnames(mobility_matrix_2011) <- sapply(colnames(mobility_matrix_2011), function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  pop_2020$region <- sapply(pop_2020$region, function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  pop_2020$area_name <- sapply(pop_2020$area_name, function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  north_east_ltla <- pop_2020$area_name[pop_2020$region == "North east"]
  ltla_index <- which(colnames(mobility_matrix_2011) %in% north_east_ltla)
  C_base <- apply(as.matrix(mobility_matrix_2011[ltla_index, ltla_index]), 2, function(col) col/sum(col)) 
  
  C_lockdown <- matrix(0.0001,M_regions,M_regions)
  diag(C_lockdown) <- 0.9992
  
  initial_seeding_day = 6
  
  death_north_east <- death_north_east[epidemic_start: (epidemic_start + ceiling(final_time/7) -1 ),]   #(adjust accordingly)
  death_len_data <- nrow(death_north_east)
  fitting_death_start <- infection_gen_time + 1
  
  #--------------- google mobility -----------------------------------------------------------------------------
  
  # mobility_change <- readRDS("data/mobility_change.rds") 
  # mobility_change <- mobility_change %>% filter(region == "North east")
  # length(unique(mobility_change$area_name))
  # 
  gmobility <- read.csv("data/2020_GB_Region_Mobility_Report.csv")
  
  gmobility$sub_region_1 <- sapply(gmobility$sub_region_1, function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  gmobility$sub_region_2 <- sapply(gmobility$sub_region_2, function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  
  gmobility1 <- gmobility %>% filter(sub_region_1 %in% c("Hartlepool","Middlesbrough","Redcar and cleveland","Stockton-on-tees","Darlington","County durham","Northumberland"))
  gmobility2 <- gmobility %>% filter(sub_region_2 %in% c("Metropolitan borough of gateshead","Newcastle upon tyne district","North tyneside district","South tyneside","Sunderland district"))
  gmobility2 <- gmobility2 %>% mutate(sub_region_2 = case_when(sub_region_2 == "Metropolitan borough of gateshead" ~ "Gateshead",
                                                               sub_region_2 == "Newcastle upon tyne district" ~ "Newcastle upon tyne",
                                                               sub_region_2 == "North tyneside district" ~ "North tyneside",
                                                               sub_region_2 == "South tyneside" ~ "South tyneside",
                                                               sub_region_2 == "Sunderland district" ~ "Sunderland"))
  
  gmobility2$sub_region_1 <- gmobility2$sub_region_2
  
  gmobility <- rbind(gmobility1,gmobility2)
  gmobility <- gmobility %>% select(c("sub_region_1","date","retail_and_recreation_percent_change_from_baseline","grocery_and_pharmacy_percent_change_from_baseline",
                                      "transit_stations_percent_change_from_baseline","workplaces_percent_change_from_baseline","residential_percent_change_from_baseline"))
  
 NE_ltla <- colnames(north_east_matrix)
  
  gm <- array(NA, dim = c(length(unique(gmobility$date)), ncol(gmobility), M_regions))
 
   for (i in 1:M_regions){
    
    a <- gmobility %>% filter(sub_region_1 == NE_ltla[i])
    a$date <- as.Date(a$date)
    complete_dates <- seq(min(a$date), max(a$date), by = "day")
    missing_dates <- as.Date(setdiff(complete_dates, a$date))
    if (length(missing_dates) !=0){
      new_rows <- data.frame(
      sub_region_1 = NE_ltla[i],  # Or set a default region name if applicable
      date = missing_dates,
      retail_and_recreation_percent_change_from_baseline = NA,
      grocery_and_pharmacy_percent_change_from_baseline = NA,
      transit_stations_percent_change_from_baseline = NA,
      workplaces_percent_change_from_baseline = NA,
      residential_percent_change_from_baseline = NA)
      a <- rbind(a,new_rows)
      a <- a[order(a$date), ]
    }
    gm[,,i] <- as.matrix( a %>% filter(sub_region_1 == NE_ltla[i]))
  }
  
  for (i in 1:M_regions){
    a <- as.matrix((gm[,3:ncol(gmobility),i]))
    a <- apply(a, c(1, 2), as.numeric) 
    for (j in 1:ncol(a)){
      if (is.na(a[1,j])){a[1,j] = 1}
      for (k in 2:length(complete_dates)){
        if (is.na(a[k,j])){
          a[k,j] <- mean(c(as.numeric(a[max(which(!is.na(a[1:(k-1),j]))),j]),  as.numeric(a[min(which(!is.na(a[(k+1):length(complete_dates),j])))])))
        }
      }
    }
    gm[,3:ncol(gmobility),i] <- a
  }
  
  
  
  
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
  
  if (min(mobility_change_region$date) > inf_start_date){
    new_dates <- seq.Date(from = inf_start_date, to = min(mobility_change_region$date) -1, by = "day")
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
  
  # colors <- rainbow(M_regions)  # Generates M_regions distinct colors
  # matplot(gm, type = "l", lty = 1, col = colors, xlab = "Time", ylab = "Average Mobility")
  # legend("topright", legend = paste("Region", 1:M_regions), col = colors, lty = 1, cex = 0.8)
  
  #------------- lockdown effect --------------------------------------------------------------------------------
  
  lockdown1_started <- as.Date("2020-03-23", format = "%Y-%m-%d")
  lockdown1_lifted <- as.Date("2020-05-10", format = "%Y-%m-%d")
  
  lockdown2_started <- as.Date("2020-11-05", format = "%Y-%m-%d")  
  lockdown2_lifted <- as.Date("2020-12-02", format = "%Y-%m-%d") 
  
  lockdown_index <-  data.frame(date = seq.Date(from = inf_start_date,as.Date("2020-12-31",format = "%Y-%m-%d"),by="day"))
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
                    death= as.matrix(death_regions),
                    SI=dist$si,
                    f=dist$f,
                    pop=pop,
                    C_base=C_base,
                    C_lockdown = C_lockdown,
                    day_week_index = day_week_index,
                    I = as.matrix(lockdown_index),
                    first_lockdown_end = first_lockdown_end,
                    second_lockdown_end = second_lockdown_end,
                    fitting_death_start = fitting_death_start)
  
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
