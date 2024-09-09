stan_data_arrangements <- function(death_threshold, script_directory){

setwd(script_directory) 
load("data/final_pop_2020_ltla.Rdata")
load("data/england_death_2020.Rdata")       ## weekly data
load("data/uk_regions_mobility_matrix.Rdata")

#-------- death data rrangements -------------------------------------------------------
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

pop_2020$region <- sapply(pop_2020$region, function(x){
  paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})

death_data <- death_data %>% select(all_of(pop_2020$area_name))  
death_regions$total_death <- apply(death_data,1, sum)

#--------- start date of the epidemic ----------------------------------------------------------
fitting_death_start <- which(cumsum(death_regions$total_death) > death_threshold )[1]       # week when cumulative death exceeds above 10
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
   day_week_index[t] = ceiling(t/7) -(fitting_death_start - epidemic_start - 1)
}


day_week_index[1:((fitting_death_start-epidemic_start)*7)] = 1
week <- day_week_index[length(day_week_index)]  

#------- region_pop ----------------------------------------------------------------- 

pop_region <- data.frame(regions = unique(pop_2020$region))

population_region <- rep(0,nrow(pop_region))
for (i in 1:nrow(pop_region)){
  population_region[i] = sum(pop_2020$population[pop_2020$region == pop_region$regions[i]])
}

pop_region$population <- population_region

#--------- data arrangement ------------------------------------------------------------------

M_regions <-ncol(death_regions)-1     # number of region

C_base <- mobmatrix_region_norm   # mobility matrix 
C_lockdown <- diag(M_regions)

initial_seeding_day = 6
pop = pop_region$population[1:M_regions]          #pop_2020$population[1:M_regions]

death_regions <- death_regions[epidemic_start:nrow(death_regions),1:M_regions]   #(adjust accordingly)
death_len_data <- nrow(death_regions)
fitting_death_start <- infection_gen_time + 1

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
            ave_workplaces_percent_change_from_baseline = mean(ave_workplaces_percent_change_from_baseline, na.rm = TRUE),.groups = "drop")

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

gmobility <- array(NA, dim = c(final_time,5,M_regions))
for (i in 1:M_regions){
  gmobility[,,i] <- as.matrix(mobility_change_region[(((i-1)*(final_time))+1):(i*(final_time)),3:7])
}

for (m in 1:M_regions) {
  for (i in 1:4){
    for (t in 1:(final_time)){
      if (is.nan(gmobility[t,i,m])){
        
        prev_val <- if (t > 1) gmobility[max(which(!is.nan(gmobility[1:(t-1),i,m] ))),i,m]
        next_val <- if (t < final_time) gmobility[(min(which(!is.nan(gmobility[(t+1):final_time,i,m]))) + t),i,m]
        gmobility[t,i,m] <- (prev_val + next_val)/2
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

first_lockdown_end <- which(lockdown_index$g_mobility_index == 2)[1]
second_lockdown_end <- which(lockdown_index$g_mobility_index == 3)[1]

lockdown_index$L1 <- ifelse(lockdown_index$g_mobility_index == 1,1,0)
lockdown_index$L2 <- ifelse(lockdown_index$g_mobility_index == 2,1,0)
lockdown_index$L3 <- ifelse(lockdown_index$g_mobility_index == 3,1,0)

lockdown_index <- lockdown_index %>% select(- g_mobility_index, -date)



source("data/distributions.R")
dist <- distributions(final_time)

stan_data <- list(M_regions= M_regions,
                  final_time=final_time,
                  W = week,
                  gmobility = gmobility,
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
