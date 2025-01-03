# the estimated infection at own region and at different regions due to the mobility
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
library(stringdist)
library(this.path)

script_directory <- this.path::this.dir()
setwd(script_directory)

previous_gmobility20 <- readRDS("mobility_change.rds") 
load("final_pop_2020_ltla.Rdata")

gmobility_2021 <- read.csv("2021_GB_Region_Mobility_Report.csv")
gmobility_2021 <- gmobility_2021 %>% select(!c("country_region_code","place_id","country_region","metro_area","iso_3166_2_code",
                             "census_fips_code","residential_percent_change_from_baseline")) %>% filter(date <= ymd(20210131)) %>% filter(sub_region_1 != "")

ltla_2020 <- sort(pop_2020$area_name)

Scotland_region <- c("Clackmannanshire","Glasgow City","Dumfries and Galloway","East Ayrshire Council","East Lothian Council","East Renfrewshire Council","Falkirk","Fife",
                     "Highland Council","Inverclyde","Midlothian","Moray","North Ayrshire Council","Perth and Kinross","Scottish Borders","South Ayrshire Council","South Lanarkshire",
                     "Stirling","Aberdeen City","Aberdeenshire","Argyll and Bute Council","Edinburgh","Renfrewshire","West Dunbartonshire Council","West Lothian",
                     "Angus Council","Dundee City Council","North Lanarkshire","East Dunbartonshire Council","Orkney","Shetland Islands","Eilean Siar","Na h-Eileanan an Iar")

Northern_Ireland_region <- c("Belfast","Ballymena","Magherafelt","Ards and North Down","Dungannon","Coleraine","Cookstown","Larne","Newry, Mourne and Down",
                             "Armagh City, Banbridge and Craigavon","Limvady","Derry and Strabane", "Lisburn and Castlereagh","Carrickfergus",
                             "Ballymoney","Fermanagh and Omagh","Moyle","Antrim and Newtownabbey","Causeway Coast and Glens","Mid and East Antrim","Mid Ulster")

Wales_region <- c("Isle of Anglesey","Gwynedd","Conwy Principal Area","Denbighshire","Flintshire","Wrexham","Ceredigion","Pembrokeshire","Carmarthenshire",
                  "Swansea","Neath Port Talbot Principle Area","Bridgend County Borough","Vale of Glamorgan","Cardiff","Rhondda Cynon Taff","Caerphilly County Borough","Blaenau Gwent",
                  "Torfaen Principal Area","Monmouthshire","Powys","Merthyr Tydfil","Newport","Merthyr Tydfil County Borough","Wrexham Principal Area")

gmobility_2021 <- gmobility_2021 %>% filter(!sub_region_1 %in% c(Scotland_region, Northern_Ireland_region, Wales_region))
gmobility_2021 <- gmobility_2021 %>% select(!sub_region_2)

colnames(gmobility_2021)[1] <- "area_name"
gmobility_2021 <- gmobility_2021 %>% mutate(area_name = ifelse(area_name == "Bristol, City of","Bristol City",  area_name))

previous_gmobility20 <- previous_gmobility20 %>% distinct(area_name, .keep_all = TRUE)

gmobility_2021 <- gmobility_2021 %>%
  left_join(select(previous_gmobility20, area_name, region), by = "area_name")

gmobility_2021$region[gmobility_2021$area_name %in% c("Cheshire East","Cheshire West and Chester","Cumbria")] = "North west"
gmobility_2021$region[gmobility_2021$area_name %in% c("County Durham","Darlington")] = "North east"
gmobility_2021$region[gmobility_2021$area_name %in% c("Derbyshire")] = "East midlands"
gmobility_2021$region[gmobility_2021$area_name %in% c("Devon")] = "South west"
gmobility_2021$region[gmobility_2021$area_name %in% c("East Riding of Yorkshire")] = "Yorkshire and the humber"
gmobility_2021$region[gmobility_2021$area_name %in% c("East Sussex")] = "South east"
gmobility_2021$region[gmobility_2021$area_name %in% c("West Midlands")] = "West midlands"

#### gmobility 2020 #####---------------------------------------------------------------------------------------------------

previous_gmobility20 <- readRDS("mobility_change.rds") 
load("final_pop_2020_ltla.Rdata")

gmobility_2020 <- read.csv("2020_GB_Region_Mobility_Report.csv")
gmobility_2020 <- gmobility_2020 %>% select(!c("country_region_code","place_id","country_region","metro_area","iso_3166_2_code",
                                               "census_fips_code","residential_percent_change_from_baseline"))


gmobility_2020 <- gmobility_2020 %>% filter(!sub_region_1 %in% c(Scotland_region, Northern_Ireland_region, Wales_region,""))
gmobility_2020 <- gmobility_2020 %>% select(!sub_region_2)

colnames(gmobility_2020)[1] <- "area_name"
gmobility_2020 <- gmobility_2020 %>% mutate(area_name = ifelse(area_name == "Bristol, City of","Bristol City",  area_name))

previous_gmobility20  <- previous_gmobility20  %>% distinct(area_name, .keep_all = TRUE)

gmobility_2020 <- gmobility_2020 %>%
  left_join(select(previous_gmobility20, area_name, region), by = "area_name")

gmobility_2020$region[gmobility_2020$area_name %in% c("Cheshire East","Cheshire West and Chester","Cumbria")] = "North west"
gmobility_2020$region[gmobility_2020$area_name %in% c("County Durham","Darlington")] = "North east"
gmobility_2020$region[gmobility_2020$area_name %in% c("Derbyshire")] = "East midlands"
gmobility_2020$region[gmobility_2020$area_name %in% c("Devon")] = "South west"
gmobility_2020$region[gmobility_2020$area_name %in% c("East Riding of Yorkshire")] = "Yorkshire and the humber"
gmobility_2020$region[gmobility_2020$area_name %in% c("East Sussex")] = "South east"
gmobility_2020$region[gmobility_2020$area_name %in% c("West Midlands")] = "West midlands"

# -----------check everything is fine -----------------------------------------------------------------------------------

sum(sort(unique(gmobility_2020$area_name))== sort(unique(gmobility_2020$area_name))) == length(unique(gmobility_2020$area_name))

gmobility <- rbind(gmobility_2020,gmobility_2021)

gmobility_region <- gmobility %>%
  group_by(region, date) %>%
  summarize(avg_retail_and_recreation_percent_change_from_baseline = mean(retail_and_recreation_percent_change_from_baseline, na.rm = TRUE),
            ave_grocery_and_pharmacy_percent_change_from_baseline = mean(grocery_and_pharmacy_percent_change_from_baseline, na.rm = TRUE),
            ave_transit_stations_percent_change_from_baseline = mean(transit_stations_percent_change_from_baseline, na.rm = TRUE),
            ave_parks_percent_change_from_baseline = mean(parks_percent_change_from_baseline, na.rm = TRUE),
            ave_workplaces_percent_change_from_baseline = mean(workplaces_percent_change_from_baseline, na.rm = TRUE),.groups = "drop")

saveRDS(gmobility_region, file = "gmobility_region.rds")




##### to check the similarity

# gmobility_2021_area <- unique(gmobility_2021$sub_region_1)
# a <- stringdistmatrix(unique(gmobility_2021$sub_region_1), Wales_region, method = "jaccard")
# b <- apply(a,2, function(x) min(x))




