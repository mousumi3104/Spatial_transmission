library(stringdist)
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
# 
script_directory <- this.path::this.dir()
setwd(script_directory)

mobility_2020 <- read.csv("2020_GB_Region_Mobility_Report.csv") #google mobility
unitary_counties <- unique(mobility_2020$sub_region_1)  #counties and unitary authorities (google mobility)

counties_unitary <- read.csv("~/downloads/Counties_and_Unitary_Authorities_December_2021_EN_BFC_2022_2329911750846456281.csv")   # check with google mobility data is from ONS website

wales_region <- c("Isle of Anglesey","Gwynedd","Conwy","Denbighshire","Flintshire","Wrexham","Ceredigion","Pembrokeshire","Carmarthenshire",
                  "Swansea","Neath Port Talbot","Bridgend","The Vale of Glamorgan","Cardiff","Rhondda Cynon Taf","Caerphilly","Blaenau Gwent",
                  "Torfaen","Monmouthshire","Powys","Merthyr Tydfil","Newport")



scotland_region <- c("Clackmannanshire","Glasgow City","Dumfries and Galloway","East Ayrshire","East Lothian","East Renfrewshire","Falkirk","Fife",
                     "Highland","Inverclyde","Midlothian","Moray","North Ayrshire","Perth and Kinross","Scottish Borders","South Ayrshire","South Lanarkshire",
                     "Stirling","Aberdeen City","Aberdeenshire","Argyll and Bute","City of Edinburgh","Renfrewshire","West Dunbartonshire","West Lothian",
                     "Angus","Dundee City","North Lanarkshire","East Dunbartonshire","Orkney Islands","Shetland Islands","Eilean Siar")



northern_ireland_region <- c("Belfast","Ballymena","Magherafelt","North Down","Dungannon","Coleraine","Cookstown","Larne","Strabane","Newry and Mourne",
                             "Omagh","Armagh","Castlereagh","Down","Antrim","Limvady","Derry","Craigavon","Banbridge","Ards","Lisburn","Carrickfergus",
                             "Ballymoney","Fermanagh","Moyle","Newtownabbey")

##find similar word
distance_matrix_scotland <- stringdistmatrix(scotland_region,unitary_counties, method = "jw")
matches_scotland <- expand.grid(Name1 = scotland_region, Name2 = unitary_counties)
matches_scotland$Distance <- as.vector(distance_matrix_scotland)

distance_matrix_wales <- stringdistmatrix(wales_region,unitary_counties, method = "jw")
matches_wales <- expand.grid(Name1 = wales_region, Name2 = unitary_counties)
matches_wales$Distance <- as.vector(distance_matrix_wales)

distance_matrix_ireland <- stringdistmatrix(northern_ireland_region,unitary_counties, method = "jw")
matches_ireland <- expand.grid(Name1 = northern_ireland_region, Name2 = unitary_counties)
matches_ireland$Distance <- as.vector(distance_matrix_ireland)

# Filter matches with a low distance (threshold)
threshold_scotland <- 0.13  # Adjust based on similarity needs
threshold_wales <- 0.142
threshold_ireland <- 0.13 
  
similar_matches_scotland <- matches_scotland[matches_scotland$Distance < threshold_scotland, ]
similar_matches_wales <- matches_wales[matches_wales$Distance < threshold_wales, ]
similar_matches_ireland <- matches_ireland[matches_ireland$Distance < threshold_ireland, ]


mobility_2020$country = NA
mobility_2020$country[mobility_2020$sub_region_1 %in% similar_matches_scotland$Name1] <- "Scotland"
mobility_2020$country[mobility_2020$sub_region_1 %in% similar_matches_wales$Name1] <- "Wales"
mobility_2020$country[mobility_2020$sub_region_1 %in% similar_matches_ireland$Name1] <- "Northern Ireland"


########-----------  deleting the rows with subregions in scotland, wales and noryher ireland-------------------------------------------------------------------------#######

mobility_2020 <- mobility_2020 %>% filter(!country %in% c("Scotland","Wales","Northern Ireland"))
gmob_subregion1 <- unique(mobility_2020$sub_region_1)
gmob_subregion2 <- unique(mobility_2020$sub_region_2)
gmob_subregion <- union(gmob_subregion1,gmob_subregion2)

#---- population data ---------------------------------------------------------------------------------

load("final_pop_2020_ltla.Rdata") 
pop_region <- pop_2020$area_name

#------------------------------------------------------------------------------------------------------

gmob_subregion_extra <- setdiff(gmob_subregion, pop_region)
pop_region_extra <- setdiff(pop_region, gmob_subregion)
gmob_subregion2_extra <-setdiff(gmob_subregion2,pop_region)






