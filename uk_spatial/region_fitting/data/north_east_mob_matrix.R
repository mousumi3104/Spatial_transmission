#mobility matrix for London and Yorkshire and Hamber
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

  load("final_pop_2020_ltla.Rdata")
  load("uk_ltla_mobility_matrix.Rdata")       ## weekly data
  load("absolute_matrix_2011.Rdata")
  load("uk_regions_mobility_matrix.Rdata")
  load("england_death_2020.Rdata")       ## weekly data
  
  #-------- death data rrangements -------------------------------------------------------
  #------- regions index -----------------------------------------------------------------
  pop_2020$region <- sapply(pop_2020$region, function(x){
    paste0(toupper(substring(x,1,1)),tolower(substring(x,2)))})
  

  northeast_index <- which(pop_2020$region == "North east")

  #------- population ----------------------------------------------------------------- 
  
  pop_ne_ltla <- data.frame(ltla = pop_2020$area_name[pop_2020$region == "North east"] , pop = pop_2020$population[pop_2020$region == "North east"] )
  
 #-------- mobility matrix -------------------------------------------------------------

  mob_ne <- mobility_matrix_2011[1:12,2:13]
  mob_ne <- apply(mob_ne, 2, function(x) x/sum(x))
  
  #------ death data -------------------------------------------------------------------
  death_ne <- data.frame(north_east = apply(death_data[,north_east_index],1,sum))

  
  
  
  
  