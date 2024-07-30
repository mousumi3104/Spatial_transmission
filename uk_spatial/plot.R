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
library(this.path)

script_directory <- this.path::this.dir()
setwd(script_directory)

load("result/natonal_fitting.Rdata")

