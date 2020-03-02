#code to create various temp and precip metrics from daymet data
#the various indices will relate to mean pirgd over time period

#load packages
library(rgdal)
library(sp)
library(raster)

#set wd
setwd('/Volumes/SSD/climate_effects')