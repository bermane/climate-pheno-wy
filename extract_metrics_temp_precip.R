#code to create various temp and precip metrics from daymet data
#the various indices will relate to mean pirgd over time period

#load packages
library(rgdal)
library(sp)
library(raster)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load mean pirgd values over 2002-2018
pirgd <- raster('./emodis/pirgd/pirgd_mean_wy_2002_2018.tif')

#Create climate indicies
#we want to create an index for each year
#and stack them together and export

#set variable to run and number of days prior to pirgd to average
vari <- "tmax"
window <- 30 #takes the mean of this number of days prior to pirgd

yr <- 2002

#load daymet data file list
files <- list.files(paste0('./', vari, '/mosaic'), ".tif$")
#remove 2001
files <- files[!files %in% grep("2001", files, value = T)]

#read in raster
ras <- raster(paste0('./', vari, '/mosaic/', vari, '_wy_' , yr, '.tif'))





