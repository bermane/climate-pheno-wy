#code to calculate the percentage of LC classes of interest for methods section of paper

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(dynatopmodel)
library(foreach)
library(spatial.tools)
library(gdalUtils)
library(lubridate)
library(RCurl)

#set wd
setwd('/Volumes/SSD/climate_effects')


#reload landcover to set classes
lc <- raster('landcover/landcover_wy_laea_2016.tif')

#segment values we want to keep
lc_keep <- c(41, 42, 43, #deci, evergreen, mixed forests
             52, #shrub
             71, #herbaceous
             90, 95) #woody wetlands, herbaceous wetlands

#load values
lc_vals <- getValues(lc)

#set vals
shrub <- 52
herb <- 71
deci <- c(41, 43)
ever <- 42
wetl <- c(90, 95)

#perc shrub
length(lc_vals[lc_vals == shrub])/length(lc_vals)*100

#perc herb
length(lc_vals[lc_vals == herb])/length(lc_vals)*100

#perc ever
length(lc_vals[lc_vals == ever])/length(lc_vals)*100

#perc deci
length(lc_vals[lc_vals == deci])/length(lc_vals)*100

#perc wetl
length(lc_vals[lc_vals == wetl])/length(lc_vals)*100
