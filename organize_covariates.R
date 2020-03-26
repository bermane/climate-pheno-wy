#code to organize, crop, and output finalized covariates for analysis

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load 2016 land cover map
lc <- raster('/Volumes/SSD/climate_effects/landcover/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img')

#load dlc data to reproject land cover map
pirgd <- stack('./dlc/maxIRGdate.grd')

#reproject and crop land cover map
lc <- lc %>% projectRaster(lc, crs = crs(pirgd), method = "ngb")




#load snowmelt timing map
smt <- raster('./snowmelt/Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2002_V2.tif')
smt <- smt %>% projectRaster(pirgd, method = "bilinear") %>% crop(pirgd)
