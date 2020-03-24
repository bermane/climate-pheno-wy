#load packages
library(rgdal)
library(sp)
library(raster)
library(dynatopmodel)
library(tidyverse)

#set wd
setwd('/Volumes/SSD/climate_effects')

###LOAD DATA###
#start with one year for 2002

#load dlc data
pirgd <- stack('./dlc/maxIRGdate.grd')

#load WY shapefile and transform to proj of dlc
wy <- readOGR('./reference/wyoming.shp')
wy <- spTransform(wy, crs(pirgd))

#load snowmelt timing map
smt <- raster('./snowmelt/Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2002_V2.tif')
smt <- smt %>% projectRaster(pirgd, method = "bilinear") %>% crop(pirgd)

#load elevation data
dem <- raster('./dem/dem_90m/w001001.adf')
dem <- dem %>% projectRaster(pirgd, method = "bilinear")

#before doing aany transformations, calculate twi and solar
#calculate twi, export to disk and load again
#twi <- upslope.area(dem, log = TRUE, atb = TRUE, deg = 0.1, fill.sinks = TRUE)
#writeRaster(twi[[2]], "./dem/twi_wy.tif", format = "GTiff")
twi <- raster('./dem/twi_wy.tif')
twi <- twi %>% projectRaster(pirgd, method = "bilinear")

#load solar insolation. calculated in arcmap
insol <- raster('./solar/wy_solar_insolation.tif')
insol <- insol %>% projectRaster(pirgd, method = "bilinear")

