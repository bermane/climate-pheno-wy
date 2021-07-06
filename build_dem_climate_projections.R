#this code builds a dem in the projection and resolution of MACA climate projections
#so we can separate the observations by elevation categories

#load packages
library(rgdal)
library(sp)
library(raster)
library(dynatopmodel)
library(tidyverse)
library(data.table)
library(lme4)
library(tictoc)
library(chillR)
library(scales)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load sample climate layer
clim <- raster('wy_projections/raw/pr_wy_projections_CNRM-CM5_r1i1p1_rcp45_macav2metdata_2040_2069_daily.tif', band = 1)

#list dem files
dem_files <- list.files('dem/srtm/raw', full.names = T)

#mosaic rasters
dem_list <- lapply(1:length(dem_files),function(x) raster(dem_files[x]))
dem_list$fun <- mean

dem_mosaic <- do.call(mosaic, dem_list)

dem_ag <- dem_mosaic %>% aggregate(fact = 150)
dem_ag <- dem_ag %>% projectRaster(., clim) %>% crop(., clim)
writeRaster(dem_ag, filename = 'dem/srtm/dem_wy_maca_proj.tif', format = "GTiff")


