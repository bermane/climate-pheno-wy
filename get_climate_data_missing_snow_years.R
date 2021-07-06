#this code downloads climate projection data for various scenarios

#load packages
library(ggplot2)
library(cft)
library(raster)
library(rgdal)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load pirgd data to transform into shapefile
pirgd <- raster('/Volumes/SSD/climate_effects/dlc/maxIRGdate_wy_laea_2000_2019.tif', band = 1)

#create spatial object
wy <- extent(pirgd)
wy <- as(wy, 'SpatialPolygons')
crs(wy) <- crs(pirgd)
rm(pirgd)

#list climate models we are interested in
mods <- c('MIROC-ESM-CHEM', 'IPSL-CM5A-MR', 'CNRM-CM5', 'inmcm4', 'HadGEM2-ES365')

#list parameters we are interested in
params <- c('tasmax', 'tasmin', 'pr', 'rhsmax', 'rhsmin')

#download the mid-century and end-of-century values separetely, will make it easier to average over time periods
#for middle scenario we also want data from 2020-2039 so that we can compute statistics over whole time period

#missing early
cftdata(aoi = wy, area_name = 'wy_projections', years = c(2018, 2020), local_dir = getwd(),
               models = 'HadGEM2-ES365',
               parameters = params, ncores = 8)

#missing mid
cftdata(aoi = wy, area_name = 'wy_projections', years = c(2038, 2040), local_dir = getwd(),
        models = mods,
        parameters = params, ncores = 8)

#missing end
cftdata(aoi = wy, area_name = 'wy_projections', years = c(2068, 2070), local_dir = getwd(),
        models = mods,
        parameters = params, ncores = 8)


#check spatial overlay
pirgd <- raster('/Volumes/SSD/climate_effects/dlc/maxIRGdate_wy_laea_2000_2019.tif', band = 1)
clim <- raster('/Volumes/SSD/climate_effects/wy_projections/tasmax_wy_projections_HadGEM2-ES365_r1i1p1_rcp45_macav2metdata_2040_2069_daily.tif', band = 1)
clim <- projectRaster(from = clim, to = pirgd)

plot(pirgd)
plot(clim, add = T)
