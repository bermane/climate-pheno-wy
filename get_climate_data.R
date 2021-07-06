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

#mid century
cftdata(aoi = wy, area_name = 'wy_projections', years = c(2040, 2069), local_dir = getwd(),
               models = mods,
               parameters = params, ncores = 2)

#end of century
cftdata(aoi = wy, area_name = 'wy_projections', years = c(2070, 2099), local_dir = getwd(),
        models = mods,
        parameters = params, ncores = 2)

#early century for middle scenario
cftdata(aoi = wy, area_name = 'wy_projections', years = c(2020, 2039), local_dir = getwd(),
        models = 'HadGEM2-ES365',
        parameters = params, ncores = 2)

#baseline for middle scenario. start 1999 because want snow data for year 2000.
cftdata(aoi = wy, area_name = 'wy_projections', years = c(1999, 2019), local_dir = getwd(),
        models = 'HadGEM2-ES365',
        parameters = params, ncores = 2)

#baseline for extreme scenario to compare. start 1999 because want snow data for year 2000.
cftdata(aoi = wy, area_name = 'wy_projections', years = c(1999, 2019), local_dir = getwd(),
        models = 'MIROC-ESM-CHEM',
        parameters = params, ncores = 2)

#baseline for rest of scenarios. start 1999 because want snow data for year 2000.
cftdata(aoi = wy, area_name = 'wy_projections', years = c(1999, 2019), local_dir = getwd(),
        models = c('IPSL-CM5A-MR', 'CNRM-CM5', 'inmcm4'),
        parameters = params, ncores = 2)

#check spatial overlay
pirgd <- raster('/Volumes/SSD/climate_effects/dlc/maxIRGdate_wy_laea_2000_2019.tif', band = 1)
clim <- raster('/Volumes/SSD/climate_effects/wy_projections/tasmax_wy_projections_HadGEM2-ES365_r1i1p1_rcp45_macav2metdata_2040_2069_daily.tif', band = 1)
clim <- projectRaster(from = clim, to = pirgd)

plot(pirgd)
plot(clim, add = T)
