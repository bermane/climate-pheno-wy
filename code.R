#load relevant packages
library(rgdal)
library(sp)
library(raster)
library(dynatopmodel)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load wyoming shapefile
wy <- readOGR('./reference/wyoming.shp')
wy <- spTransform(wy, "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

###LOAD DATA###
#start with one year for 2002

#load emodis data
#load sosd and pgsd, crop and mask bad values
sosd <- raster('./emodis/sosd/SOST2002_wUSTeM250m_v1.bsq')
sosd <- crop(sosd, wy)
sosd[sosd == -1000 | sosd == 1000] <- NA
pgsd <- raster('./emodis/pgsd/MAXT2002_wUSTeM250m_v1.bsq')
pgsd <- crop(pgsd, wy)
pgsd[pgsd == -1000 | pgsd == 1000] <- NA

#calculate pirgd
pirgd <- sosd + (pgsd - sosd)/2

#remove sosd and pgsd
rm(sosd, pgsd)

#load snowmelt timing map
smt <- raster('./snowmelt/Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2002_V2.tif')
smt <- projectRaster(smt, pirgd, method = "bilinear")

#load elevation data
dem <- raster('./dem/dem_90m/w001001.adf')

#before doing aany transformations, calculate twi and solar
#calculate twi, export to disk and load again
#twi <- upslope.area(dem, log = TRUE, atb = TRUE, deg = 0.1, fill.sinks = TRUE)
#writeRaster(twi[[2]], "./dem/twi_wy.tif", format = "GTiff")
twi <- raster('./dem/twi_wy.tif')
twi <- projectRaster(twi, pirgd, method = "bilinear")

#load solar insolation. calculated in arcmap
insol <- raster('./solar/wy_solar_insolation.tif')
insol <- projectRaster(insol, pirgd, method = "bilinear")
