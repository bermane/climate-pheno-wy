#code to stack and crop emodis PIRGd calculated from SOSd and PGSd
#calculate mean pirgd in each pixel over time period
#and other summary statistics for all of WY and smaller SW sub-region
#writes files to disk to be able to use for future analyses

#load packages
library(rgdal)
library(sp)
library(raster)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load wyoming shapefile
wy <- readOGR('./reference/wyoming.shp')
#transform to CRS of emodis
wy <- spTransform(wy, "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

#load emodis data
#for some reason having problems changing values of raster stack... try one by one
#set empty pirgd stack and load sosd and pgsd file names
pirgd <- brick()
sosd_files <- list.files('./emodis/sosd/', '.bsq', full.names = T)
pgsd_files <- list.files('./emodis/pgsd/', '.bsq', full.names = T)

for(i in 1:length(sosd_files)){
#load sosd and pgsd raster
sosd <- raster(sosd_files[i])
pgsd <- raster(pgsd_files[i])

#crop them to wy
sosd<- crop(sosd, wy)
pgsd <- crop(pgsd, wy)

#mask bad values
sosd[sosd %in% c(-1000, 1000)] <- NA
pgsd[pgsd %in% c(-1000, 1000)] <- NA

#calculate pirgd and save to raster stack
pirgd_hold <- sosd + (pgsd - sosd)/2

if(i == 1){pirgd <- brick(pirgd_hold)} else{pirgd <- addLayer(pirgd, pirgd_hold)}

#remove raster files
rm(sosd, pgsd, pirgd_hold)
}

#change names of layers to years
names(pirgd) <- 2002:2018

#write pirgd stack to disk
#writeRaster(pirgd, './emodis/pirgd/pirgd_wy_2002_2018.tif', format = "GTiff")

#remove and reload brick for mean pirgd processing
rm(pirgd)
pirgd <- brick('./emodis/pirgd/pirgd_wy_2002_2018.tif')

#mask to WY before calc summary stats
pirgd <- mask(pirgd, wy)

#calc mean pirgd for each pixel
pirgd_m <- calc(pirgd, function(x) {mean(x,na.rm = T)})

#calc number of na values in annual layers
pirgd_na <- calc(pirgd, function(x) {sum(is.na(x))})

#only save pirgd mean value if less than half of years are missing
#we have 17 years total
pirgd_m[pirgd_na > 8] <- NA

#round to whole DOY
pirgd_m <- round(pirgd_m)

#write pirgd mean raster
#writeRaster(pirgd_m, './emodis/pirgd/pirgd_mean_wy_2002_2018.tif', format = "GTiff", overwrite = T)

#plot pirgd mean raster
plot(pirgd_m, main = "Mean PIRGd 2002-2018 in Wyoming")
