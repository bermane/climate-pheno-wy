#code to mosaic together daymet tiles and crop to study area

#load packages
library(rgdal)
library(sp)
library(raster)

#load WY shapefile and transform to proj of daymet data
wy <- readOGR('/Volumes/SSD/climate_effects/reference/wyoming.shp')

###RUN FOR TEMP AND PRECIP AND ALL YEARS
#set initial variables
clim <- c("vp")
years <- 2000:2019

#run for all variables
for(vari in clim){
  
  #set vari wd
  setwd(paste0('/Volumes/SSD/climate_effects/',vari,'/raw'))
  
  #load all raw tile file names
  vari_f <- list.files()
  
  #check for any duplicate files
  if(anyDuplicated(vari_f) != 0){print(paste('Duplicate files!', vari, 'not processed!!!'))}
  else{
    print(paste('No duplicate files for', vari, '!!! Processing!!!'))
  
  #run for all years  
  for(yr in years){
  
#test for first year
#load annual files
vari_yr <- grep(yr, vari_f, value = T)

#check for all 12 tiles
if(length(vari_yr) != 15){print(paste('Not all tiles present!', vari, 'for', yr, 'not processed!!!'))}
else{
  print(paste('All tiles present!', vari, 'for', yr, 'processing!!!'))

#stack rasters in list
ras_list <- lapply(1:length(vari_yr), function(x){stack(vari_yr[x])})

#set parameters for mosaic function
names(ras_list)[1:2] <- c('x', 'y')
ras_list$fun <- mean
#ras_list$na.rm <- TRUE

#mosaic raster stacks
mos <- do.call(mosaic, ras_list)

#crop to WY
wy <- spTransform(wy, crs(mos))
mos <- crop(mos, wy)

#write to disk
writeRaster(mos, paste0('/Volumes/SSD/climate_effects/', vari ,'/mosaic/', vari, '_wy_', yr, '.tif'), format = "GTiff")

rm(mos, ras_list, vari_yr)
}
  }
  }
 
  rm(vari_f) 
}
