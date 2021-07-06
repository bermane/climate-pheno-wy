#this code reads in the migration data as .RData files and exports
#as shp files so can be read into QGIS and put on a map

library(raster)
library(sp)
library(rgdal)
library(tidyverse)

########################
###RAW ORGANIZED DATA###
########################

#read in file list w dir
files <- list.files('/Volumes/SSD/climate_effects/mig_data/Organized Data', pattern = '.RData', full.names = T)

#read in file names
names <- list.files('/Volumes/SSD/climate_effects/mig_data/Organized Data', pattern = '.RData', full.names = F)

#load reference for crs
load(files[1])
data_ex <- data.sub
rm(data.sub)

for(i in 1:length(files)){
  
  #read in RData file
  load(files[i])
  
  #make sure variable has correct name
  if(exists('gps.data')) {
    data.sub <- gps.data
    rm(gps.data)
  }
  
  #change to spdf if necessary
  if(class(data.sub)[1] != 'SpatialPointsDataFrame'){
    data.sub <- SpatialPointsDataFrame(coords = matrix(c(data.sub$lon, data.sub$lat), ncol = 2),
                                                   proj4string = crs(data_ex),
                                                   data = data.frame(AnimalID = data.sub$AnimalID,
                                                                     Timestamp = data.sub$Timestamp,
                                                                     Population = data.sub$Population,
                                                                     AY_ID = data.sub$AY_ID,
                                                                     GlobalID = data.sub$GlobalID))
  }
  
  #write as shp
  writeOGR(data.sub, 
           dsn = str_c('/Volumes/SSD/climate_effects/mig_data/Organized Data/shp/',
                                 str_replace(names[i], 'RData', 'shp')),
           layer = "migration",
           driver="ESRI Shapefile",
           overwrite_layer = T)
  
  #clean up
  rm(data.sub)
}

########################
###PROCESSED DATA###
########################

#read in file list w dir
files <- list.files('/Volumes/SSD/climate_effects/mig_data/processed', pattern = '.RData', full.names = T)

#read in file names
names <- list.files('/Volumes/SSD/climate_effects/mig_data/processed', pattern = '.RData', full.names = F)

for(i in 1:length(files)){
  
  #read in RData file
  load(files[i])
  
  if(exists('samples')){
    #write as shp
    writeOGR(samples, 
             dsn = str_c('/Volumes/SSD/climate_effects/mig_data/processed/shp/',
                         str_replace(names[i], 'RData', 'shp')),
             layer = "migration",
             driver="ESRI Shapefile",
             overwrite_layer = T)
    
    #clean up
    rm(samples)
  }
  
  if(exists('mig_start_end')){
    #write as shp
    writeOGR(mig_start_end, 
             dsn = str_c('/Volumes/SSD/climate_effects/mig_data/processed/shp/',
                         str_replace(names[i], 'RData', 'shp')),
             layer = "migration",
             driver="ESRI Shapefile",
             overwrite_layer = T)
    
    #clean up
    rm(mig_start_end)
  }
}

#save all start points for study area figure
#read in file list w dir
files <- list.files('/Volumes/SSD/climate_effects/mig_data/processed', pattern = 'seg_start', full.names = T)

#set count
count <- 1

#loop through files
for(i in files){
  #load data
  load(i)
  
  #remove end points (even rows)
  del <- seq(1, length(mig_start_end), 2)
  mig_hold <- mig_start_end[del,]
  crs(mig_hold) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
  
  if(count == 1){
    mig_start <- mig_hold
  }else{
    mig_start <- rbind(mig_start, mig_hold)
  }
  
  count <- count + 1
}

#write as shp
writeOGR(mig_start, 
         dsn = '/Volumes/SSD/climate_effects/mig_data/processed/shp/mig_start_all_herds.shp',
         layer = "migration",
         driver="ESRI Shapefile",
         overwrite_layer = T)

