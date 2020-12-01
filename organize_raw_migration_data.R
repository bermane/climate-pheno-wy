#this code organizes raw deer migration data for varying sources
#to match the formatting and file type of Ellen's cleaned data
#(although this data is not cleaned in the same way)

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(rmapshaper)
library(foreign)

#set wd
setwd('/Volumes/SSD/climate_effects')

#############################################################
###ORGANIZE RAW MIGRATION DATA INTO SAME FORMAT AS ELLEN's###
#############################################################

#load ellen's data to see format
load('/Volumes/SSD/climate_effects/mig_data/Data from Ellen/Mule_WY_AtlanticRim_cleaned.RData')
dat_ref <- data.sub
rm(data.sub)

#load data table as df
dat_ref_df <- dat_ref@data

#load WY shapefile and change crs
wy <- readOGR('/Volumes/SSD/climate_effects/reference/wyoming.shp') %>% spTransform(., crs(dat_ref))

########################
###ATLANTIC RIM NORTH###
########################

#load file list
files <- list.files('/Volumes/SSD/climate_effects/mig_data/AR_North', full.names = T)
file_names <- list.files('/Volumes/SSD/climate_effects/mig_data/AR_North')

#create count
count <- 1

#loop through all individuals #UTM ZONE 13N
for(i in 1:length(files)){
  
  #load individual
  ind <- read.dbf(files[i])
  
  #create id column
  ind$AnimalID <- str_split(file_names[i], '_')[[1]][1] %>% str_replace(., 'gps', '') %>% as.integer
  
  #create timestamp
  ind$Timestamp <- str_c(ind$date, ' ', ind$hour) %>%
    as.POSIXct(., tz = 'GMT', format = '%m/%d/%Y %H')
  
  #create DOP
  if('status_tex' %in% colnames(ind)) {ind$DOP <- ind$status_tex}else{
    ind$DOP <- ind$fix_attemp
  }
  
  #create population label
  ind$Population <- 'AtlanticRimNorth'
  
  #create animal year id
  ind$AY_ID <- str_c(str_split(file_names[i], '_')[[1]][1] %>% str_replace(., 'gps', '') %>% as.integer,
                     '_', lubridate::year(ind$Timestamp))
  
  #create global id
  ind$GlobalID <- str_c('Mule_WY_', ind$Population, '_', ind$AY_ID)
  
  #check for problems with longitude

    #create spdf if i = 1 otherwise bind spdf
    if(count == 1){
      
      #create spdf
      dat_spdf <- SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                         proj4string = crs('+init=epsg:32613'),
                                         data = data.frame(AnimalID = ind$AnimalID,
                                                           Timestamp = ind$Timestamp,
                                                           DOP = ind$DOP,
                                                           Population = ind$Population,
                                                           AY_ID = ind$AY_ID,
                                                           GlobalID = ind$GlobalID))
      #add to count
      count <- count + 1
    }else{
      dat_spdf <- rbind(dat_spdf,
                        SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                               proj4string = crs('+init=epsg:32613'),
                                               data = data.frame(AnimalID = ind$AnimalID,
                                                                 Timestamp = ind$Timestamp,
                                                                 DOP = ind$DOP,
                                                                 Population = ind$Population,
                                                                 AY_ID = ind$AY_ID,
                                                                 GlobalID = ind$GlobalID)))  
      
    }
}

#change projection to lat lon
dat_spdf <- spTransform(dat_spdf, CRSobj = crs(dat_ref))

#check plot in WY
plot(wy)
plot(dat_spdf, add = T)

#check data frame
dat <- dat_spdf@data

#change name to match ellen's
data.sub <- dat_spdf

#save RData file with spdf
save(data.sub, file = "mig_data/Organized Data/Mule_WY_AtlanticRimNorth.RData")

#clean up
rm(dat, dat_spdf, data.sub, ind, i, count, file_names, files)

########################
###ATLANTIC RIM SOUTH###
########################

#load file list
files <- list.files('/Volumes/SSD/climate_effects/mig_data/AR_South', full.names = T)
file_names <- list.files('/Volumes/SSD/climate_effects/mig_data/AR_South')

#create count
count <- 1

#loop through all individuals #UTM ZONE 13N
for(i in 1:length(files)){
  
  #load individual
  ind <- read.dbf(files[i])
  
  #create id column
  ind$AnimalID <- str_split(file_names[i], '_')[[1]][1] %>% str_replace(., 'gps', '') %>% as.integer
  
  #create timestamp
  ind$Timestamp <- str_c(ind$date, ' ', ind$hour) %>%
    as.POSIXct(., tz = 'GMT', format = '%m/%d/%Y %H')
  
  #create DOP
  if('status_tex' %in% colnames(ind)) {
    ind$DOP <- ind$status_tex
    }else if('fix_attemp' %in% colnames(ind)){
    ind$DOP <- ind$fix_attemp
    }else ind$DOP <- ind$fix_succes
  
  #create population label
  ind$Population <- 'AtlanticRimSouth'
  
  #create animal year id
  ind$AY_ID <- str_c(str_split(file_names[i], '_')[[1]][1] %>% str_replace(., 'gps', '') %>% as.integer,
                     '_', lubridate::year(ind$Timestamp))
  
  #create global id
  ind$GlobalID <- str_c('Mule_WY_', ind$Population, '_', ind$AY_ID)
  
    #create spdf if i = 1 otherwise bind spdf
    if(count == 1){
      
      #create spdf
      dat_spdf <- SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                         proj4string = crs('+init=epsg:32613'),
                                         data = data.frame(AnimalID = ind$AnimalID,
                                                           Timestamp = ind$Timestamp,
                                                           DOP = ind$DOP,
                                                           Population = ind$Population,
                                                           AY_ID = ind$AY_ID,
                                                           GlobalID = ind$GlobalID))
      #add to count
      count <- count + 1
    }else{
      dat_spdf <- rbind(dat_spdf,
                        SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                               proj4string = crs('+init=epsg:32613'),
                                               data = data.frame(AnimalID = ind$AnimalID,
                                                                 Timestamp = ind$Timestamp,
                                                                 DOP = ind$DOP,
                                                                 Population = ind$Population,
                                                                 AY_ID = ind$AY_ID,
                                                                 GlobalID = ind$GlobalID)))  
      
    }
}

#change projection to lat lon
dat_spdf <- spTransform(dat_spdf, CRSobj = crs(dat_ref))

#check plot in WY
plot(wy)
plot(dat_spdf, add = T)

#check data frame
dat <- dat_spdf@data

#change name to match ellen's
data.sub <- dat_spdf

#save RData file with spdf
save(data.sub, file = "mig_data/Organized Data/Mule_WY_AtlanticRimSouth.RData")

#clean up
rm(dat, dat_spdf, data.sub, ind, i, count, file_names, files)

#########################
###Platte Valley North###
#########################

#load file list
files <- list.files('/Volumes/SSD/climate_effects/mig_data/PV_North', full.names = T)
file_names <- list.files('/Volumes/SSD/climate_effects/mig_data/PV_North')

#create count
count <- 1

#loop through all individuals
for(i in 1:length(files)){
  
  #load individual
  ind <- read.dbf(files[i])
  
  #create id column
  check <- str_split(file_names[i], '_')[[1]][1] 
  if(str_detect(check, 'd') == T){
    ind$AnimalID <- str_replace(check, 'd', '') %>% as.integer
  }else ind$AnimalID <- str_replace(check, 'gps', '') %>% as.integer

  #create timestamp
  ind$Timestamp <- str_c(ind$date, ' ', ind$hour) %>%
    as.POSIXct(., tz = 'GMT', format = '%m/%d/%Y %H')
  
  #create DOP
  if('status_tex' %in% colnames(ind)) {
    ind$DOP <- ind$status_tex
  }else if('fix_attemp' %in% colnames(ind)){
    ind$DOP <- ind$fix_attemp
  }else ind$DOP <- ind$fix_succes
  
  #create population label
  ind$Population <- 'PlatteValleyNorth'
  
  #create animal year id
  ind$AY_ID <- str_c(ind$AnimalID[1],
                     '_', lubridate::year(ind$Timestamp))
  
  #create global id
  ind$GlobalID <- str_c('Mule_WY_', ind$Population, '_', ind$AY_ID)
  
    #create spdf if i = 1 otherwise bind spdf
    if(count == 1){
      
      #create spdf
      dat_spdf <- SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                         proj4string = crs('+init=epsg:32613'),
                                         data = data.frame(AnimalID = ind$AnimalID,
                                                           Timestamp = ind$Timestamp,
                                                           DOP = ind$DOP,
                                                           Population = ind$Population,
                                                           AY_ID = ind$AY_ID,
                                                           GlobalID = ind$GlobalID))
      #add to count
      count <- count + 1
    }else{
      dat_spdf <- rbind(dat_spdf,
                        SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                               proj4string = crs('+init=epsg:32613'),
                                               data = data.frame(AnimalID = ind$AnimalID,
                                                                 Timestamp = ind$Timestamp,
                                                                 DOP = ind$DOP,
                                                                 Population = ind$Population,
                                                                 AY_ID = ind$AY_ID,
                                                                 GlobalID = ind$GlobalID)))  
      
    }
}

#change projection to lat lon
dat_spdf <- spTransform(dat_spdf, CRSobj = crs(dat_ref))

#check plot in WY
plot(wy)
plot(dat_spdf, add = T)

#check data frame
dat <- dat_spdf@data

#change name to match ellen's
data.sub <- dat_spdf

#make sure to record issues with longitude in external excel file!!

#save RData file with spdf
save(data.sub, file = "mig_data/Organized Data/Mule_WY_PlatteValleyNorth.RData")

#clean up
rm(dat, dat_spdf, data.sub, ind, i, count, file_names, files, check)

#########################
###Platte Valley South###
#########################

#load file list
files <- list.files('/Volumes/SSD/climate_effects/mig_data/PV_South', full.names = T)
file_names <- list.files('/Volumes/SSD/climate_effects/mig_data/PV_South')

#create count
count <- 1

#loop through all individuals
for(i in 1:length(files)){
  
  #load individual
  ind <- read.dbf(files[i])
  
  #create id column
  check <- str_split(file_names[i], '_')[[1]][1] 
  if(str_detect(check, 'd') == T){
    ind$AnimalID <- str_replace(check, 'd', '') %>% as.integer
  }else ind$AnimalID <- str_replace(check, 'gps', '') %>% as.integer
  
  #create timestamp
  ind$Timestamp <- str_c(ind$date, ' ', ind$hour) %>%
    as.POSIXct(., tz = 'GMT', format = '%m/%d/%Y %H')
  
  #create DOP
  if('status_tex' %in% colnames(ind)) {
    ind$DOP <- ind$status_tex
  }else if('fix_attemp' %in% colnames(ind)){
    ind$DOP <- ind$fix_attemp
  }else ind$DOP <- ind$fix_succes
  
  #create population label
  ind$Population <- 'PlatteValleySouth'
  
  #create animal year id
  ind$AY_ID <- str_c(ind$AnimalID[1],
                     '_', lubridate::year(ind$Timestamp))
  
  #create global id
  ind$GlobalID <- str_c('Mule_WY_', ind$Population, '_', ind$AY_ID)
  
    #create spdf if i = 1 otherwise bind spdf
    if(count == 1){
      
      #create spdf
      dat_spdf <- SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                         proj4string = crs('+init=epsg:32613'),
                                         data = data.frame(AnimalID = ind$AnimalID,
                                                           Timestamp = ind$Timestamp,
                                                           DOP = ind$DOP,
                                                           Population = ind$Population,
                                                           AY_ID = ind$AY_ID,
                                                           GlobalID = ind$GlobalID))
      #add to count
      count <- count + 1
    }else{
      dat_spdf <- rbind(dat_spdf,
                        SpatialPointsDataFrame(coords = matrix(c(ind$easting, ind$northing), ncol = 2),
                                               proj4string = crs('+init=epsg:32613'),
                                               data = data.frame(AnimalID = ind$AnimalID,
                                                                 Timestamp = ind$Timestamp,
                                                                 DOP = ind$DOP,
                                                                 Population = ind$Population,
                                                                 AY_ID = ind$AY_ID,
                                                                 GlobalID = ind$GlobalID)))  
      
    }
}

#change projection to lat lon
dat_spdf <- spTransform(dat_spdf, CRSobj = crs(dat_ref))

#check plot in WY
plot(wy)
plot(dat_spdf, add = T)

#check data frame
dat <- dat_spdf@data

#change name to match ellen's
data.sub <- dat_spdf

#save RData file with spdf
save(data.sub, file = "mig_data/Organized Data/Mule_WY_PlatteValleySouth.RData")

#clean up
rm(dat, dat_spdf, data.sub, ind, i, count, file_names, files, check)

########################
###SweetWaterGreenMtn###
########################

#load shape file of all data
ind <- readOGR('mig_data/SweetwaterGreenMtnData/33CleanDirectDownload/GreenMSweetWfilter3.shp')

#change date column to POSIXct
ind@data$date <- as.POSIXct(ind@data$date, tz = 'GMT', format = '%Y-%m-%d %H:%M:%OS')

#remove second date column
ind@data <- ind@data[, !(colnames(ind@data) %in% 'Date_1')]

#change column names
colnames(ind@data) <- c('AnimalID', 'Timestamp')

#add additional columns
ind@data$Population <- 'SweetWaterGreenMtn'
ind@data$AY_ID <- str_c(ind@data$AnimalID,
                        '_', lubridate::year(ind@data$Timestamp))
ind@data$GlobalID <- str_c('Mule_WY_', ind@data$Population, '_', ind@data$AY_ID)

#crs is already longlat
#change name to match ellen's
data.sub <- ind

#save RData file with spdf
save(data.sub, file = "mig_data/Organized Data/Mule_WY_SweetWaterGreenMtn.RData")

#clean up
rm(ind, data.sub)

##########
###DEER###
##########

#waiting to hear back from Katey about structure of data!!

ind <- readRDS('mig_data/DEERP_mig+premig_clean.rds')
#ind$AnimalID <- ind$AID
#ind$Timestamp <- as.POSIXct(ind, tz = 'GMT', format = '%m/%d/%Y %H')
#ind$DOP
#ind$Population
#ind$AY_ID
#ind$GlobalID


