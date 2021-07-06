#sample code to use NSD to identify a migration route for each individual and calculate distance
#between first and last point

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(rmapshaper)
library(foreign)
library(adehabitatLT)
library(sf)

#SET WORKING DIR
setwd('/Volumes/SSD/climate_effects/mig_data/Organized Data/')

#load WY shapefile
#wy <- readOGR('reference/wyoming.shp')

########################
###SELECT INDIVIDUALS###
########################

#FIRST LOAD WYOMING RANGE DATA
#load data
load('Mule_WY_WyomingRange_cleaned.RData')

#load as df
gps_df <- gps.data@data

#change date to POSIXct object
gps.data$Timestamp <- as.POSIXct(gps.data$Timestamp, tz = 'GMT')

#load individual
#need to look through trajectories manually to assess viability of migration
dat <- gps.data[gps.data$AnimalID == unique(gps.data$AnimalID)[1],]

#load as df
dat_df <- dat@data

#convert to ltraj object
dat <- as.ltraj(coordinates(dat), date = dat$Timestamp, id = dat$AnimalID, typeII=TRUE)

#plot NSD
plotltr(dat, which="R2n")

#find ideal number of segments
lav <- lavielle(dat, Lmin=10, Kmax=10, type="mean", which="R2n")
chooseseg(lav)

#try optimal number of segments
seg <- findpath(lav, 5)

#pick one spring migration segment for this individual
#check df to ensure it is spring!!
mig <- seg[[1]]

#find distance between first and last point
#create points object
pts <- st_sfc(st_point(as.numeric(mig[1, c('x','y')])),
              st_point(as.numeric(mig[NROW(mig), c('x','y')])),
              crs = as.character(crs(gps.data)))

#calculate distance between first and last point. > 15 KM? If so save individual and mig segment
st_distance(pts)/1000

#NEXT LOAD ATLANTIC RIM NORTH FROM HALL FOR ELLEN TO CHECK
#clear workspace
rm(list=ls())

#load data
load('Mule_WY_AtlanticRimNorth.RData')

#load individual
#need to look through trajectories manually to assess viability of migration
#play around with which animal gets loaded here
dat <- data.sub[data.sub$AnimalID == unique(data.sub$AnimalID)[5],]

#load as df
dat_df <- dat@data

#convert to ltraj object
dat <- as.ltraj(coordinates(dat), date = dat$Timestamp, id = dat$AnimalID, typeII=TRUE)

#plot NSD
plotltr(dat, which="R2n")

#can I use this data? Or should only use data where there is a clear start and end of migration? 
#i.e. collaring starts in winter range and goes at least unless animal is in summer range
