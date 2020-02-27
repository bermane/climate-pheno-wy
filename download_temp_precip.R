#load packages
library(daymetr)
library(sp)
library(rgdal)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load wy shapefile to get bounding box
wy <- readOGR('./reference/wyoming.shp')

#reproject to lat lon
wy <- spTransform(wy, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#change to vector and change order to match daymetr format needed
xy <- as.vector(extent(wy))
xy <- xy[c(4,1,3,2)]

download_daymet_ncss(xy, start = 2001, end = 2018, param = "tmax",
                     frequency = "daily", path = "./temp")
