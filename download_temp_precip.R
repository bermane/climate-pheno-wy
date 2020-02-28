#load packages
library(daymetr)
library(sp)
library(rgdal)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load wy shapefile to get bounding box
wy <- readOGR('./reference/wyoming.shp')

#reproject to lat lon
wy <- spTransform(wy, '+proj=longlat +datum=WGS84 +no_defs')

#change to vector and change order to match daymetr format needed
xy <- as.vector(extent(wy))
xy <- xy[c(4,1,3,2)]

#tiles needed: queried from website if bb doesn't work
tiles <- c(12095, 12096, 12097, 12098,
           11915, 11916, 11917, 11918,
           11735, 11736, 11737, 11738)

for(tt in tiles){
  for(yr in 2001:2018){
    for(param in c("tmax", "prcp")) {
      download_daymet_tiles(tiles = tt, start = yr, end = yr, param = param, path = "./temp")
    }
  }
}

#redownload corrupt file
download_daymet_tiles(tiles = 12095, start = 2017, end = 2017, param = "tmax", path = "./temp")