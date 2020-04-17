#load packages
library(daymetr)
library(sp)
library(rgdal)
library(tidyverse)

#set wd
setwd('/Volumes/SSD/climate_effects')

#tiles needed: queried from website
tiles <- c(12095, 12096, 12097, 12098,
           11915, 11916, 11917, 11918,
           11735, 11736, 11737, 11738)

#extra tiles to overlap with PIRGd map
tiles_ex <- c(12099, 11919, 11739)

for(tt in tiles_ex){
  for(yr in 2000:2019){
    for(param in c("tmin", "tmax", "prcp")) {
      download_daymet_tiles(tiles = tt, start = yr, end = yr, param = param, path = str_c("./", param))
    }
  }
}

#redownload corrupt file
download_daymet_tiles(tiles = , start = 2017, end = 2017, param = "tmax", path = "./temp")
