#this code loads pre-processed data from the organize_covariates code and others
#and extracts the data on the sampling grid. Some processing is done to extract the correct
#variables and the result is a workspace image with a data table that can then be used for the analysis

#load packages
library(rgdal)
library(sp)
library(raster)
library(dynatopmodel)
library(tidyverse)
library(data.table)
library(lme4)
library(tictoc)
library(chillR)
library(scales)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load source functions
source("reference/HighstatLibV10.R")

###############
###LOAD DATA###
###############

#MAKE SURE TO ONLY LOAD 2001-2018!!!

#load sampling points
grid <- readOGR('reference/sampling_points_three_equal_classes.shp')

#load PIRGd
pirgd <- stack('dlc/maxIRGdate_wy_laea_2000_2019.tif')
pirgd <- pirgd[[2:19]]

#extract PIRGd and remove raster
dat <- raster::extract(pirgd, grid, df = T)
rm(pirgd)

#change col names
colnames(dat) <- c("id", 2001:2018)

#melt into df
dat <- reshape2::melt(dat, id.vars = "id")

#change ID values to be unique for each data point
dat$id <- rep(1:length(unique(dat$id)), times = 18)

#rename columns
dat <- dat %>% rename(c(year = variable, pirgd = value))
dat$year <- as.numeric(levels(dat$year))[dat$year]

#load elevation and twi
dem <- raster('dem/srtm/dem_wy_laea.tif')
twi <- raster('dem/srtm/twi_wy_laea.tif')

#extract and remove rasters
dat <- cbind(dat, dem = raster::extract(dem, grid) %>% rep(times = 18) %>% round(2),
             twi = raster::extract(twi, grid) %>% rep(times = 18) %>% round(2))
rm(dem, twi)

#load precip data
#first create precip columns then fill values in loop
dat$rain_oct_apr <- NA
dat$rain_jan_apr <- NA

#also need to calculate tmax thresholds to separate rain and snow using 
#linear regression of tmax vs. elevation best across all ecoregions
#from Rajagopal and Harpold (2016)
tmax_thresh <- 3.469 - (1.667*(dat$dem[dat$year == 2001]/1000))

#loop through years
for(yr in unique(dat$year)){
  
  #load previous yr and current year precip data
  prcp_yr <- brick(str_c('prcp/prcp_wy_laea_', yr - 1, '.tif'))
  prcp_yr1 <- brick(str_c('prcp/prcp_wy_laea_', yr, '.tif'))
  
  #load previous yr and current year max temp data to separate rain and snow
  tmax_yr <- brick(str_c('tmax/tmax_wy_laea_', yr - 1, '.tif'))
  tmax_yr1 <- brick(str_c('tmax/tmax_wy_laea_', yr, '.tif'))
  
  #extract points from oct to end of pirgd yr (will be faster then doing each separate)
  prcp <- raster::extract(prcp_yr[[274:365]], grid)
  prcp <- cbind(prcp, raster::extract(prcp_yr1[[1:227]], grid))
  tmax <- raster::extract(tmax_yr[[274:365]], grid)
  tmax <- cbind(tmax, raster::extract(tmax_yr1[[1:227]], grid))
  
  #rescale values divide by 100
  prcp <- prcp / 100
  tmax <- tmax / 100
  
  #calc total, snow, and rain prcp values from oct to apr
  prcp_oct_apr <- prcp[, 1:212]
  tmax_oct_apr <- tmax[, 1:212]
  rain_oct_apr <- prcp_oct_apr
  rain_oct_apr[tmax_oct_apr <= tmax_thresh] <- 0
  
  #calc rain from jan - apr
  rain_jan_apr <- prcp[, 93:212]
  tmax_jan_apr <- tmax[, 93:212]
  rain_jan_apr[tmax_jan_apr <= tmax_thresh] <- 0
  
  #input into dat
  dat$rain_oct_apr[dat$year == yr] <- rowSums(rain_oct_apr)
  dat$rain_jan_apr[dat$year == yr] <- rowSums(rain_jan_apr)
  
  #clean up
  rm(prcp_oct_apr, tmax_oct_apr, rain_oct_apr, rain_jan_apr)
  
}

#clean up
rm(prcp, prcp_yr, prcp_yr1, tmax, tmax_jan_apr, tmax_yr, tmax_yr1, tmax_thresh, yr)

#save dat as dat_hold so can load proper workspace
dat_hold <- dat

#load workspace so can add new column
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#cbind dat
dat <- cbind(dat, dat_hold$rain_oct_apr, dat_hold$rain_jan_apr)

#change col name
colnames(dat)[54:55] <- c('rain_oct_apr', 'rain_jan_apr')

#clean up
rm(dat_hold)

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#save back up
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_backup.RData')



