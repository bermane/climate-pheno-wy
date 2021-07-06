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

#we need to separate dem into three groups: low, med, high
#use quantile function to define the cutoffs
elev <- quantile(dat$dem, probs = c(1/3, 2/3))

#load gdd data
#first create columns then fill values in loop
dat$gdd_jan_apr <- NA
dat$gdd_mar_may <- NA
dat$gdd_apr_may <- NA

#loop through years
for(yr in unique(dat$year)){
  
  #load current year min/max temp data
  #we will calculate the mean daily temp and use that as the threshold for GDD
  tmax_yr1 <- brick(str_c('tmax/tmax_wy_laea_', yr, '.tif'))
  tmin_yr1 <- brick(str_c('tmin/tmin_wy_laea_', yr, '.tif'))
  
  #extract points until last valid pirgd value
  #only extract until the end of may
  tmax <- raster::extract(tmax_yr1[[1:151]], grid)
  tmin <- raster::extract(tmin_yr1[[1:151]], grid)
  
  #calc mean temperature and rescale
  tmean <- (tmax + tmin) / 2
  tmean <- tmean / 100
  
  #clean up
  rm(tmax_yr1, tmin_yr1, tmax, tmin)
  
  #set negative growing degrees to zero. use 5 degree threshold
  tmean <- tmean - 5
  tmean[tmean < 0] <- 0
  
  #calc gdd for different time periods
  gdd_jan_apr <- tmean[, 1:120]
  gdd_mar_may <- tmean[, 60:151]
  gdd_apr_may <- tmean[, 91:151]
  
  #input into dat
  dat$gdd_jan_apr[dat$year == yr] <- rowSums(gdd_jan_apr)
  dat$gdd_mar_may[dat$year == yr] <- rowSums(gdd_mar_may)
  dat$gdd_apr_may[dat$year == yr] <- rowSums(gdd_apr_may)
  
  #clean up
  rm(gdd_jan_apr, gdd_mar_may, gdd_apr_may, tmean)
}

#load rain data
#first create columns then fill values in loop
dat$rain_mar_may <- NA
dat$rain_apr_may <- NA

#also need to calculate tmax thresholds to separate rain and snow using 
#linear regression of tmax vs. elevation best across all ecoregions
#from Rajagopal and Harpold (2016)
tmax_thresh <- 3.469 - (1.667*(dat$dem[dat$year == 2001]/1000))

#loop through years
for(yr in unique(dat$year)){
  
  #load precip and temp data
  prcp <- brick(str_c('prcp/prcp_wy_laea_', yr, '.tif'))
  tmax <- brick(str_c('tmax/tmax_wy_laea_', yr, '.tif'))
  
  #extract points from mar to may
  prcp <- raster::extract(prcp[[60:151]], grid)
  tmax <- raster::extract(tmax[[60:151]], grid)
  
  #rescale values divide by 100
  prcp <- prcp / 100
  tmax <- tmax / 100
  
  #calc rain prcp values from mar to may
  rain_mar_may <- prcp
  rain_mar_may[tmax <= tmax_thresh] <- 0
  
  #calc rain from apr_may
  rain_apr_may <- rain_mar_may[,32:92]
  
  #input into dat
  dat$rain_mar_may[dat$year == yr] <- rowSums(rain_mar_may)
  dat$rain_apr_may[dat$year == yr] <- rowSums(rain_apr_may)
  
  #clean up
  rm(prcp, tmax, rain_mar_may, rain_apr_may)
}

#save dat as dat_hold so can load proper workspace
dat_hold <- dat

#temp save
save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_gdd_rain_hold_2020_10_09.RData')

#load workspace so can add new column
#load temp save
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_gdd_rain_hold_2020_10_09.RData')

#load main image
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_09_17.RData')

#bind new columns to data
dat$gdd_jan_apr <- dat_hold$gdd_jan_apr
dat$gdd_mar_may <- dat_hold$gdd_mar_may
dat$gdd_apr_may <- dat_hold$gdd_apr_may  
dat$rain_mar_may <- dat_hold$rain_mar_may
dat$rain_apr_may <- dat_hold$rain_apr_may

#clean up
rm(dat_hold, elev, yr)

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_10_18.RData')



