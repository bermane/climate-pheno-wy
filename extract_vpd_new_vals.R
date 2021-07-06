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

#load vapor pressure data
#first create columns then fill values in loop
dat$vpd_tmax_mean_jan_apr <- NA
dat$vpd_tmin_mean_jan_apr <- NA
dat$vpd_tavg_mean_jan_apr <- NA
dat$vpd_tmax_mean_elev <- NA
dat$vpd_tmin_mean_elev <- NA
dat$vpd_tavg_mean_elev <- NA

#loop through years
for(yr in unique(dat$year)){
  
  #load previous yr and current year vp and temp data
  vp <- brick(str_c('vp/vp_wy_laea_', yr, '.tif'))
  tmax <- brick(str_c('tmax/tmax_wy_laea_', yr, '.tif'))
  tmin <- brick(str_c('tmin/tmin_wy_laea_', yr, '.tif'))
  
  #extract points from jan - may
  vp <- raster::extract(vp[[1:151]], grid)
  tmax <- raster::extract(tmax[[1:151]], grid)
  tmin <- raster::extract(tmin[[1:151]], grid)
  
  #revalue temps
  tmax <- tmax / 100
  tmin <- tmin / 100
  
  #calc vpd variables
  vpd_tmax <- (610.7 * exp((17.38 * tmax)/(tmax + 239))) - vp
  vpd_tmin <- (610.7 * exp((17.38 * tmin)/(tmin + 239))) - vp
  vpd_tavg <- (((610.7 * exp((17.38 * tmax)/(tmax + 239))) + (610.7 * exp((17.38 * tmin)/(tmin + 239)))) / 2) - vp
  
  #calc vpd at low elevation for mar-apr
  vpd_tmax_low_mar_apr <- vpd_tmax[dat$dem[dat$year == yr] < elev[1], 60:120]
  vpd_tmin_low_mar_apr <- vpd_tmin[dat$dem[dat$year == yr] < elev[1], 60:120]
  vpd_tavg_low_mar_apr <- vpd_tavg[dat$dem[dat$year == yr] < elev[1], 60:120]
  
  #calc vp at med elevation for mmar-mmay
  vpd_tmax_med_mmar_mmay <- vpd_tmax[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2], 75:135]
  vpd_tmin_med_mmar_mmay <- vpd_tmin[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2], 75:135]
  vpd_tavg_med_mmar_mmay <- vpd_tavg[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2], 75:135]
  
  #calc vp at hi elevation for apr-may
  vpd_tmax_hi_apr_may <- vpd_tmax[dat$dem[dat$year == yr] >= elev[2], 91:151]
  vpd_tmin_hi_apr_may <- vpd_tmin[dat$dem[dat$year == yr] >= elev[2], 91:151]
  vpd_tavg_hi_apr_may <- vpd_tavg[dat$dem[dat$year == yr] >= elev[2], 91:151]
  
  #input into dat
  dat$vpd_tmax_mean_jan_apr[dat$year == yr] <- rowMeans(vpd_tmax[, 1:120])
  dat$vpd_tmin_mean_jan_apr[dat$year == yr] <- rowMeans(vpd_tmin[, 1:120])
  dat$vpd_tavg_mean_jan_apr[dat$year == yr] <- rowMeans(vpd_tavg[, 1:120])
  dat$vpd_tmax_mean_elev[dat$year == yr & dat$dem < elev[1]] <- rowMeans(vpd_tmax_low_mar_apr)
  dat$vpd_tmax_mean_elev[dat$year == yr & dat$dem >= elev[1] & dat$dem < elev[2]] <- rowMeans(vpd_tmax_med_mmar_mmay)
  dat$vpd_tmax_mean_elev[dat$year == yr & dat$dem >= elev[2]] <- rowMeans(vpd_tmax_hi_apr_may)
  dat$vpd_tmin_mean_elev[dat$year == yr & dat$dem < elev[1]] <- rowMeans(vpd_tmin_low_mar_apr)
  dat$vpd_tmin_mean_elev[dat$year == yr & dat$dem >= elev[1] & dat$dem < elev[2]] <- rowMeans(vpd_tmin_med_mmar_mmay)
  dat$vpd_tmin_mean_elev[dat$year == yr & dat$dem >= elev[2]] <- rowMeans(vpd_tmin_hi_apr_may)
  dat$vpd_tavg_mean_elev[dat$year == yr & dat$dem < elev[1]] <- rowMeans(vpd_tavg_low_mar_apr)
  dat$vpd_tavg_mean_elev[dat$year == yr & dat$dem >= elev[1] & dat$dem < elev[2]] <- rowMeans(vpd_tavg_med_mmar_mmay)
  dat$vpd_tavg_mean_elev[dat$year == yr & dat$dem >= elev[2]] <- rowMeans(vpd_tavg_hi_apr_may)
  
  #clean up
  rm(vp, tmax, tmin, vpd_tmax, vpd_tmin, vpd_tavg,
     vpd_tmax_low_mar_apr, vpd_tmin_low_mar_apr, vpd_tavg_low_mar_apr,
     vpd_tmax_med_mmar_mmay, vpd_tmin_med_mmar_mmay, vpd_tavg_med_mmar_mmay,
     vpd_tmax_hi_apr_may, vpd_tmin_hi_apr_may, vpd_tavg_hi_apr_may)
}

#save dat as dat_hold so can load proper workspace
dat_hold <- dat

#temp save
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_vpd_hold_2020_09_14.RData')

#load workspace so can add new column
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_09_04.RData')

#change vpd columns to match new values
dat$vpd_tmax_mean_jan_apr <- dat_hold$vpd_tmax_mean_jan_apr
dat$vpd_tmin_mean_jan_apr <- dat_hold$vpd_tmin_mean_jan_apr
dat$vpd_tavg_mean_jan_apr <- dat_hold$vpd_tavg_mean_jan_apr  
dat$vpd_tmax_mean_elev <- dat_hold$vpd_tmax_mean_elev
dat$vpd_tmin_mean_elev <- dat_hold$vpd_tmin_mean_elev
dat$vpd_tavg_mean_elev <- dat_hold$vpd_tavg_mean_elev

#change col name
#colnames(dat)[59:64] <- c('vpd_tmax_mean_jan_apr', 'vpd_tmin_mean_jan_apr', 'vpd_tavg_mean_jan_apr',
#                          'vpd_tmax_mean_elev', 'vpd_tmin_mean_elev', 'vpd_tavg_mean_elev')

#clean up
rm(dat_hold, elev, yr)

#save data after processing
save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_09_17.RData')



