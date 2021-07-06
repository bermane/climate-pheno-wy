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

#load drought data, negative is drought
#we want to test average from previous year
#first create drought column then fill values in loop
dat$pdsi_jan_apr_min <- NA

#loop through years
for(yr in unique(dat$year)){
  
  #load previous yr and current year drought data
  drought_yr1 <- stack(str_c('drought/pdsi_wy_laea_', yr, '.tif'))
  
  #create stack of this year jan april data
  drought_jan_apr <- stack(drought_yr1[[1:4]])
  
  #calc mean across all months
  d_jan_apr_min <- drought_jan_apr %>% calc(., fun = min, na.rm = T)
  
  #extract points
  d_jan_apr_min <- raster::extract(d_jan_apr_min, grid) %>% round(2)
  
  #ensure there are the correct number of points -- data for each location
  #if there are add data
  if(length(d_jan_apr_min) == length(unique(dat$id))) {
    dat$pdsi_jan_apr_min[dat$year == yr] <- d_jan_apr_min
  }else{
    print(str_c('Not enough data points for min ', yr))} #if not don't
  
  #clean up
  rm(drought_jan_apr, d_jan_apr_min)
  removeTmpFiles(h = 0.000000001)
}
rm(yr, drought_yr1)

#save dat as dat_hold so can load proper workspace
dat_hold <- dat

#load workspace so can add new column
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#cbind dat
dat <- cbind(dat, dat_hold$pdsi_jan_apr_min)

#change col name
colnames(dat)[53] <- 'pdsi_jan_apr_min'

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#save back up
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_backup.RData')



