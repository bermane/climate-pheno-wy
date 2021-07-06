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

#load annual homer shrubland components data
#first create columns then fill values in loop
dat$herb_homer_ann <- NA
dat$sage_homer_ann <- NA
dat$shrub_homer_ann <- NA

#load raster stacks
herb <- stack('homer_annual/herb_wy_laea_2001_2018_no_2012.tif')
sage <- stack('homer_annual/sage_wy_laea_2001_2018_no_2012.tif')
shrub <- stack('homer_annual/shrub_wy_laea_2001_2018_no_2012.tif')

#no data for 2012
year <- unique(dat$year) %>% .[!. %in% 2012]

#loop through years
for(i in 1:length(year)){
  
  #input into dat
  dat$herb_homer_ann[dat$year == year[i]] <- raster::extract(herb[[i]], grid)
  dat$sage_homer_ann[dat$year == year[i]] <- raster::extract(sage[[i]], grid)
  dat$shrub_homer_ann[dat$year == year[i]] <- raster::extract(shrub[[i]], grid)
  
}

#clean up
rm(herb, sage, shrub, i, year)

#change NA values to zero except for 2012 for which there is no data
dat$herb_homer_ann[is.na(dat$herb_homer_ann) == 1 & dat$year != 2012] <- 0
dat$sage_homer_ann[is.na(dat$sage_homer_ann) == 1 & dat$year != 2012] <- 0
dat$shrub_homer_ann[is.na(dat$shrub_homer_ann) == 1 & dat$year != 2012] <- 0

#save dat as dat_hold so can load proper workspace
dat_hold <- dat

#load workspace so can add new column
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#cbind dat
dat <- cbind(dat, dat_hold$herb_homer_ann, dat_hold$sage_homer_ann, dat_hold$shrub_homer_ann)

#change col name
colnames(dat)[56:58] <- c('herb_homer_ann', 'sage_homer_ann', 'shrub_homer_ann')

#clean up
rm(dat_hold)

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_homer_ann.RData')

#save back up
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_homer_ann_backup.RData')



