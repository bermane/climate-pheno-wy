#this code processes raw deer migration data for varying sources
#first sampling individuals from each herd and than
#using Visvalingam’s algorithm for line simplification and resampling of specific
#individuals at 4 km points

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

#load WY shapefile
wy <- readOGR('reference/wyoming.shp')

#########################
###LINE SIMPLIFICATION###
#########################

#load one of ellen's data files
load('/Volumes/SSD/climate_effects/mig_data/Data from Ellen/Mule_WY_AtlanticRim_cleaned.RData')

#create spdf of single individual
dat <- data.sub[data.sub$AnimalID == 1,]

#convert to spLines
dat_lines <- spLines(dat)

#simplify migration route
dat_s <- ms_simplify(dat_lines, keep = 0.05)
