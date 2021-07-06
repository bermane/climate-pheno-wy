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

#load pirgd sd
pirgd_sd <- stack('dlc/sd_maxIRGdate_wy_laea_2000_2019.tif')
pirgd_sd <- pirgd_sd[[2:19]]

#extract spring scale and remove raster
dat <- cbind(dat, pirgd_sd = raster::extract(pirgd_sd, grid) %>% as.vector %>% round(2))
rm(pirgd_sd)

#load spring scale
ss <- stack('dlc/springScale_wy_laea_2000_2019.tif')
ss <- ss[[2:19]]

#extract spring scale and remove raster
dat <- cbind(dat, ss = raster::extract(ss, grid) %>% as.vector %>% round(2))
rm(ss)

#load landcover
lc <- raster('landcover/five_class_landcover_wy_laea_2016.tif')

#extract landcover and remove raster
dat <- cbind(dat, lc = raster::extract(lc, grid) %>% rep(times = 18))
rm(lc)

#load landcover change (1 = no change)
#lc_ch <- raster('landcover/landcover_change_wy_laea_2016.tif')

#extract landcover change and remove raster
#dat <- cbind(dat, lc_change = raster::extract(lc_ch, grid) %>% rep(times = 18))
#rm(lc_ch)

#load fire data
#fire <- raster('fire/burn_wy_laea_1984_2017.tif')

#extract fire and remove raster
#dat <- cbind(dat, fire = raster::extract(fire, grid) %>% rep(times = 18))
#rm(fire)

#change NA values to 9999 so we can subset below
#dat$fire[is.na(dat$fire) == 1] <- 9999

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

#plot elevation with vertical lines
hist(dat$dem, breaks = 30, main = "Histogram of elevation with equal quantile breaks at 1627 and 2281 m",
     xlab = "Elevation (m)")
abline(v = c(elev[1], elev[2]), col = c('blue', 'red'), lwd = 3)

#load snowmelt timing
snowmelt <- stack('snowmelt/snowmelt_timing_wy_laea_2001_2018.tif')

#extract and remove raster
dat <- cbind(dat, snowmelt = raster::extract(snowmelt, grid) %>% as.vector)

#zero values are NA
dat$snowmelt[dat$snowmelt == 0] <- NA
rm(snowmelt)

#load shrub component layers
herb_homer <- raster('shrublands/herbaceous/ann_herbaceous_wy_laea_2016.tif')
sage_homer <- raster('shrublands/sagebrush/sagebrush_wy_laea_2016.tif')
shrub_homer <- raster('shrublands/shrub/shrub_wy_laea_2016.tif')

#extract and remove rasters
dat <- cbind(dat, herb_homer = raster::extract(herb_homer, grid) %>% rep(times = 18),
             sage_homer = raster::extract(sage_homer, grid) %>% rep(times = 18),
             shrub_homer = raster::extract(shrub_homer, grid) %>% rep(times = 18))
rm(herb_homer, sage_homer, shrub_homer)

#change NA values to zero
dat$herb_homer[is.na(dat$herb_homer) == 1] <- 0
dat$sage_homer[is.na(dat$sage_homer) == 1] <- 0
dat$shrub_homer[is.na(dat$shrub_homer) == 1] <- 0

#load drought data, negative is drought
#we want to test average from previous year
#first create drought column then fill values in loop
dat$pdsi_apr_apr_mean <- NA
dat$pdsi_apr_apr_min <- NA
dat$pdsi_apr_apr_med <- NA
dat$pdsi_mar_apr_mean <- NA
dat$pdsi_mar_apr_min <- NA
dat$pdsi_mar_apr_med <- NA
dat$pdsi_jan_pirgd_mean <- NA
dat$pdsi_jan_pirgd_min <- NA
dat$pdsi_jan_pirgd_med <- NA
dat$pdsi_jan_apr_min <- NA

#loop through years
for(yr in unique(dat$year)){
  
  #load previous yr and current year drought data
  drought_yr <- stack(str_c('drought/pdsi_wy_laea_', yr - 1, '.tif'))
  drought_yr1 <- stack(str_c('drought/pdsi_wy_laea_', yr, '.tif'))
  
  #create stack of april to april data
  drought_apr_apr <- stack(drought_yr[[4:12]], drought_yr1[[1:3]])
  
  #calc mean across all months
  d_apr_apr_mean <- drought_apr_apr %>% calc(., fun = mean, na.rm = T)
  d_apr_apr_min <- drought_apr_apr %>% calc(., fun = min, na.rm = T)
  d_apr_apr_med <- drought_apr_apr %>% calc(., fun = median, na.rm = T)
  
  #extract points
  d_apr_apr_mean <- raster::extract(d_apr_apr_mean, grid) %>% round(2)
  d_apr_apr_min <- raster::extract(d_apr_apr_min, grid) %>% round(2)
  d_apr_apr_med <- raster::extract(d_apr_apr_med, grid) %>% round(2)
  
  #ensure there are the correct number of points -- data for each location
  #if there are add data
  if(length(d_apr_apr_mean) == length(unique(dat$id))) {
    dat$pdsi_apr_apr_mean[dat$year == yr] <- d_apr_apr_mean
    }else{
    print(str_c('Not enough data points for mean ', yr))} #if not don't
  
  if(length(d_apr_apr_min) == length(unique(dat$id))) {
    dat$pdsi_apr_apr_min[dat$year == yr] <- d_apr_apr_min
  }else{
    print(str_c('Not enough data points for max ', yr))} #if not don't
  
  if(length(d_apr_apr_med) == length(unique(dat$id))) {
    dat$pdsi_apr_apr_med[dat$year == yr] <- d_apr_apr_med
  }else{
    print(str_c('Not enough data points for mean ', yr))} #if not don't

  #clean up
  rm(drought_apr_apr, d_apr_apr_mean, d_apr_apr_min, d_apr_apr_med)
  
  #create stack of this year march april data
  drought_mar_apr <- stack(drought_yr1[[3:4]])
  
  #calc mean across all months
  d_mar_apr_mean <- drought_mar_apr %>% calc(., fun = mean, na.rm = T)
  d_mar_apr_min <- drought_mar_apr %>% calc(., fun = min, na.rm = T)
  d_mar_apr_med <- drought_mar_apr %>% calc(., fun = median, na.rm = T)
  
  #extract points
  d_mar_apr_mean <- raster::extract(d_mar_apr_mean, grid) %>% round(2)
  d_mar_apr_min <- raster::extract(d_mar_apr_min, grid) %>% round(2)
  d_mar_apr_med <- raster::extract(d_mar_apr_med, grid) %>% round(2)
  
  #ensure there are the correct number of points -- data for each location
  #if there are add data
  if(length(d_mar_apr_mean) == length(unique(dat$id))) {
    dat$pdsi_mar_apr_mean[dat$year == yr] <- d_mar_apr_mean
  }else{
    print(str_c('Not enough data points for mean ', yr))} #if not don't
  
  if(length(d_mar_apr_min) == length(unique(dat$id))) {
    dat$pdsi_mar_apr_min[dat$year == yr] <- d_mar_apr_min
  }else{
    print(str_c('Not enough data points for max ', yr))} #if not don't
  
  if(length(d_mar_apr_med) == length(unique(dat$id))) {
    dat$pdsi_mar_apr_med[dat$year == yr] <- d_mar_apr_med
  }else{
    print(str_c('Not enough data points for mean ', yr))} #if not don't
  
  #clean up
  rm(drought_mar_apr, d_mar_apr_mean, d_mar_apr_min, d_mar_apr_med)
  
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
  
  #load this years data, only through august
  drought_jan_pirgd <- stack(drought_yr1[[1:8]])
  
  #extract points
  drought_jan_pirgd <- raster::extract(drought_jan_pirgd, grid)
  
  #loop through each data point and calc
  for(i in 1:nrow(dat[dat$year == yr,])){
    if(is.na(dat$pirgd[dat$year == yr][i]) == 0){
      dat$pdsi_jan_pirgd_mean[dat$year == yr][i] <- mean(drought_jan_pirgd[i, 1:(month(as.Date(dat$pirgd[dat$year == yr][i] - 1, 
                                                                              format = "%j", origin = "1.1.2014")) - 1)],
                               na.rm = T) %>% round(2)
      dat$pdsi_jan_pirgd_min[dat$year == yr][i] <- min(drought_jan_pirgd[i, 1:(month(as.Date(dat$pirgd[dat$year == yr][i] - 1, 
                                                                            format = "%j", origin = "1.1.2014")) - 1)],
                              na.rm = T) %>% round(2)
      dat$pdsi_jan_pirgd_med[dat$year == yr][i] <- median(drought_jan_pirgd[i, 1:(month(as.Date(dat$pirgd[dat$year == yr][i] - 1, 
                                                                               format = "%j", origin = "1.1.2014")) - 1)],
                              na.rm = T) %>% round(2)
    }
  }
  
  #clean up
  rm(drought_jan_pirgd, drought_yr, drought_yr1)
  removeTmpFiles(h = 0.000000001)
}
rm(yr, i)

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#load solar data
solar <- stack('/Volumes/SSD/climate_effects/solar/solar_wy_laea_jan_aug_4_day.tif')

#extract points and remove raster
sol <- raster::extract(solar, grid)
rm(solar)

#pre allocate solar columns
dat$solar_jan_pirgd <- NA
#dat$solar_snowmelt_pirgd <- NA
dat$solar_jan_apr <- NA

#calc solar values
#values are in 4 day intervals
for(i in 1:nrow(dat)){
  
  #end of april is DOY 120 so solar_jan_apr is just sum until sol[,30]
  dat$solar_jan_apr[i] <- sum(sol[dat$id[i], 1:30])
  
  if(is.na(dat$pirgd[i]) == 0) {
  
    dat$solar_jan_pirgd[i] <- sum(sol[dat$id[i],1:(dat$pirgd[i] %/% 4)]) + 
      ((sol[dat$id[i],(dat$pirgd[i] %/% 4) + 1]) * ((dat$pirgd[i] %% 4) / 4)) %>% round(2)
    
    #only calc solar_snowmelt if snowmelt < pirgd
    #if(is.na(dat$snowmelt[i]) == F){
      #if(dat$snowmelt[i] <= dat$pirgd[i]){
        
        #dat$solar_snowmelt_pirgd[i] <- sum(sol[dat$id[i], ((dat$snowmelt[i] %/% 4) + 2):(dat$pirgd[i] %/% 4)]) +
          #((sol[dat$id[i],(dat$pirgd[i] %/% 4) + 1]) * ((dat$pirgd[i] %% 4) / 4)) #add on for pirgd end
        
        #if else add on for snowmelt start
        #if((dat$snowmelt[i] %% 4) == 0){
          #dat$solar_snowmelt_pirgd[i] <- dat$solar_snowmelt_pirgd[i] +
            #(sol[dat$id[i], (dat$snowmelt[i] %/% 4) + 1]) + ((sol[dat$id[i], (dat$snowmelt[i] %/% 4)]) * .25)
        #}else{
          #dat$solar_snowmelt_pirgd[i] <- dat$solar_snowmelt_pirgd[i] +
            #((sol[dat$id[i], (dat$snowmelt[i] %/% 4) + 1]) * ((5 - (dat$pirgd[i] %% 4)) / 4))
        #}
    #}
    
    #}
  }
}

#clean up
rm(sol)

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#load precip data
#first create precip columns then fill values in loop
dat$snow_oct_apr <- NA
dat$rain_oct_pirgd <- NA
dat$rain_oct_apr <- NA
#dat$rain_snowmelt_pirgd <- NA
dat$rain_jan_pirgd <- NA
dat$rain_jan_apr <- NA
dat$rain_low_mar_apr <- NA
dat$rain_med_mmar_mmay <- NA
dat$rain_hi_apr_may <- NA
dat$tot_prcp_oct_apr <- NA

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
  snow_oct_apr <- prcp_oct_apr
  snow_oct_apr[tmax_oct_apr > tmax_thresh] <- 0
  rain_oct_apr <- prcp_out_apr
  rain_oct_apr[tmax_oct_apr <= tmax_thresh] <- 0
  
  #calc rain from jan - apr
  rain_jan_apr <- prcp[, 93:212]
  tmax_jan_apr <- tmax[, 93:212]
  rain_jan_apr[tmax_jan_apr <= tmax_thresh] <- 0

  #input into dat
  dat$snow_oct_apr[dat$year == yr] <- rowSums(snow_oct_apr)
  dat$tot_prcp_oct_apr[dat$year == yr] <- rowSums(prcp_oct_apr)
  dat$rain_oct_apr[dat$year == yr] <- rowSums(rain_oct_apr)
  dat$rain_jan_apr[dat$year == yr] <- rowSums(rain_jan_apr)
  
  #clean up
  rm(prcp_oct_apr, tmax_oct_apr, snow_oct_apr, rain_oct_apr, rain_jan_apr)
  
  #calc rain at low elevation for mar-apr
  rain_low_mar_apr <- prcp[dat$dem[dat$year == yr] < elev[1], 152:212]
  tmax_low_mar_apr <- tmax[dat$dem[dat$year == yr] < elev[1], 152:212]
  rain_low_mar_apr[tmax_low_mar_apr <= tmax_thresh[dat$dem[dat$year == yr] < elev[1]]] <- 0
  
  #input into dat
  dat$rain_low_mar_apr[dat$year == yr & dat$dem < elev[1]] <- rowSums(rain_low_mar_apr)
  
  #clean up
  rm(rain_low_mar_apr, tmax_low_mar_apr)
  
  #calc rain at med elevation for mmar-mmay
  rain_med_mmar_mmay <- prcp[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2], 167:227]
  tmax_med_mmar_mmay <- tmax[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2], 167:227]
  rain_med_mmar_mmay[tmax_med_mmar_mmay <= tmax_thresh[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2]]] <- 0
  
  #input into dat
  dat$rain_med_mmar_mmay[dat$year == yr & dat$dem >= elev[1] & dat$dem < elev[2]] <- rowSums(rain_med_mmar_mmay)
  
  #clean up
  rm(rain_med_mmar_mmay, tmax_med_mmar_mmay)
  
  #calc rain at hi elevation for apr-may
  rain_hi_apr_may <- prcp[dat$dem[dat$year == yr] >= elev[2], 183:243]
  tmax_hi_apr_may <- tmax[dat$dem[dat$year == yr] >= elev[2], 183:243]
  rain_hi_apr_may[tmax_hi_apr_may <= tmax_thresh[dat$dem[dat$year == yr] >= elev[2]]] <- 0
  
  #input into dat
  dat$rain_hi_apr_may[dat$year == yr & dat$dem >= elev[2]] <- rowSums(rain_hi_apr_may)
  
  #clean up
  rm(rain_hi_apr_may, tmax_hi_apr_may)
  
  #iterate over each row in each to calc values that use pirgd
  for(i in 1:nrow(dat[dat$year == yr,])){
    if(is.na(dat$pirgd[dat$year == yr][i]) == 0) {
      
      #calc rain oct 1 - pirgd
      rain_oct_pirgd <- prcp[i, 1:(dat$pirgd[dat$year == yr][i] + 92)]
      tmax_oct_pirgd <- tmax[i, 1:(dat$pirgd[dat$year == yr][i] + 92)]
      rain_oct_pirgd[tmax_oct_pirgd <= tmax_thresh[i]] <- 0
      
      #input into dat
      dat$rain_oct_pirgd[dat$year == yr][i] <- sum(rain_oct_pirgd)
      
      #clean up
      rm(rain_oct_pirgd, tmax_oct_pirgd)
      
      #calc rain jan 1 - pirgd
      rain_jan_pirgd <- prcp[i, 93:(dat$pirgd[dat$year == yr][i] + 92)]
      tmax_jan_pirgd <- tmax[i, 93:(dat$pirgd[dat$year == yr][i] + 92)]
      rain_jan_pirgd[tmax_jan_pirgd <= tmax_thresh[i]] <- 0
      
      #input into dat
      dat$rain_jan_pirgd[dat$year == yr][i] <- sum(rain_jan_pirgd)
      
      #clean up
      rm(rain_jan_pirgd, tmax_jan_pirgd)
      
      #calc rain snowmelt - pirgd
      #only if snowmelt is before pirgd
      #if((is.na(dat$snowmelt[dat$year == yr][i]) == F)){
      #if(dat$snowmelt[dat$year == yr][i] <= dat$pirgd[dat$year == yr][i]){
        #rain_snowmelt_pirgd <- prcp[i, (dat$snowmelt[dat$year == yr][i] + 92):(dat$pirgd[dat$year == yr][i] + 92)]
        #tmax_snowmelt_pirgd <- tmax[i, (dat$snowmelt[dat$year == yr][i] + 92):(dat$pirgd[dat$year == yr][i] + 92)]
        #rain_snowmelt_pirgd[tmax_snowmelt_pirgd <= tmax_thresh[i]] <- 0
        
        #input into dat
        #dat$rain_snowmelt_pirgd[dat$year == yr][i] <- sum(rain_snowmelt_pirgd)
        
        #clean up
        #rm(rain_snowmelt_pirgd, tmax_snowmelt_pirgd)
      #}}
    }
  }
  
  #clean up
  rm(prcp, prcp_yr, prcp_yr1, tmax_yr, tmax_yr1, tmax)
}

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#load temp data for GDD
#first create GDD columns then fill values in loop
dat$gdd_jan_pirgd <- NA
#dat$gdd_snowmelt_pirgd <- NA
dat$gdd_jan_earlyavg_pirgd <- NA
dat$gdd_low_apr <- NA
dat$gdd_med_mapr_mmay <- NA
dat$gdd_hi_may <- NA

#loop through years
for(yr in unique(dat$year)){
  
  #load previous yr and current year min/max temp data
  #we will calculate the mean daily temp and use that as the threshold for GDD
  tmax_yr1 <- brick(str_c('tmax/tmax_wy_laea_', yr, '.tif'))
  tmin_yr1 <- brick(str_c('tmin/tmin_wy_laea_', yr, '.tif'))
  
  #extract points until last valid pirgd value
  tmax <- raster::extract(tmax_yr1[[1:227]], grid)
  tmin <- raster::extract(tmin_yr1[[1:227]], grid)
  
  #calc mean temperature and rescale
  tmean <- (tmax + tmin) / 2
  tmean <- tmean / 100
  
  #clean up
  rm(tmax_yr1, tmin_yr1, tmax, tmin)
  
  #set negative growing degrees to zero. use 5 degree threshold
  tmean <- tmean - 5
  tmean[tmean < 0] <- 0
  
  #calc gdd jan - DOY 95 (earliest average in state out of samples)
  gdd_jan_earlyavg_pirgd <- tmean[, 1:95]
  
  #input into dat
  dat$gdd_jan_earlyavg_pirgd[dat$year == yr] <- rowSums(gdd_jan_earlyavg_pirgd)
  
  #clean up
  rm(gdd_jan_earlyavg_pirgd)
  
  #calc gdd low elevation apr
  gdd_low_apr <- tmean[dat$dem[dat$year == yr] < elev[1], 91:120]
  
  #input into dat
  dat$gdd_low_apr[dat$year == yr & dat$dem < elev[1]] <- rowSums(gdd_low_apr)
  
  #clean up
  rm(gdd_low_apr)
  
  #calc gdd med elevation mapr - mmay
  gdd_med_mapr_mmay <- tmean[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2], 106:135]
  
  #input into dat
  dat$gdd_med_mapr_mmay[dat$year == yr & dat$dem >= elev[1] & dat$dem < elev[2]] <- rowSums(gdd_med_mapr_mmay)
  
  #clean up
  rm(gdd_med_mapr_mmay)
  
  #calc gdd hi elevation may
  gdd_hi_may <- tmean[dat$dem[dat$year == yr] >= elev[2], 121:151]
  
  #input into dat
  dat$gdd_hi_may[dat$year == yr & dat$dem >= elev[2]] <- rowSums(gdd_hi_may)
  
  #clean up
  rm(gdd_hi_may)
  
  #iterate over each row in each to calc values that use pirgd
  for(i in 1:nrow(dat[dat$year == yr,])){
    if(is.na(dat$pirgd[dat$year == yr][i]) == 0) {
      
      #calc gdd jan 1 - pirgd
      gdd_jan_pirgd <- tmean[i, 1:(dat$pirgd[dat$year == yr][i])]
      
      #input into dat
      dat$gdd_jan_pirgd[dat$year == yr][i] <- sum(gdd_jan_pirgd)
      
      #clean up
      rm(gdd_jan_pirgd)
      
      #calc gdd snowmelt - pirgd
      #only if snowmelt is before pirgd
      #if((is.na(dat$snowmelt[dat$year == yr][i]) == F)){
      #if(dat$snowmelt[dat$year == yr][i] <= dat$pirgd[dat$year == yr][i]){
        #gdd_snowmelt_pirgd <- tmean[i, (dat$snowmelt[dat$year == yr][i]):(dat$pirgd[dat$year == yr][i])]
        
        #input into dat
        #dat$gdd_snowmelt_pirgd[dat$year == yr][i] <- sum(gdd_snowmelt_pirgd)
        
        #clean up
        #rm(gdd_snowmelt_pirgd)
      #}}
    }
  }
  #clean up
  rm(tmean)
  
 }

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#load temp data for CHILLING
#first create CHILLING columns then fill values in loop
dat$chill_oct_apr <- NA
dat$chill_jan_apr <- NA
dat$chill_jan_may <- NA

#convert grid to lat/lon to extract lat later
grid_lat <- spTransform(grid, "+init=epsg:4326")

#only store latitude
grid_lat <- grid_lat@coords[,2]

#loop through years
for(yr in unique(dat$year)){
  
  #load previous yr and current year min/max temp data
  #we will calculate the mean daily temp and use that as the threshold for GDD
  tmax_yr <- brick(str_c('tmax/tmax_wy_laea_', yr - 1, '.tif'))
  tmin_yr <- brick(str_c('tmin/tmin_wy_laea_', yr - 1, '.tif'))
  tmax_yr1 <- brick(str_c('tmax/tmax_wy_laea_', yr, '.tif'))
  tmin_yr1 <- brick(str_c('tmin/tmin_wy_laea_', yr, '.tif'))
  
  #extract points for oct to may
  tmax <- raster::extract(tmax_yr[[274:365]], grid)
  tmax <- cbind(tmax, raster::extract(tmax_yr1[[1:151]], grid))
  tmin <- raster::extract(tmin_yr[[274:365]], grid)
  tmin <- cbind(tmin, raster::extract(tmin_yr1[[1:151]], grid))
  rm(tmax_yr, tmin_yr, tmax_yr1, tmin_yr1)
  
  #revalue temps
  tmax <- tmax / 100
  tmin <- tmin / 100
  
  #loop through and run code for each data point
  for(i in 1:nrow(dat[dat$year == yr,])){
    
    #create dataframe to output hourly temperatures
    df <- data.frame(Tmin = tmin[i,], Tmax = tmax[i,], Year = c(rep(yr - 1, times = 92), rep(yr, times = 151)),
                     JDay = c(274:365, 1:151))
    
    #extract hourly temps
    hr <- make_hourly_temps(latitude = grid_lat[i], year_file = df)
    hr <- hr[,-(1:4)]
    
    #convert to chilling hours, 1 if <= 7.2, 0 if > 7.2
    hr[hr <= 7.2] <- 1
    hr[hr > 7.2] <- 0
    
    #calc chilling hours for oct - apr
    hr_oct_apr <- rowSums(hr[1:212,]) %>% sum
    
    #calc chilling hours for jan - apr
    hr_jan_apr <- rowSums(hr[93:212,]) %>% sum
    
    #calc chilling hours for jan - may
    hr_jan_may <- rowSums(hr[93:243,]) %>% sum
    
    #into into dat
    dat$chill_oct_apr[dat$year == yr][i] <- hr_oct_apr
    dat$chill_jan_apr[dat$year == yr][i] <- hr_jan_apr
    dat$chill_jan_may[dat$year == yr][i] <- hr_jan_may
    
    #clean up
    rm(hr_oct_apr, hr_jan_apr, hr_jan_may, hr, df)
    
  }
  #clean up
  rm(tmax, tmin)
}

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#load vapor pressure data
#first create columns then fill values in loop
dat$vp_min_jan_apr <- NA
dat$vp_min_low_mar_apr <- NA
dat$vp_min_med_mmar_mmay <- NA
dat$vp_min_hi_apr_may <- NA
dat$vp_avg_jan_apr <- NA
dat$vp_avg_low_mar_apr <- NA
dat$vp_avg_med_mmar_mmay <- NA
dat$vp_avg_hi_apr_may <- NA

#loop through years
for(yr in unique(dat$year)){
  
  #load previous yr and current year vp data
  vp <- brick(str_c('vp/vp_wy_laea_', yr, '.tif'))
  
  #extract points from jan - may
  vp <- raster::extract(vp[[1:151]], grid)
  
  #calc total vp from jan - apr
  vp_jan_apr <- vp[, 1:120]
  
  #calc vp at low elevation for mar-apr
  vp_low_mar_apr <- vp[dat$dem[dat$year == yr] < elev[1], 60:120]
  
  #calc vp at med elevation for mmar-mmay
  vp_med_mmar_mmay <- vp[dat$dem[dat$year == yr] >= elev[1] & dat$dem[dat$year == yr] < elev[2], 75:135]
  
  #calc vp at hi elevation for apr-may
  vp_hi_apr_may <- vp[dat$dem[dat$year == yr] >= elev[2], 91:151]
  
  #input into dat
  dat$vp_min_jan_apr[dat$year == yr] <- apply(vp_jan_apr, 1, FUN = min)
  dat$vp_min_low_mar_apr[dat$year == yr & dat$dem < elev[1]] <- apply(vp_low_mar_apr, 1, FUN = min)
  dat$vp_min_med_mmar_mmay[dat$year == yr & dat$dem >= elev[1] & dat$dem < elev[2]] <- apply(vp_med_mmar_mmay, 1, FUN = min)
  dat$vp_min_hi_apr_may[dat$year == yr & dat$dem >= elev[2]] <- apply(vp_hi_apr_may, 1, FUN = min)
  dat$vp_avg_jan_apr[dat$year == yr] <- rowMeans(vp_jan_apr)
  dat$vp_avg_low_mar_apr[dat$year == yr & dat$dem < elev[1]] <- rowMeans(vp_low_mar_apr)
  dat$vp_avg_med_mmar_mmay[dat$year == yr & dat$dem >= elev[1] & dat$dem < elev[2]] <- rowMeans(vp_med_mmar_mmay)
  dat$vp_avg_hi_apr_may[dat$year == yr & dat$dem >= elev[2]] <- rowMeans(vp_hi_apr_may)
  
  #clean up
  rm(vp, vp_jan_apr, vp_low_mar_apr, vp_med_mmar_mmay, vp_hi_apr_may)
}

#save data after processing
#save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#load RAP data
#first create columns then fill values in loop
dat$ann_forb_rap <- NA #band 1
dat$bare_ground_rap <- NA #band 2
dat$perenn_forb_rap <- NA #band 4
dat$shrub_rap <- NA #band 5
dat$tree_rap <- NA #band 6
dat$ann_perenn_forb_rap <- NA #band 1+4

#loop through years
for(yr in unique(dat$year)){
  
  #load RAP layers
  rap <- brick(str_c('rap/rap_wy_laea_', yr, '.tif'))
  
  #input into dat
  dat$ann_forb_rap[dat$year == yr] <- raster::extract(rap[[1]], grid)
  dat$bare_ground_rap[dat$year == yr] <- raster::extract(rap[[2]], grid)
  dat$perenn_forb_rap[dat$year == yr] <- raster::extract(rap[[4]], grid)
  dat$shrub_rap[dat$year == yr] <- raster::extract(rap[[5]], grid)
  dat$tree_rap[dat$year == yr] <- raster::extract(rap[[6]], grid)
  dat$ann_perenn_forb_rap[dat$year == yr] <- (dat$ann_forb_rap[dat$year == yr] +  dat$perenn_forb_rap[dat$year == yr])
  
  #clean up
  rm(rap)
}

#change NA values to zero
dat$herb_homer[is.na(dat$herb_homer) == 1] <- 0
dat$sage_homer[is.na(dat$sage_homer) == 1] <- 0
dat$shrub_homer[is.na(dat$shrub_homer) == 1] <- 0

#clean up
rm(t)
rm(i)
rm(tmax_thresh)
rm(yr)
rm(elev)
rm(grid_lat)

#save data after processing
save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat.RData')

#save back up
save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_backup.RData')



