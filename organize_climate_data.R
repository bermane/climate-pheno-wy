#this code organizes and stores relevant data metrics to be used for future scenario modelling
#it calculates our model variables of interest for each year in the climate time series
#so we can than take the mean values to define the time period

#load packages
library(raster)
library(rgdal)
library(tidyverse)
library(lubridate)
library(SPEI)
library(doParallel)

#set wd
setwd('/Volumes/SSD/climate_effects')

#set raster options
rasterOptions()
rasterOptions(maxmemory = 1e+10)
#rasterOptions()

#climate models we downloaded
mods <- c('MIROC-ESM-CHEM', 'IPSL-CM5A-MR', 'CNRM-CM5', 'inmcm4', 'HadGEM2-ES365')

#parameters we downloaded
params <- c('tasmax', 'tasmin', 'pr')

#year periods we downloaded (start year only). There is a third year period for middle scenario
years <- c(2040, 2070)

#scenario periods we downloaded
rcp <- c('rcp45', 'rcp85')

#we need a dem in the climate data project to define the elevation categories
#the thresholds from model building are: 1627 and 2281 m


################################
###CALCULATE MONTHLY AVERAGES###
################################    

#register parallel backend
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

#loop through models (possibly in parallel?)
foreach::foreach(m = mods) %dopar% {
  
  #load packages
  library(raster)
  library(rgdal)
  library(tidyverse)
  library(lubridate)
  library(SPEI)
  
  #loop through time periods
  for(yr in years){
    
    #loop through rcp scenarios
    for(r in rcp){
      
      #first we need mean monthly temps
      #load daily tmin and tmax
      tmin <- brick(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')))
      tmax <- brick(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*', yr, '*.tif')))
      
      #get a stack of daily average temps
      tavg <- overlay(tmin, tmax, fun = function(x,y){(x+y)/2},
                      filename = str_c('wy_projections/temp/tavg_', m, '_', r, '_', yr, '.tif'),
                      format = 'GTiff')
      
      #convert to monthly average temps
      #create a time series of dates
      #strip start and end year
      start_yr <- str_extract(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')),
                              "[:digit:]{4}_[:digit:]{4}") %>% str_sub(start = 1, end = 4)
      end_yr <- str_extract(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')),
                            "[:digit:]{4}_[:digit:]{4}") %>% str_sub(start = 6, end = 9)
      
      #create time series of years and months
      ts <- seq(mdy(str_c('01-01-', start_yr)), mdy(str_c('12-31-', end_yr)), by = "day") %>%
        format(., format = '%Y%m') %>% as.numeric
      
      #create index to match number of unique values
      index <- 1:length(unique(ts))
      names(index) <- unique(ts)
      
      #record ts values with unique months
      ts <- dplyr::recode(ts, !!!index)
      
      #calculate monthly average temps
      tavg_m <- stackApply(tavg, indices = ts, fun = mean, na.rm = T,
                           filename = str_c('wy_projections/temp/tavg_monthly_', m, '_', r, '_', yr, '.tif'),
                           format = 'GTiff')
      
      #load precipitation brick
      pr <- brick(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*', yr, '*.tif')))
      
      #compute monthly totals
      pr_m <- stackApply(pr, indices = ts, fun = sum, na.rm = T,
                         filename = str_c('wy_projections/temp/pr_monthly_', m, '_', r, '_', yr, '.tif'),
                         format = 'GTiff')
      
      
    }
  }
}

#stop parallel cluster
parallel::stopCluster(cl)

#######################################################################
###still need to run first part with middle model data starting 2020###
#######################################################################

######################################
###CALCULATE ANNUAL MODEL VARIABLES###
######################################  

#just for PIRGd model. will do spring scale later
#current final model variables:
#gdd_jan_apr
#rain_mar_may
#snow_oct_apr
#pdsi_mar_apr_min

#first run rcp45 2040 complete took ~28 hours

#register parallel backend
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

#loop through models
foreach::foreach(m = mods) %dopar% {
  
  #load packages
  library(raster)
  library(rgdal)
  library(tidyverse)
  library(lubridate)
  library(SPEI)
  library(scPDSI)
  
  #loop through time periods
  #for(yr in years){
  for(yr in 2070){
    
    #loop through rcp scenarios
    #for(r in rcp){
    for(r in "rcp85"){
      
      #strip start and end year
      start_yr <- str_extract(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')),
                              "[:digit:]{4}_[:digit:]{4}") %>% str_sub(start = 1, end = 4)
      end_yr <- str_extract(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')),
                            "[:digit:]{4}_[:digit:]{4}") %>% str_sub(start = 6, end = 9)
      
      #calc number of years
      num_years <- as.numeric(end_yr) - as.numeric(start_yr) + 1
      
      #load monthly average temperature time series
      tavg_m <- brick(str_c('wy_projections/temp/tavg_monthly_', m, '_', r, '_', yr, '.tif'))
      
      #load monthly total precipitation time series
      pr_m <- brick(str_c('wy_projections/temp/pr_monthly_', m, '_', r, '_', yr, '.tif'))
      
      #transform to matrices
      tavg_mat <- getValues(tavg_m)
      pr_mat <- getValues(pr_m)
      
      #transform to degrees C
      tavg_mat <- tavg_mat - 273.15
      
      #create vector of latitudes
      #create blank to fill because of missing values
      lat <- as.numeric(rep(NA, times = NROW(tavg_mat)))
      lat[is.na(tavg_mat[,1]) == F] <- rasterToPoints(tavg_m[[1]], spatial = F) %>% .[,'y']
      
      #calculate potential evapotranspiration
      pet <- sapply(1:NROW(tavg_mat), function(i){
        thornthwaite(tavg_mat[i,], lat[i], na.rm = T)
      })
      
      #transpose rows and columns
      pet <- t(pet)
      
      #calculate pdsi timeseries named dat
      #will manipulate it to be the main data frame for model run
      dat <- sapply(1:NROW(pet), function(i){
        #calculate scpdsi
        x <- pdsi(pr_mat[i,], pet[i,], start = yr) %>% .$X
        
        #if all values zero return NA
        if(sum(x == 0) == length(x)){as.numeric(rep(NA, num_years))} else{
          
          #convert pdsi to matrix of years and months
          x <- matrix(x, nrow = num_years, byrow = T)
          
          #extract the min value from mar and apr
          apply(x, 1, function(x) min(x[3:4]))
        }
      })
      
      #transpose rows and columns
      dat <- t(dat)
      
      #change to data frame
      dat <- as.data.frame(dat)
      
      #change col names to years
      colnames(dat) <- as.numeric(start_yr):as.numeric(end_yr)
      
      #create id col
      dat <- cbind(data.frame(id = 1:NROW(dat), dat))
      
      #melt df
      dat <- reshape2::melt(dat, id.vars = 'id', variable.name = 'year',
                            value.name = 'pdsi')
      
      #fix year column
      dat$year <- as.character(dat$year) %>% 
        str_replace(., 'X', '') %>% 
        as.numeric
      
      #load elevation layer
      dem <- raster('dem/srtm/dem_wy_maca_proj.tif')
      
      #extract and remove raster
      dat <- cbind(dat, dem = raster::getValues(dem) %>% rep(times = num_years) %>% round(2))
      rm(dem)
      
      #calculate rain_mar_may, snow_oct_apr, and gdd_jan_apr
      #load daily max temp, mean temp, and prcp
      tmax <- brick(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*', yr, '*.tif')))
      tavg <- brick(str_c('wy_projections/temp/tavg_', m, '_', r, '_', yr, '.tif'))
      pr <- brick(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*', yr, '*.tif')))
      
      #also need to calculate tmax thresholds to separate rain and snow using 
      #linear regression of tmax vs. elevation best across all ecoregions
      #from Rajagopal and Harpold (2016)
      tmax_thresh <- 3.469 - (1.667*(dat$dem[dat$year == start_yr]/1000))
      
      #create vector of dates to use to extract data
      #create time series of years, months, days
      ts <- seq(mdy(str_c('01-01-', start_yr)), mdy(str_c('12-31-', end_yr)), by = "day") %>%
        format(., format = '%Y%m%d')
      
      #allocate columns for rain and gdd
      dat$rain <- NA
      dat$snow <- NA
      dat$gdd <- NA
      
      #loop through years
      for(y in unique(dat$year)){
        
        #find time series date IDs that match variables of interest
        #need oct-may to cover all variables. only oct-dec if isn't first year
        if(y != start_yr){
          index <- c(which(str_detect(ts, str_c(y-1, "10|", y-1, "11|", y-1, "12"))), 
                          which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04|", y, "05"))))
        } else index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04|", y, "05")))
        
        #extract raster data for time period
        tmax_yr <- getValues(tmax[[index]])
        tavg_yr <- getValues(tavg[[index]])
        pr_yr <- getValues(pr[[index]])
        
        #change temp to celcius
        tmax_yr <- tmax_yr - 273.15
        tavg_yr <- tavg_yr - 273.15
        
        #create date time series for specific year
        ts_yr <- ts[index]
        
        #calc rain from mar-may
        rain <- pr_yr[, which(ts_yr == str_c(y, '0301')):which(ts_yr == str_c(y, '0531'))]
        temp <- tmax_yr[, which(ts_yr == str_c(y, '0301')):which(ts_yr == str_c(y, '0531'))]
        rain[temp <= tmax_thresh] <- 0
        
        #input into dat
        dat$rain[dat$year == y] <- rowSums(rain)
        
        #clean up
        rm(rain, temp)
        
        #calc snow from oct-apr, NA for first year of time-series
        if(y != start_yr){
          snow <- pr_yr[, which(ts_yr == str_c(y-1, '1001')):which(ts_yr == str_c(y, '0430'))]
          temp <- tmax_yr[, which(ts_yr == str_c(y-1, '1001')):which(ts_yr == str_c(y, '0430'))]
          snow[temp > tmax_thresh] <- 0
          
          #input into dat
          dat$snow[dat$year == y] <- rowSums(snow)
          
          #clean up
          rm(snow, temp)
        } 
        
        #set negative gdd to zero. use 5 degree threshold
        tavg_yr <- tavg_yr - 5
        tavg_yr[tavg_yr < 0] <- 0
        
        #calc gdd from jan to apr
        gdd <- tavg_yr[, which(ts_yr == str_c(y, '0101')):which(ts_yr == str_c(y, '0430'))]
        
        #input into dat
        dat$gdd[dat$year == y] <- rowSums(gdd)
        
        #clean up
        rm(gdd)
        
      } #end of loop over years
      
      #calc average values over time period
      dat_avg <- dat %>%
        group_by(id) %>%
        summarise(pdsi = mean(pdsi, na.rm = T), rain = mean(rain, na.rm = T), snow = mean(snow, na.rm = T), 
                  gdd = mean(gdd, na.rm = T))
      
      #change NaN values to NA
      dat_avg[is.na(dat_avg)] <- NA
      
      #write out csvs
      #check after to make sure averages are calculating correctly!!
      write.csv(dat, file = str_c('wy_projections/dat/dat_', m, '_', r, '_', yr, '.csv'),
                row.names = F)
      write.csv(dat_avg, file = str_c('wy_projections/dat/dat_avg_', m, '_', r, '_', yr, '.csv'),
                row.names = F)
      
    } #end of loop over rcps
  } #end of loop over years
} #end of loop over mods

#stop parallel cluster
parallel::stopCluster(cl)

