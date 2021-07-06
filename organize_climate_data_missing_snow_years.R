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
years <- c(2038, 2068)

#scenario periods we downloaded
rcp <- c('rcp45', 'rcp85')

#we need a dem in the climate data project to define the elevation categories
#the thresholds from model building are: 1627 and 2281 m

#years 2018 for the middle model
# mods <- 'HadGEM2-ES365'
# years <- 2018

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
  
  #don't loop through time periods anymore. put them together
  #for(yr in years){
  for(yr in years){
    
    #loop through rcp scenarios
    #for(r in rcp){
    for(r in rcp){
      
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

# ##############################################################
# ###CALCULATE MISSING YEARS OF SNOW DATA FOR MIDDLE SCENARIO###
# ##############################################################
# 
# #only run for middle scenario, two rcps
# m <- 'HadGEM2-ES365'
# #r <- 'rcp45'
# r <- 'rcp85'
# 
# #calculate snow_oct_apr for missing years
# #2020, 2040, 2070
# 
# #create vector of dates to use to extract data
# #create time series of years, months, days
# ts_1999 <- seq(mdy(str_c('01-01-', 1999)), mdy(str_c('12-31-', 2019)), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# ts_2020 <- seq(mdy(str_c('01-01-', 2020)), mdy(str_c('12-31-', 2039)), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# ts_2040 <- seq(mdy(str_c('01-01-', 2040)), mdy(str_c('12-31-', 2069)), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# ts_2070 <- seq(mdy(str_c('01-01-', 2070)), mdy(str_c('12-31-', 2099)), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# #create index of the bands we want for each year combo
# index_1999 <- which(str_detect(ts_1999, str_c(2019, "10|", 2019, "11|", 2019, "12")))
# 
# index_2020_1 <- which(str_detect(ts_2020, str_c(2020, "01|", 2020, "02|", 2020, "03|", 2020, "04|", 2020, "05")))
# 
# index_2020_2 <- which(str_detect(ts_2020, str_c(2039, "10|", 2039, "11|", 2039, "12")))
# 
# index_2040_1 <- which(str_detect(ts_2040, str_c(2040, "01|", 2040, "02|", 2040, "03|", 2040, "04|", 2040, "05")))
# 
# index_2040_2 <- which(str_detect(ts_2040, str_c(2069, "10|", 2069, "11|", 2069, "12")))
# 
# index_2070 <- which(str_detect(ts_2070, str_c(2070, "01|", 2070, "02|", 2070, "03|", 2070, "04|", 2070, "05")))
# 
# #load daily max temp and prcp for the days we want
# tmax_2020 <- stack(stack(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*1999*.tif')), bands = index_1999),
#                    stack(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*2020*.tif')), bands = index_2020_1)) %>%
#   getValues
# 
# pr_2020 <- stack(stack(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*1999*.tif')), bands = index_1999),
#                  stack(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*2020*.tif')), bands = index_2020_1)) %>%
#   getValues
# 
# tmax_2040 <- stack(stack(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*2020*.tif')), bands = index_2020_2),
#                    stack(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*2040*.tif')), bands = index_2040_1)) %>%
#   getValues
# 
# pr_2040 <- stack(stack(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*2020*.tif')), bands = index_2020_2),
#                  stack(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*2040*.tif')), bands = index_2040_1)) %>%
#   getValues
# 
# tmax_2070 <- stack(stack(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*2040*.tif')), bands = index_2040_2),
#                    stack(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*2070*.tif')), bands = index_2070)) %>%
#   getValues
# 
# pr_2070 <- stack(stack(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*2040*.tif')), bands = index_2040_2),
#                  stack(Sys.glob(str_c('wy_projections/raw/pr*', m, '*', r, '*2070*.tif')), bands = index_2070)) %>%
#   getValues
# 
# #set temp to degree C
# tmax_2020 <- tmax_2020 - 273.15
# tmax_2040 <- tmax_2040 - 273.15
# tmax_2070 <- tmax_2070 - 273.15
# 
# #load elevation layer
# dem <- raster('dem/srtm/dem_wy_maca_proj.tif')
# 
# #extract and remove raster
# dem <- raster::getValues(dem) %>% round(2)
# 
# #also need to calculate tmax thresholds to separate rain and snow using 
# #linear regression of tmax vs. elevation best across all ecoregions
# #from Rajagopal and Harpold (2016)
# tmax_thresh <- 3.469 - (1.667*(dem/1000))
# 
# #allocate columns for rain and gdd
# #dat$rain <- NA
# #dat$snow <- NA
# #dat$gdd <- NA
# 
# #need to remove may from time series since over calculated
# pr_2020 <- pr_2020[,1:(ncol(pr_2020)-31)]
# pr_2040 <- pr_2040[,1:(ncol(pr_2040)-31)]
# pr_2070 <- pr_2070[,1:(ncol(pr_2070)-31)]
# tmax_2020 <- tmax_2020[,1:(ncol(tmax_2020)-31)]
# tmax_2040 <- tmax_2040[,1:(ncol(tmax_2040)-31)]
# tmax_2070 <- tmax_2070[,1:(ncol(tmax_2070)-31)]
# 
# #calc snow from oct-apr
# snow_2020 <- pr_2020
# snow_2020[tmax_2020 > tmax_thresh] <- 0
# snow_2040 <- pr_2040
# snow_2040[tmax_2040 > tmax_thresh] <- 0
# snow_2070 <- pr_2070
# snow_2070[tmax_2070 > tmax_thresh] <- 0
# 
# #input into dat
# dat_2020 <- data.frame(snow = rowSums(snow_2020), year = 2020)
# dat_2040 <- data.frame(snow = rowSums(snow_2040), year = 2040)
# dat_2070 <- data.frame(snow = rowSums(snow_2070), year = 2070)
# 
# #write csv
# write.csv(dat_2020, file = str_c('wy_projections/dat/dat_', m, '_', r, '_', 2020, '_snow.csv'),
#           row.names = F)
# write.csv(dat_2040, file = str_c('wy_projections/dat/dat_', m, '_', r, '_', 2040, '_snow.csv'),
#           row.names = F)
# write.csv(dat_2070, file = str_c('wy_projections/dat/dat_', m, '_', r, '_', 2070, '_snow.csv'),
#           row.names = F)

##########################################
###CALCULATE VPD FOR SPRING SCALE MODEL###
##########################################  

#the only additional variable for the SS model is VPD
#vpd_tavg_mean_jan_apr

#we want four year categories for middle model test
#mods <- 'HadGEM2-ES365'
#years <- c('1999', '2020', '2040', '2070')

#remove HadGEM2 model since already calculated
#mods <- mods[!(mods %in% 'HadGEM2-ES365')]

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
  
  for(r in rcp){
    for(yr in years){
      
      #strip start and end year
      start_yr <- str_extract(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')),
                              "[:digit:]{4}_[:digit:]{4}") %>% str_sub(start = 1, end = 4)
      end_yr <- str_extract(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')),
                            "[:digit:]{4}_[:digit:]{4}") %>% str_sub(start = 6, end = 9)
      
      #calc number of years
      num_years <- as.numeric(end_yr) - as.numeric(start_yr) + 1
      
      #load daily max temp, min temp, max rh, min rh
      tmax <- brick(Sys.glob(str_c('wy_projections/raw/tasmax*', m, '*', r, '*', yr, '*.tif')))
      tmin <- brick(Sys.glob(str_c('wy_projections/raw/tasmin*', m, '*', r, '*', yr, '*.tif')))
      rhmax <- brick(Sys.glob(str_c('wy_projections/raw/rhsmax*', m, '*', r, '*', yr, '*.tif')))
      rhmin <- brick(Sys.glob(str_c('wy_projections/raw/rhsmin*', m, '*', r, '*', yr, '*.tif')))
      
      #create vector of dates to use to extract data
      #create time series of years, months, days
      ts <- seq(mdy(str_c('01-01-', start_yr)), mdy(str_c('12-31-', end_yr)), by = "day") %>%
        format(., format = '%Y%m%d')
      
      #create output dataframe
      #id col is 1:ncell(tmax) and for number of years
      dat <- data.frame(id = rep(1:ncell(tmax), times = num_years),
                        year = rep(start_yr:end_yr, each = ncell(tmax)))
      
      #allocate vpd col
      dat$vpd <- NA
      
      #loop through years
      for(y in unique(dat$year)){
        
        #find time series date IDs that match variables of interest
        #need jan-apr
        index <- c(which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04"))))
        
        #extract raster data for time period
        tmax_yr <- getValues(tmax[[index]])
        tmin_yr <- getValues(tmin[[index]])
        rhmax_yr <- getValues(rhmax[[index]])
        rhmin_yr <- getValues(rhmin[[index]])
        
        #change temp to celcius
        tmax_yr <- tmax_yr - 273.15
        tmin_yr <- tmin_yr - 273.15
        
        #calc average rh
        rhavg_yr <- (rhmax_yr + rhmin_yr) / 2
        
        #calculate saturation vapor pressure same method as used for daymet
        svp <- (((610.7 * exp((17.38 * tmax_yr)/(tmax_yr + 239))) + (610.7 * exp((17.38 * tmin_yr)/(tmin_yr + 239)))) / 2)
        
        #calc daily vpd using formula of svp and rh
        vpd <- svp * (1 - (rhavg_yr / 100))
        
        #input into dat
        dat$vpd[dat$year == y] <- round(rowMeans(vpd), 2)
        
        #clean up
        rm(tmax_yr, tmin_yr, rhmax_yr, rhmin_yr, rhavg_yr, svp, vpd, index)
        
      } #end of loop over years
      
      #calc average values over time period
      dat_avg <- dat %>%
        group_by(id) %>%
        summarise(vpd = mean(vpd, na.rm = T))
      
      #change NaN values to NA
      dat_avg[is.na(dat_avg)] <- NA
      
      #write out csvs
      #check after to make sure averages are calculating correctly!!
      write.csv(dat, file = str_c('wy_projections/dat/dat_vpd_', m, '_', r, '_', yr, '.csv'),
                row.names = F)
      write.csv(dat_avg, file = str_c('wy_projections/dat/dat_avg_vpd_', m, '_', r, '_', yr, '.csv'),
                row.names = F)
      
    } #end of loop over rcps
  } #end of loop over years
} #end of loop over mods

#stop parallel cluster
parallel::stopCluster(cl)
