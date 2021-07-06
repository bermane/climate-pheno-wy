#this code looks at PIRGd model outputs and assesses the cause of different values based
#on the input covariates

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(lubridate)

#set wd
setwd('/Volumes/SSD/climate_effects')

#################################
###LOAD AND RESAMPLE LANDCOVER###
#################################

#first we need to resample landcover to include it in the model
#load climate raster template
clim_ras <- raster('wy_projections/raw/pr_wy_projections_CNRM-CM5_r1i1p1_rcp45_macav2metdata_2040_2069_daily.tif',
                   band = 1)

#load landcover map
lc <- raster('landcover/five_class_landcover_wy_laea_2016.tif')
lc[lc == 3] <- NA
lc[lc == 5] <- NA
lc[lc == 4] <- 3

#aggregate and resample landcover map to match projections
lc <- lc %>% raster::aggregate(., fact = 16, fun = modal) %>% 
  projectRaster(., clim_ras, method = 'ngb') %>% crop(., clim_ras)

#extract lc values as vector
lc_v <- getValues(lc)

#load middle and low scenario RCP 8.5 values
dat_gem_2041 <- read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp85_2040.csv') 
dat_gem_2071 <- read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp85_2070.csv')
dat_inmcm_2041 <- read.csv('wy_projections/dat/dat_inmcm4_rcp85_2040.csv')
dat_inmcm_2071 <- read.csv('wy_projections/dat/dat_inmcm4_rcp85_2070.csv')

#only keep year 2041 and 2071 to convert into rasters (as example). 2040/2070 missing snow values
dat_gem_2041 <- dat_gem_2041[dat_gem_2041$year == 2041,]
dat_gem_2071 <- dat_gem_2071[dat_gem_2071$year == 2071,]
dat_inmcm_2041 <- dat_inmcm_2041[dat_inmcm_2041$year == 2041,]
dat_inmcm_2071 <- dat_inmcm_2071[dat_inmcm_2071$year == 2071,]

#add landcover to dataframes
dat_gem_2041$lc <- lc_v
dat_gem_2071$lc <- lc_v
dat_inmcm_2041$lc <- lc_v
dat_inmcm_2071$lc <- lc_v

#load pirgd values for same years/models. these will be base of raster stacks
ras_gem_2041 <- raster('wy_projections/pirgd/pirgd_ann_HadGEM2-ES365_rcp85_2040.tif', band = 2)
ras_gem_2071 <- raster('wy_projections/pirgd/pirgd_ann_HadGEM2-ES365_rcp85_2070.tif', band = 2)
ras_inmcm_2041 <- raster('wy_projections/pirgd/pirgd_ann_inmcm4_rcp85_2040.tif', band = 2)
ras_inmcm_2071 <- raster('wy_projections/pirgd/pirgd_ann_inmcm4_rcp85_2070.tif', band = 2)

#ensure rasters have same number of cells as data frames
ncell(ras_gem_2041)
nrow(dat_gem_2041)

#create raster stack of covariates on top of pirgd
covar <- c('rain', 'snow', 'gdd', 'pdsi', 'lc')

#loop through co-variates
for(vari in covar){
  
  #load raster template
  ras <- clim_ras
  
  #gem_2040
  ras[] <- dat_gem_2041[,vari] %>% round(2)
  ras_gem_2041 <- addLayer(ras_gem_2041, ras)
  
  #gem_2070
  ras[] <- dat_gem_2071[,vari] %>% round(2)
  ras_gem_2071 <- addLayer(ras_gem_2071, ras)
  
  #inmcm_2040
  ras[] <- dat_inmcm_2041[,vari] %>% round(2)
  ras_inmcm_2041 <- addLayer(ras_inmcm_2041, ras)
  
  #inmcm_2070
  ras[] <- dat_inmcm_2071[,vari] %>% round(2)
  ras_inmcm_2071 <- addLayer(ras_inmcm_2071, ras)
  
}

#clean up
rm(ras)

#change layer names
names(ras_gem_2041) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')
names(ras_gem_2071) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')
names(ras_inmcm_2041) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')
names(ras_inmcm_2071) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')

#create output df of values from different pixel locations
pix <- data.frame(model = character(), location = character(), pirgd = numeric(), 
                  rain = numeric(), snow = numeric(), gdd = numeric(), pdsi = numeric(),
                  landcover = numeric())

#click gem_2041
plot(ras_gem_2041[[1]])
vals <- click(ras_gem_2041)

#rbind to df
pix <- rbind(pix,
             data.frame(model = 'HadGEM 2041 RCP 8.5',
                        location = c('SE (high pirgd)', 'SW (low pirgd)', 
                                     'NW (high pirgd)', 'NW (mid pirgd)'),
                        pirgd = vals[,'pirgd'],
                        rain = vals[,'rain'],
                        snow = vals[,'snow'],
                        gdd = vals[,'gdd'],
                        pdsi = vals[,'pdsi'],
                        landcover = vals[,'lc']))

#click gem_2071
plot(ras_gem_2071[[1]])
vals <- click(ras_gem_2071)

#rbind to df
pix <- rbind(pix,
             data.frame(model = 'HadGEM 2071 RCP 8.5',
                        location = c('SE (high pirgd)', 'SW (low pirgd)', 
                                     'NW (high pirgd)', 'NW (mid pirgd)'),
                        pirgd = vals[,'pirgd'],
                        rain = vals[,'rain'],
                        snow = vals[,'snow'],
                        gdd = vals[,'gdd'],
                        pdsi = vals[,'pdsi'],
                        landcover = vals[,'lc']))

#click inmcm_2041
plot(ras_inmcm_2041[[1]])
vals <- click(ras_inmcm_2041)

#rbind to df
pix <- rbind(pix,
             data.frame(model = 'inmcm 2041 RCP 8.5',
                        location = c('SE (high pirgd)', 'SW (low pirgd)', 
                                     'NW (high pirgd)', 'NW (mid pirgd)'),
                        pirgd = vals[,'pirgd'],
                        rain = vals[,'rain'],
                        snow = vals[,'snow'],
                        gdd = vals[,'gdd'],
                        pdsi = vals[,'pdsi'],
                        landcover = vals[,'lc']))

#click inmcm_2071
plot(ras_inmcm_2071[[1]])
vals <- click(ras_inmcm_2071)

#rbind to df
pix <- rbind(pix,
             data.frame(model = 'inmcm 2071 RCP 8.5',
                        location = c('SE (high pirgd)', 'SW (low pirgd)', 
                                     'NW (high pirgd)', 'NW (mid pirgd)'),
                        pirgd = vals[,'pirgd'],
                        rain = vals[,'rain'],
                        snow = vals[,'snow'],
                        gdd = vals[,'gdd'],
                        pdsi = vals[,'pdsi'],
                        landcover = vals[,'lc']))

# export pix as csv
# write.csv(pix, file = 'output/pixel_extremes_2020_12_29.csv',
#          row.names = F, col.names = F)

##############################
###CALCULATE CUMULATIVE GDD###
##############################
  
#gem_2041
#load daily mean temp
tavg <- brick(str_c('wy_projections/temp/tavg_HadGEM2-ES365_rcp85_2040.tif'))

#create vector of dates to use to extract data
#create time series of years, months, days
ts <- seq(mdy('01-01-2040'), mdy('12-31-2069'), by = "day") %>%
  format(., format = '%Y%m%d')

#create index of which bands to extract
y <- 2041
index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))

#subset raster dataset
tavg <- tavg[[index]]

#plot pirgd map
plot(ras_gem_2041[[1]])

#load gdd values from click
vals <- click(tavg)

#from tavg variables, GDD equals tavg - 273.15 - 5
#set negatives to zero
vals <- round(vals - 273.15 - 5, 2)
vals[vals < 0] <- 0

#calc cumsum
vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t

#set plotting
par(mfrow = c(2, 2))

#calculate cumsum
plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))

#reset par
dev.off()

#gem_2071
#load daily mean temp
tavg <- brick(str_c('wy_projections/temp/tavg_HadGEM2-ES365_rcp85_2070.tif'))

#create vector of dates to use to extract data
#create time series of years, months, days
ts <- seq(mdy('01-01-2070'), mdy('12-31-2099'), by = "day") %>%
  format(., format = '%Y%m%d')

#create index of which bands to extract
y <- 2071
index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))

#subset raster dataset
tavg <- tavg[[index]]

#plot pirgd map
plot(ras_gem_2071[[1]])

#load gdd values from click
vals <- click(tavg)

#from tavg variables, GDD equals tavg - 273.15 - 5
#set negatives to zero
vals <- round(vals - 273.15 - 5, 2)
vals[vals < 0] <- 0

#calc cumsum
vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t

#set plotting
par(mfrow = c(2, 2))

#calculate cumsum
plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))

#reset par
dev.off()

#inmcm_2041
#load daily mean temp
tavg <- brick(str_c('wy_projections/temp/tavg_inmcm4_rcp85_2040.tif'))

#create vector of dates to use to extract data
#create time series of years, months, days
ts <- seq(mdy('01-01-2040'), mdy('12-31-2069'), by = "day") %>%
  format(., format = '%Y%m%d')

#create index of which bands to extract
y <- 2041
index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))

#subset raster dataset
tavg <- tavg[[index]]

#plot pirgd map
plot(ras_inmcm_2041[[1]])

#load gdd values from click
vals <- click(tavg)

#from tavg variables, GDD equals tavg - 273.15 - 5
#set negatives to zero
vals <- round(vals - 273.15 - 5, 2)
vals[vals < 0] <- 0

#calc cumsum
vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t

#set plotting
par(mfrow = c(2, 2))

#calculate cumsum
plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))

#reset par
dev.off()

#inmcm_2071
#load daily mean temp
tavg <- brick(str_c('wy_projections/temp/tavg_inmcm4_rcp85_2070.tif'))

#create vector of dates to use to extract data
#create time series of years, months, days
ts <- seq(mdy('01-01-2070'), mdy('12-31-2099'), by = "day") %>%
  format(., format = '%Y%m%d')

#create index of which bands to extract
y <- 2071
index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))

#subset raster dataset
tavg <- tavg[[index]]

#plot pirgd map
plot(ras_inmcm_2071[[1]])

#load gdd values from click
vals <- click(tavg)

#from tavg variables, GDD equals tavg - 273.15 - 5
#set negatives to zero
vals <- round(vals - 273.15 - 5, 2)
vals[vals < 0] <- 0

#calc cumsum
vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t

#set plotting
par(mfrow = c(2, 2))

#calculate cumsum
plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
     main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))

#reset par
dev.off()

################################
###LOOK AT PIRGd BY ELEVATION###
################################

#load pirgd and dem
pirgd <- raster('dlc/mean_maxIRGdate_wy_laea_2001_2018.tif')
dem <- raster('dem/srtm/dem_wy_laea.tif')   

#create dem categores
dem_lo <- dem
dem_lo[dem_lo >= 1645] <- NA

dem_mi <- dem
dem_mi[dem_mi < 1645 | dem_mi >= 2298] <- NA

dem_hi <- dem
dem_hi[dem_hi < 2298] <- NA

#load WY shapefile
wy <- readOGR('reference/wyoming.shp')

#reproject to PIRGd
wy <- spTransform(wy, crs(pirgd))

#create pirgd layers by elevation
pirgd %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
               main = "PIRGd at all elevations") 
plot(wy, add = T)

pirgd %>% mask(., dem_lo) %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
                                   main = "PIRGd below 1645 m")
plot(wy, add = T)

pirgd %>% mask(., dem_mi) %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
                                   main = 'PIRGd between 1645 and 2298 m')
plot(wy, add = T)

pirgd %>% mask(., dem_hi) %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
                                   main = "PIRGd above 2298 m")
plot(wy, add = T)

#plot pirgd by elevation
plot(pirgd[dem < 1645])




