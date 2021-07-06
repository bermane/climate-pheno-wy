#this code predicts future PIRGd using future climate scenario data and the PIRGd model built
#using the run_pirgd_model function
#the first part calculates average PIRGd over the time periods using the average annual climate data
#the second part calculates annual PIRGd over the time periods, which can then be used to calculate
#average pirgd (currently the preferred method in subsequent code)

#load packages
library(rgdal)
library(sp)
library(raster)
library(dynatopmodel)
library(tidyverse)
library(data.table)
library(lme4)
library(tictoc)
library(scales)
library(nlme)
library(berryFunctions)
library(effects)
library(ggplot2)
library(car)
library(performance)
library(doParallel)
library(parallel)
library(foreach)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load source functions
source("reference/HighstatLibV10.R")

#load workspace
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_10_18.RData')

##################################
###REMOVE MISSING/UNWANTED DATA###
##################################
#put the variables with different elevation levels into single columns
dat$rain_elev <- rowSums(data.frame(dat$rain_low_mar_apr, dat$rain_med_mmar_mmay, dat$rain_hi_apr_may), na.rm = T)
dat$gdd_elev <- rowSums(data.frame(dat$gdd_low_apr, dat$gdd_med_mapr_mmay, dat$gdd_hi_may), na.rm = T)
dat$vp_min_elev <- rowSums(data.frame(dat$vp_min_low_mar_apr, dat$vp_min_med_mmar_mmay, dat$vp_min_hi_apr_may), na.rm = T)
dat$vp_avg_elev <- rowSums(data.frame(dat$vp_avg_low_mar_apr, dat$vp_avg_med_mmar_mmay, dat$vp_avg_hi_apr_may), na.rm = T)

#change evergreen to 3 instead of 4
dat$lc[dat$lc == 4] <- 3

#load into new dataframe
dat2 <- dat

#remove columns with missing data due to elevation levels
dat2 <- select(dat2, -c(rain_low_mar_apr, rain_med_mmar_mmay, rain_hi_apr_may,
                        gdd_low_apr, gdd_med_mapr_mmay, gdd_hi_may,
                        vp_min_low_mar_apr, vp_min_med_mmar_mmay, vp_min_hi_apr_may,
                        vp_avg_low_mar_apr, vp_avg_med_mmar_mmay, vp_avg_hi_apr_may))

#remove all veg layers since we are moving forward without them
dat2 <- select(dat2, -c(herb_homer_ann, sage_homer_ann, shrub_homer_ann, herb_homer,
                        sage_homer, shrub_homer, ann_forb_rap, bare_ground_rap,
                        perenn_forb_rap, shrub_rap, tree_rap, ann_perenn_forb_rap))

#remove vapor pressure variables
dat2 <- select(dat2, -c(vp_min_elev, vp_avg_elev, vp_min_jan_apr, vp_avg_jan_apr))

#remove time1-PIRGd variables
dat2 <- select(dat2, -c(solar_jan_pirgd, rain_oct_pirgd, rain_jan_pirgd,
                        gdd_jan_pirgd, pdsi_jan_pirgd_mean, pdsi_jan_pirgd_min,
                        pdsi_jan_pirgd_med))

#remove vpd by elev. keep in values from jan-apr
dat2 <- select(dat2, -c(vpd_tmax_mean_elev, vpd_tmin_mean_elev, vpd_tavg_mean_elev))


#remove any remaining missing data points for PIRGd or snowmelt
#2 points for pirgd
dat2 <- dat2[is.na(dat2$pirgd) == 0,]
#dat2 <- dat2[is.na(dat2$ss) == 0,]
#200 points for snowmelt
dat2 <- dat2[is.na(dat2$snowmelt) == 0,]

#set zero and decimal gdd and rain values to 1 from variables we will log transform
dat2$gdd_jan_apr[dat2$gdd_jan_apr < 1] <- 1
dat2$rain_mar_may[dat2$rain_mar_may < 1] <- 1

#any other missing points?
#may need to remove the variables that include snowmelt and pirgd...
#16% of samples have PIRGd before snowmelt...
sum(is.na(dat2))

#id col and lc as factors
dat2$id <- as.factor(dat2$id)
dat2$lc <- as.factor(dat2$lc)
levels(dat2$lc) <- c('shrub', 'herb', 'evergreen')

###########################
###RUN FINAL PIRGd MODEL###
###########################

#run model without rescaled values
m_pirgd <- lme(pirgd ~ lc + log(gdd_jan_apr) +
                    poly(pdsi_mar_apr_min, 2, raw = TRUE) + log(rain_mar_may) + 
                    snow_oct_apr +
                    log(rain_mar_may):snow_oct_apr +
                    poly(pdsi_mar_apr_min, 2, raw = TRUE):lc, 
                  random = ~ 1 | id, data = dat2)

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

###################################################
###LOAD AVERAGE CLIMATE DATA AND RUN PREDICTIONS###
###################################################

#SUBSET LIST OF FILES TO RUN PREDICTIONS ON

# #load list of climate data csv files
# clim_files <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat_avg*'),
#                          full.names = T)
# 
# #load file names only
# clim_names <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat_avg*'),
#                          full.names = F)

#load list of climate data csv files
clim_files <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat_avg*2068*'),
                         full.names = T)

#load file names only
clim_names <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat_avg*2068*'),
                         full.names = F)

#remove files with vpd in name
clim_files <- clim_files[str_detect(clim_files, 'vpd') == F]
clim_names <- clim_names[str_detect(clim_names, 'vpd') == F]

#register parallel backend
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

#loop through files to generate predictions
foreach::foreach(i = 1:length(clim_files)) %dopar% {
  
  #load packages
  library(rgdal)
  library(sp)
  library(raster)
  library(dynatopmodel)
  library(tidyverse)
  library(data.table)
  library(lme4)
  library(tictoc)
  library(scales)
  library(nlme)
  library(berryFunctions)
  library(effects)
  library(ggplot2)
  library(car)
  library(performance)
  
  #load sample climate data
  dat_clim <- read.csv(clim_files[i])
  
  #remove id column
  dat_clim <- dat_clim[, -(colnames(dat_clim) %in% "id")]
  
  #add lc column and factor
  dat_clim$lc <- lc_v
  dat_clim$lc <- as.factor(dat_clim$lc)
  levels(dat_clim$lc) <- c('shrub', 'herb', 'evergreen')
  
  #recode column names to match model variables
  index <- c('gdd_jan_apr', 'pdsi_mar_apr_min', 'rain_mar_may', 'lc', 'snow_oct_apr')
  names(index) <- c('gdd', 'pdsi', 'rain', 'lc', 'snow')
  colnames(dat_clim) <- dplyr::recode(colnames(dat_clim), !!!index)
  
  #apply over rows to handle NA's
  pred <- sapply(1:NROW(dat_clim), function(x){
    if(sum(is.na(dat_clim[x,])) == 0){
      predict(m_pirgd, dat_clim[x,], level = 0)
    } else{NA}
  })
  
  #create output raster of pirgd values
  pred_ras <- clim_ras
  pred_ras[] <- round(pred)
  
  #create output name
  out_name <- clim_names[i] %>% str_replace(., '.csv', '.tif') %>%
    str_replace(., 'dat', 'pirgd')
  
  #write raster
  writeRaster(pred_ras, filename = str_c('wy_projections/pirgd_log_gdd_log_rain/', out_name),
              format = 'GTiff')
  
  #clean up
  rm(pred_ras, out_name, index, dat_clim, pred)
  
}

#stop parallel cluster
parallel::stopCluster(cl)

#################################################
###LOAD DAILY CLIMATE DATA AND RUN PREDICTIONS###
#################################################

# #load list of daily climate data csv files
# clim_files <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat*'),
#                          full.names = T)
# 
# #load file names only
# clim_names <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat*'),
#                          full.names = F)

#load list of daily climate data csv files
clim_files <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat*2068*'),
                         full.names = T)

#load file names only
clim_names <- list.files(path = 'wy_projections/dat', pattern = glob2rx('dat*2068*'),
                         full.names = F)

#remove files with vpd in name
clim_files <- clim_files[str_detect(clim_files, 'vpd') == F]
clim_names <- clim_names[str_detect(clim_names, 'vpd') == F]

#remove files with average data only
clim_files <- clim_files[str_detect(clim_files, 'dat_avg') == F]
clim_names <- clim_names[str_detect(clim_names, 'dat_avg') == F]

#remove files with snow names
clim_files <- clim_files[str_detect(clim_files, 'snow') == F]
clim_names <- clim_names[str_detect(clim_names, 'snow') == F]

#register parallel backend
cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)

#loop through files to generate predictions
foreach::foreach(i = 1:length(clim_files)) %dopar% {
  
  #load packages
  library(rgdal)
  library(sp)
  library(raster)
  library(dynatopmodel)
  library(tidyverse)
  library(data.table)
  library(lme4)
  library(tictoc)
  library(scales)
  library(nlme)
  library(berryFunctions)
  library(effects)
  library(ggplot2)
  library(car)
  library(performance)
  
  #load sample climate data
  dat_clim <- read.csv(clim_files[i])
  
  #if middle model load missing years of snow data
  if(str_detect(clim_files[i], 'HadGEM2') & str_detect(clim_files[i], '2020|2040|2070')){
    
    #ensure loading correct year
    if(str_detect(clim_files[i], '2020')) y <- 2020
    if(str_detect(clim_files[i], '2040')) y <- 2040
    if(str_detect(clim_files[i], '2070')) y <- 2070
    
    #ensure loading correct rcp
    if(str_detect(clim_files[i], 'rcp45')) rcp <- 'rcp45'
    if(str_detect(clim_files[i], 'rcp85')) rcp <- 'rcp85'
    
    #load snow file
    snow <- read.csv(str_c('wy_projections/dat/dat_HadGEM2-ES365_', rcp, '_', y, '_snow.csv'))
    
    #add to dat_clim
    dat_clim$snow[dat_clim$year == y] <- snow$snow
    
  }
  
  #remove id column
  dat_clim <- dat_clim[, -(colnames(dat_clim) %in% "id")]
  
  #add lc column and factor
  dat_clim$lc <- rep(lc_v, times = NROW(dat_clim)/length(lc_v))
  dat_clim$lc <- as.factor(dat_clim$lc)
  levels(dat_clim$lc) <- c('shrub', 'herb', 'evergreen')
  
  #recode column names to match model variables
  index <- c('gdd_jan_apr', 'pdsi_mar_apr_min', 'rain_mar_may', 'lc', 'snow_oct_apr')
  names(index) <- c('gdd', 'pdsi', 'rain', 'lc', 'snow')
  colnames(dat_clim) <- dplyr::recode(colnames(dat_clim), !!!index)
  
  #apply over rows to handle NA's
  pred <- sapply(1:NROW(dat_clim), function(x){
    if(sum(is.na(dat_clim[x,])) == 0){
      predict(m_pirgd, dat_clim[x,], level = 0)
    } else{NA}
  })
  
  #create output raster of pirgd values
  pred_ras <- raster(clim_ras)
  
  #loop through years to add layers
  for(j in 1:length(unique(dat_clim$year))){
    if(j == 1){
      pred_ras[] <- round(pred[min(which(dat_clim$year == unique(dat_clim$year)[j])):
                                          max(which(dat_clim$year == unique(dat_clim$year)[j]))])
    }else{
      ras <- raster(clim_ras)
      ras[] <- round(pred[min(which(dat_clim$year == unique(dat_clim$year)[j])):
                            max(which(dat_clim$year == unique(dat_clim$year)[j]))])
      pred_ras <- addLayer(pred_ras, ras)
      rm(ras)
    }
  }
  
  #create output name
  out_name <- clim_names[i] %>% str_replace(., '.csv', '.tif') %>%
    str_replace(., 'dat', 'pirgd_ann')
  
  #write raster
  writeRaster(pred_ras, filename = str_c('wy_projections/pirgd_log_gdd_log_rain/', out_name),
              format = 'GTiff')
  
  #clean up
  rm(pred_ras, out_name, index, dat_clim, pred)
  
}

#stop parallel cluster
parallel::stopCluster(cl)

###############################################
###PROCESS MEAN PIRGd BASED ON ANNUAL LAYERS###
###############################################

#create file list of pirgd projection rasters
files <- list.files('wy_projections/pirgd_log_gdd_log_rain', glob2rx('pirgd_ann*.tif'),
                    full.names = T)
names <- list.files('wy_projections/pirgd_log_gdd_log_rain', glob2rx('pirgd_ann*.tif'))

#loop through projection rasters
for(i in 1:length(files)){
  
  #load projection raster, calculate rounded mean, and write file
  brick(files[i]) %>% calc(., function(x) round(mean(x, na.rm = T))) %>% 
    writeRaster(., filename = str_replace(files[i], 'pirgd_ann', 'pirgd_mean'),
                           format = 'GTiff', overwrite = T)
  
}

