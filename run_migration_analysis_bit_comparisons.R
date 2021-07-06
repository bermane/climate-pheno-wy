#this code combines all the processed migration data
#and runs the greenscape analysis with the PIRGd model outputs (log rain and log gdd)
#and contemporary data

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(lubridate)
library(scatterplot3d)
library(reshape2)
library(ggcorrplot)
library(gmodels)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load WY shapefile
wy <- readOGR('reference/wyoming.shp')

######################################
###LOAD PIRGd AND SPRING SCALE DATA###
######################################

#load projection data
#create file list of pirgd mean projection rasters
files <- list.files('wy_projections/pirgd_log_gdd_log_rain', 'pirgd_mean.*tif$',
                    full.names = T)
files_ss <- list.files('wy_projections/ss', 'ss_mean.*tif$',
                       full.names = T)

#remove files with 2018, 2020, 2038, 2068 annual date
files <- files[str_detect(files, '2020') == F]
files_ss <- files_ss[str_detect(files_ss, '2020') == F]
files <- files[str_detect(files, '2038') == F]
files_ss <- files_ss[str_detect(files_ss, '2038') == F]
files <- files[str_detect(files, '2068') == F]
files_ss <- files_ss[str_detect(files_ss, '2068') == F]
files <- files[str_detect(files, '2018') == F]
files_ss <- files_ss[str_detect(files_ss, '2018') == F]

#don't remove files with 1999 because we do want BIT
#files <- files[str_detect(files, '1999') == F]

#load into list. names are already present as file names
proj <- lapply(files, raster)
proj_ss <- lapply(files_ss, raster)

#load annual projection data from middle scenario into a list (two scenarios)
mid_1 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_1999.tif',
               'wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_2020.tif',
               'wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_2040.tif',
               'wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_2070.tif')

mid_2 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_1999.tif',
               'wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_2020.tif',
               'wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_2040.tif',
               'wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_2070.tif')

mid_1_ss <- stack('wy_projections/ss/ss_ann_HadGEM2-ES365_rcp45_1999.tif',
               'wy_projections/ss/ss_ann_HadGEM2-ES365_rcp45_2020.tif',
               'wy_projections/ss/ss_ann_HadGEM2-ES365_rcp45_2040.tif',
               'wy_projections/ss/ss_ann_HadGEM2-ES365_rcp45_2070.tif')

mid_2_ss <- stack('wy_projections/ss/ss_ann_HadGEM2-ES365_rcp85_1999.tif',
               'wy_projections/ss/ss_ann_HadGEM2-ES365_rcp85_2020.tif',
               'wy_projections/ss/ss_ann_HadGEM2-ES365_rcp85_2040.tif',
               'wy_projections/ss/ss_ann_HadGEM2-ES365_rcp85_2070.tif')
  

proj_mid <- list(mid_1, mid_2)
proj_mid_ss <- list(mid_1_ss, mid_2_ss)
rm(mid_1, mid_2, mid_1_ss, mid_2_ss)

#load contemporary data and disturbance mask.
#we need mean baseline data along with annual data for the annual projections
pirgd_ann <- brick('dlc/maxIRGdate_wy_laea_2000_2019.tif')
#pirgd <- calc(pirgd_ann, function(x) mean(x, na.rm = T))
#writeRaster(round(pirgd), 'dlc/mean_maxIRGdate_wy_laea_2000_2019.tif', format = 'GTiff', overwrite = T)
pirgd <- raster('dlc/mean_maxIRGdate_wy_laea_2000_2019.tif')
mask <- raster('dlc/maxIRGdate_mask_disturbances.tif')
ss_ann <- brick('dlc/springScale_wy_laea_2000_2019.tif')
ss <- raster('dlc/mean_springScale_wy_laea_2000_2019.tif')

#resample to match scale of projection data
pirgd <- pirgd %>% raster::aggregate(., fact = 16) %>% 
  projectRaster(., to = proj[[1]]) %>% round

# writeRaster(pirgd, filename = 'dlc/mean_maxIRGdate_wy_maca_2000_2019.tif',
#             format = 'GTiff')

pirgd_ann <- pirgd_ann %>% raster::aggregate(., fact = 16) %>% 
  projectRaster(., to = proj[[1]]) %>% round

# writeRaster(pirgd_ann, filename = 'dlc/maxIRGdate_wy_maca_2000_2019.tif',
#             format = 'GTiff')

ss <- ss %>% raster::aggregate(., fact = 16) %>% 
  projectRaster(., to = proj[[1]]) %>% round(2)

# writeRaster(ss, filename = 'dlc/mean_springScale_wy_maca_2000_2019.tif',
#             format = 'GTiff')

ss_ann <- ss_ann %>% raster::aggregate(., fact = 16) %>% 
  projectRaster(., to = proj[[1]]) %>% round(2)

# writeRaster(ss_ann, filename = 'dlc/springScale_wy_maca_2000_2019.tif',
#             format = 'GTiff')

#what if we don't resample just change the projection??
# pirgd <- pirgd %>% projectRaster(., crs = crs(proj[[1]])) %>% round
# pirgd_ann <- pirgd_ann %>% projectRaster(., crs = crs(proj[[1]])) %>% round
# mask <- mask %>% projectRaster(., crs = crs(proj[[1]]))

#do the same with the mask but set to NA if more than half the cells missing
#first set all mask values to 1
mask[is.na(mask) == F] <- 1

#aggregate using sum
mask <- mask %>% raster::aggregate(., fact = 16, fun = sum)

#set to na if value is less than a quarter(?) of original cell count
mask[mask < 64] <- NA

#reproject
mask <- mask %>% projectRaster(., to = proj[[1]]) %>% round

#mask disturbances
base <- mask(pirgd, mask)
base_ann <- mask(pirgd_ann, mask)
base_ss <- mask(ss, mask)
base_ann_ss <- mask(ss_ann, mask)
rm(pirgd, pirgd_ann, ss, ss_ann, mask)

#################################
###LOAD MIGRATION SEGMENT DATA###
#################################

#create file list of start and end points
files <- list.files('mig_data/processed', pattern = 'seg_start_end', full.names = T)

#load start and end points from each herd and rbind together
for(i in 1:length(files)){
  
  #load individual .rdata files with points
  load(files[i])
  
  if(i == 1){
    #if first, create new spdf object
    mig_dur <- mig_start_end
  }else {
    #reload crs so no errors
    crs(mig_start_end) <- crs(mig_dur)
    
    #rbind data to spdf object
    mig_dur <- rbind(mig_dur, mig_start_end)
  }
  rm(mig_start_end)
}

#create file list of migration segments
files <- list.files('mig_data/processed', pattern = 'mig_samples', full.names = T)

#load start and end points from each herd and rbind together
for(i in 1:length(files)){
  
  #load individual .rdata files with points
  load(files[i])
  
  if(i == 1){
    #if first, create new spdf object
    mig_seg <- samples
  }else {
    #reload crs so no errors
    crs(samples) <- crs(mig_seg)
    
    #rbind data to spdf object
    mig_seg <- rbind(mig_seg, samples)
  }
  rm(samples)
}

#clean up
rm(files, i)

##########################################################
###CALCULATE GREENSCAPE METRICS FOR MID AND END CENTURY###
##########################################################

#DURATION
#first allocate duration df
dur <- data.frame(global_id = character(), model = character(), duration = numeric())

#first calc dur for base scenario
#create df with extracted results and calc duration
#set extreme values to 1 or 365
dur_base <- data.frame(global_id = mig_dur$GlobalID, pirgd = raster::extract(base, mig_dur)) 
dur_base$pirgd[dur_base$pirgd <= 0] <- 1
dur_base$pirgd[dur_base$pirgd >= 366] <- 365
dur_base <- dur_base %>% group_by(global_id) %>% mutate(duration = pirgd - lag(pirgd))

#keep only even rows
dur_base <- dur_base[seq(2,NROW(dur_base),2),]

#rbind to duration df
dur <- rbind(dur, data.frame(global_id = dur_base$global_id,
                             model = 'base',
                             duration = dur_base$duration))
rm(dur_base)

#loop through all projection data
for(i in 1:length(proj)){
  
  #select data
  x <- proj[[i]]
  
  #create hold df
  hold <- data.frame(global_id = mig_dur$GlobalID,
                     model = str_replace(names(x), 'pirgd_mean_', ''),
                     pirgd = raster::extract(x, mig_dur))
  hold$pirgd[hold$pirgd <= 0] <- 1
  hold$pirgd[hold$pirgd >= 366] <- 365
  hold <- hold %>% group_by(global_id) %>% mutate(duration = pirgd - lag(pirgd))
  
  #only keep even rows
  hold <- hold[seq(2,NROW(hold),2),]
  
  #remove pirgd col
  hold <- hold[, !(colnames(hold) %in% 'pirgd')]
  
  #rbind to duration df
  dur <- rbind(dur, hold)
  
  #clean up
  rm(x, hold)
}

#ORDER
#first allocate order df
order <- data.frame(global_id = character(), model = character(), order = numeric(), n = numeric())

#first calc order for base scenario
#create df with extracted results
#set extreme values to 1 or 365
order_base <- data.frame(global_id = mig_seg$GlobalID,
                         pirgd = raster::extract(base, mig_seg))
order_base$pirgd[order_base$pirgd <= 0] <- 1
order_base$pirgd[order_base$pirgd >= 366] <- 365

#loop through global ids
for(i in unique(order_base$global_id)){
  
  #load data
  hold <- order_base[order_base$global_id == i,]
  
  #remove any rows with missing values
  hold <- hold[is.na(hold$pirgd) == F,]
  
  #add column for pirgd order
  hold$pirgd_id <- rank(hold$pirgd, ties.method = 'first')
  
  if(NROW(hold) < 3){
    #rbind to order df
    order <- rbind(order, data.frame(global_id = i,
                                     model = 'base',
                                     order = NA,
                                     n = NROW(hold)))
  }else {
    order <- rbind(order, data.frame(global_id = i,
                                     model = 'base',
                                     order = cor.test(1:NROW(hold), hold$pirgd_id, method = 'spearman')$estimate,
                                     n = NROW(hold)))
  }
  
  rm(hold)
}

#clean up
rm(i, order_base)

#loop through all projection data
for(i in 1:length(proj)){
  
  #select data
  x <- proj[[i]]
  
  #create hold df
  order_hold <- data.frame(global_id = mig_seg$GlobalID,
                     model = str_replace(names(x), 'pirgd_mean_', ''),
                     pirgd = raster::extract(x, mig_seg))
  order_hold$pirgd[order_hold$pirgd <= 0] <- 1
  order_hold$pirgd[order_hold$pirgd >= 366] <- 365
  
  #loop through global ids
  for(j in unique(order_hold$global_id)){
    
    #load data
    hold <- order_hold[order_hold$global_id == j,]
    
    #remove any rows with missing values
    hold <- hold[is.na(hold$pirgd) == F,]
    
    #add column for pirgd order
    hold$pirgd_id <- rank(hold$pirgd, ties.method = 'first')
    
    if(NROW(hold) < 3){
      #rbind to order df
      order <- rbind(order, data.frame(global_id = j,
                                       model = hold$model[1],
                                       order = NA,
                                       n = NROW(hold)))
    }else {
      #rbind to order df
      order <- rbind(order, data.frame(global_id = j,
                                       model = hold$model[1],
                                       order = cor.test(1:NROW(hold), hold$pirgd_id, method = 'spearman')$estimate,
                                       n = NROW(hold)))
    }

    rm(hold)
  }
  
  #clean up
  rm(order_hold, x)

}

#RATE
#first allocate rate df
rate <- data.frame(global_id = character(), model = character(), rate = numeric(), n = numeric())

#first calc rate for base scenario
#create df with extracted results
rate_base <- data.frame(global_id = mig_seg$GlobalID,
                         ss = raster::extract(base_ss, mig_seg))

#loop through global ids
for(i in unique(rate_base$global_id)){
  
  #load data
  hold <- rate_base[rate_base$global_id == i,]
  
  #remove any rows with missing values
  hold <- hold[is.na(hold$ss) == F,]
  
  if(NROW(hold) < 3){
    #rbind to rate df
    rate <- rbind(rate, data.frame(global_id = i,
                                     model = 'base',
                                     rate = NA,
                                     n = NROW(hold)))
  }else {
    rate <- rbind(rate, data.frame(global_id = i,
                                     model = 'base',
                                     rate = round(mean(hold$ss), 2),
                                     n = NROW(hold)))
  }
  
  rm(hold)
}

#clean up
rm(i, rate_base)

#loop through all projection data
for(i in 1:length(proj_ss)){
  
  #select data
  x <- proj_ss[[i]]
  
  #create hold df
  rate_hold <- data.frame(global_id = mig_seg$GlobalID,
                           model = str_replace(names(x), 'ss_mean_', ''),
                           ss = raster::extract(x, mig_seg))
  
  #loop through global ids
  for(j in unique(rate_hold$global_id)){
    
    #load data
    hold <- rate_hold[rate_hold$global_id == j,]
    
    #remove any rows with missing values
    hold <- hold[is.na(hold$ss) == F,]
    
    if(NROW(hold) < 3){
      #rbind to rate df
      #rbind to rate df
      rate <- rbind(rate, data.frame(global_id = j,
                                     model = hold$model[1],
                                     rate = NA,
                                     n = NROW(hold)))
    }else {
      #rbind to rate df
      rate <- rbind(rate, data.frame(global_id = j,
                                     model = hold$model[1],
                                     rate = round(mean(hold$ss), 2),
                                     n = NROW(hold)))
    }
    

    rm(hold)
  }
  
  #clean up
  rm(rate_hold, x)
  
}

#############################################################
###CALCULATE GREENSCAPE METRICS ANNUAL FOR MIDDLE SCENARIO###
#############################################################

#DURATION
#first allocate annual duration df
dur_ann <- data.frame(global_id = character(), 
                      model = character(),
                      rcp = character(),
                      year = numeric(), 
                      duration = numeric())

#first calc dur for base scenario
#loop through base years and add to dur_ann
for(i in 1:nlayers(base_ann)){
  dur_base <- data.frame(global_id = mig_dur$GlobalID, year = 1999 + i, 
                               pirgd = raster::extract(base_ann[[i]], mig_dur))
  
  #set extreme values to 1 or 365 and calc duration
  dur_base$pirgd[dur_base$pirgd <= 0] <- 1
  dur_base$pirgd[dur_base$pirgd >= 366] <- 365
  dur_base <- dur_base %>% group_by(global_id) %>% mutate(duration = pirgd - lag(pirgd))
  
  #keep only even rows
  dur_base <- dur_base[seq(2,NROW(dur_base),2),]
  
  #rbind to annual duration df
  dur_ann <- rbind(dur_ann, data.frame(global_id = dur_base$global_id,
                                       model = 'base',
                                       rcp = 'baseline',
                                       year = dur_base$year,
                                       duration = dur_base$duration))
  rm(dur_base)
}

#go through two scenarios of proj data
#loop through all years to add data
for(i in 1:nlayers(proj_mid[[1]])){
  
  #add to hold df
  hold_45 <- data.frame(global_id = mig_dur$GlobalID,
                        model = 'HadGEM2-ES365',
                        rcp = 'rcp45',
                        year = 1998 + i,
                        pirgd = raster::extract(proj_mid[[1]][[i]], mig_dur))
  
  hold_85 <- data.frame(global_id = mig_dur$GlobalID,
                        model = 'HadGEM2-ES365',
                        rcp = 'rcp85',
                        year = 1998 + i,
                        pirgd = raster::extract(proj_mid[[2]][[i]], mig_dur))  
  
  #reset extreme values and calculate duration
  hold_45$pirgd[hold_45$pirgd <= 0] <- 1
  hold_45$pirgd[hold_45$pirgd >= 366] <- 365
  hold_45 <- hold_45 %>% group_by(global_id) %>% mutate(duration = pirgd - lag(pirgd))
  
  hold_85$pirgd[hold_85$pirgd <= 0] <- 1
  hold_85$pirgd[hold_85$pirgd >= 366] <- 365
  hold_85 <- hold_85 %>% group_by(global_id) %>% mutate(duration = pirgd - lag(pirgd))

  #only keep even rows
  hold_45 <- hold_45[seq(2,NROW(hold_45),2),]
  hold_85 <- hold_85[seq(2,NROW(hold_85),2),]
  
  #remove pirgd col
  hold_45 <- hold_45[, !(colnames(hold_45) %in% 'pirgd')]
  hold_85 <- hold_85[, !(colnames(hold_85) %in% 'pirgd')]
  
  #rbind to duration df
  dur_ann <- rbind(dur_ann, hold_45, hold_85)
  
  #clean up
  rm(hold_45, hold_85)
}

#ORDER
#first allocate order df
order_ann <- data.frame(global_id = character(), 
                    model = character(),
                    rcp = character(),
                    year = numeric(),
                    order = numeric(), 
                    n = numeric())

#first calc order for base scenario
#loop through years
for(i in 1:nlayers(base_ann)){
  #create df with extracted results
  order_base <- data.frame(global_id = mig_seg$GlobalID,
                           year = 1999 + i,
                           pirgd = raster::extract(base_ann[[i]], mig_seg))
  
  #set extreme values to 1 or 365
  order_base$pirgd[order_base$pirgd <= 0] <- 1
  order_base$pirgd[order_base$pirgd >= 366] <- 365
  
  #loop through global ids
  for(j in unique(order_base$global_id)){
    
    #load data
    hold <- order_base[order_base$global_id == j,]
    
    #remove any rows with missing values
    hold <- hold[is.na(hold$pirgd) == F,]
    
    #add column for pirgd order
    hold$pirgd_id <- rank(hold$pirgd, ties.method = 'first')
    
    if(NROW(hold) < 3){
      #rbind to order df
      order_ann <- rbind(order_ann, data.frame(global_id = j,
                                               model = 'base',
                                               rcp = 'baseline',
                                               year = hold$year[1],
                                               order = NA,
                                               n = NROW(hold)))
    }else{
      #rbind to order df
      order_ann <- rbind(order_ann, data.frame(global_id = j,
                                               model = 'base',
                                               rcp = 'baseline',
                                               year = hold$year[1],
                                               order = cor.test(1:NROW(hold), hold$pirgd_id, method = 'spearman')$estimate,
                                               n = NROW(hold)))
    }
    
    rm(hold)
  }
  rm(order_base)
}

#go through two scenarios projection data
#loop through all years
for(i in 1:nlayers(proj_mid[[1]])){
  #create hold dfs
  order_hold_45 <- data.frame(global_id = mig_seg$GlobalID,
                           model = 'HadGEM2-ES365',
                           rcp = 'rcp45',
                           year = 1998 + i,
                           pirgd = raster::extract(proj_mid[[1]][[i]], mig_seg))
  
  order_hold_85 <- data.frame(global_id = mig_seg$GlobalID,
                              model = 'HadGEM2-ES365',
                              rcp = 'rcp85',
                              year = 1998 + i,
                              pirgd = raster::extract(proj_mid[[2]][[i]], mig_seg))
  
  #remove extreme values
  order_hold_45$pirgd[order_hold_45$pirgd <= 0] <- 1
  order_hold_45$pirgd[order_hold_45$pirgd >= 366] <- 365
  
  order_hold_85$pirgd[order_hold_85$pirgd <= 0] <- 1
  order_hold_85$pirgd[order_hold_85$pirgd >= 366] <- 365
  
  #loop through global ids
  for(j in unique(order_hold_45$global_id)){
    
    #load data
    hold_45 <- order_hold_45[order_hold_45$global_id == j,]
    hold_85 <- order_hold_85[order_hold_85$global_id == j,]
    
    #remove any rows with missing values
    hold_45 <- hold_45[is.na(hold_45$pirgd) == F,]
    hold_85 <- hold_85[is.na(hold_85$pirgd) == F,]
    
    #add column for pirgd order
    hold_45$pirgd_id <- rank(hold_45$pirgd, ties.method = 'first')
    hold_85$pirgd_id <- rank(hold_85$pirgd, ties.method = 'first')
    
    #rbind to order df
    order_ann <- rbind(order_ann, 
                       data.frame(global_id = j,
                                  model = hold_45$model[1],
                                  rcp = hold_45$rcp[1],
                                  year = hold_45$year[1],
                                  order = if(NROW(hold_45) >= 3) {cor.test(1:NROW(hold_45), hold_45$pirgd_id, method = 'spearman')$estimate} else NA,
                                  n = NROW(hold_45)),
                       data.frame(global_id = j,
                                  model = hold_85$model[1],
                                  rcp = hold_85$rcp[1],
                                  year = hold_85$year[1],
                                  order = if(NROW(hold_85) >= 3) {cor.test(1:NROW(hold_85), hold_85$pirgd_id, method = 'spearman')$estimate} else NA,
                                  n = NROW(hold_85)))
    rm(hold_45, hold_85)
  }
  
  #clean up
  rm(order_hold_45, order_hold_85)
}

#RATE
#first allocate order df
rate_ann <- data.frame(global_id = character(), 
                        model = character(),
                        rcp = character(),
                        year = numeric(),
                        rate = numeric(), 
                        n = numeric())

#first calc order for base scenario
#loop through years
for(i in 1:nlayers(base_ann_ss)){
  #create df with extracted results
  rate_base <- data.frame(global_id = mig_seg$GlobalID,
                           year = 1999 + i,
                           ss = raster::extract(base_ann_ss[[i]], mig_seg))
  
  #loop through global ids
  for(j in unique(rate_base$global_id)){
    
    #load data
    hold <- rate_base[rate_base$global_id == j,]
    
    #remove any rows with missing values
    hold <- hold[is.na(hold$ss) == F,]
    
    if(NROW(hold) < 3){
      #rbind to rate df
      rate_ann <- rbind(rate_ann, data.frame(global_id = j,
                                               model = 'base',
                                               rcp = 'baseline',
                                               year = hold$year[1],
                                               rate = NA,
                                               n = NROW(hold)))
    }else{
      #rbind to rate df
      rate_ann <- rbind(rate_ann, data.frame(global_id = j,
                                               model = 'base',
                                               rcp = 'baseline',
                                               year = hold$year[1],
                                               rate = round(mean(hold$ss), 2),
                                               n = NROW(hold)))
    }
    
    rm(hold)
  }
  rm(rate_base)
}

#go through two scenarios projection data
#loop through all years
for(i in 1:nlayers(proj_mid_ss[[1]])){
  #create hold dfs
  rate_hold_45 <- data.frame(global_id = mig_seg$GlobalID,
                              model = 'HadGEM2-ES365',
                              rcp = 'rcp45',
                              year = 1998 + i,
                              ss = raster::extract(proj_mid_ss[[1]][[i]], mig_seg))
  
  rate_hold_85 <- data.frame(global_id = mig_seg$GlobalID,
                              model = 'HadGEM2-ES365',
                              rcp = 'rcp85',
                              year = 1998 + i,
                              ss = raster::extract(proj_mid_ss[[2]][[i]], mig_seg))
  
  #loop through global ids
  for(j in unique(rate_hold_45$global_id)){
    
    #load data
    hold_45 <- rate_hold_45[rate_hold_45$global_id == j,]
    hold_85 <- rate_hold_85[rate_hold_85$global_id == j,]
    
    #remove any rows with missing values
    hold_45 <- hold_45[is.na(hold_45$ss) == F,]
    hold_85 <- hold_85[is.na(hold_85$ss) == F,]
    
    #rbind to rate df
    rate_ann <- rbind(rate_ann, 
                       data.frame(global_id = j,
                                  model = hold_45$model[1],
                                  rcp = hold_45$rcp[1],
                                  year = hold_45$year[1],
                                  rate = round(mean(hold_45$ss), 2),
                                  n = NROW(hold_45)),
                       data.frame(global_id = j,
                                  model = hold_85$model[1],
                                  rcp = hold_85$rcp[1],
                                  year = hold_85$year[1],
                                  rate = round(mean(hold_85$ss), 2),
                                  n = NROW(hold_85)))
    rm(hold_45, hold_85)
  }
  
  #clean up
  rm(rate_hold_45, rate_hold_85)
}

###############################################################
###TABLES AND GRAPHICS USING OVERALL MEAN VALUES NOT BY HERD###
###############################################################

#CREATE DATAFRAME OF MEAN AND CONFIDENCE INTERVALS OVERALL NOT SEGMENTING BY HERD

#DURATION
#calculate summary stats
metrics_dur <- dur %>% group_by(model) %>% 
  summarise(mean_dur = ci(duration, na.rm = T)[1], 
              low_ci_dur = ci(duration, na.rm = T)[2],
              hi_ci_dur = ci(duration, na.rm = T)[3], 
              sd_dur = ci (duration, na.rm = T)[4])

#separate model parameters into parts
metrics_dur[,c('model', 'rcp', 'period')] <- str_split_fixed(metrics_dur$model, '_', 3)
metrics_dur$model <- str_replace_all(metrics_dur$model, '[.]', '-')

#set rcp category for baseline
metrics_dur$rcp[metrics_dur$model == 'base'] <- 'Baseline'

#set rcp category for BIT (2000-2019)
metrics_dur$rcp[metrics_dur$period == '1999' & metrics_dur$rcp == 'rcp45'] <- 'BIT 4.5'
metrics_dur$rcp[metrics_dur$period == '1999' & metrics_dur$rcp == 'rcp85'] <- 'BIT 8.5'

#ORDER
#calculate summary stats
metrics_ord <- order %>% group_by(model) %>% 
  summarise(mean_ord = ci(order, na.rm = T)[1], 
            low_ci_ord = ci(order, na.rm = T)[2],
            hi_ci_ord = ci(order, na.rm = T)[3], 
            sd_ord = ci (order, na.rm = T)[4])

#separate model parameters into parts
metrics_ord[,c('model', 'rcp', 'period')] <- str_split_fixed(metrics_ord$model, '_', 3)
metrics_ord$model <- str_replace_all(metrics_ord$model, '[.]', '-')

#set rcp category for baseline
metrics_ord$rcp[metrics_ord$model == 'base'] <- 'Baseline'

#set rcp category for BIT (2000-2019)
metrics_ord$rcp[metrics_ord$period == '1999' & metrics_ord$rcp == 'rcp45'] <- 'BIT 4.5'
metrics_ord$rcp[metrics_ord$period == '1999' & metrics_ord$rcp == 'rcp85'] <- 'BIT 8.5'

#RATE
#calculate summary stats
metrics_rate <- rate %>% group_by(model) %>% 
  summarise(mean_rate = ci(rate, na.rm = T)[1], 
            low_ci_rate = ci(rate, na.rm = T)[2],
            hi_ci_rate = ci(rate, na.rm = T)[3], 
            sd_rate = ci (rate, na.rm = T)[4])

#separate model parameters into parts
metrics_rate[,c('model', 'rcp', 'period')] <- str_split_fixed(metrics_rate$model, '_', 3)
metrics_rate$model <- str_replace_all(metrics_rate$model, '[.]', '-')

#set rcp category for baseline
metrics_rate$rcp[metrics_rate$model == 'base'] <- 'Baseline'

#set rcp category for BIT (2000-2019)
metrics_rate$rcp[metrics_rate$period == '1999' & metrics_rate$rcp == 'rcp45'] <- 'BIT 4.5'
metrics_rate$rcp[metrics_rate$period == '1999' & metrics_rate$rcp == 'rcp85'] <- 'BIT 8.5'

#merge dfs
metrics <- merge(metrics_dur, metrics_ord, by = c("model", "rcp", "period"))
metrics <- merge(metrics, metrics_rate, by = c("model", "rcp", "period"))

#clean up metrics dfs
rm(metrics_dur, metrics_ord, metrics_rate)

#ORGANIZE DATA FOR PLOTTING

#set universal theme for out
theme_out <- theme(plot.title = element_text(size = 20, hjust = 0.5), axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                   axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                   legend.position = '')

#create df for mid century only and clean up
metrics_2040 <- metrics[metrics$period != 2070,]
metrics_2040$model[metrics_2040$model == 'base'] <- 'Baseline'
metrics_2040$rcp[metrics_2040$rcp == 'BIT 4.5'] <- 'rcp45'
metrics_2040$rcp[metrics_2040$rcp == 'BIT 8.5'] <- 'rcp85'

#separate by rcp
metrics_2040_rcp_45 <- metrics_2040[metrics_2040$rcp != 'rcp85',]
metrics_2040_rcp_45$model <- ordered(metrics_2040_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
metrics_2040_rcp_85 <- metrics_2040[metrics_2040$rcp != 'rcp45',]
metrics_2040_rcp_85$model <- ordered(metrics_2040_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
#set group
metrics_2040_rcp_45$group <- 'AMid-Century \nRCP 4.5'
metrics_2040_rcp_85$group <- 'BMid-Century \nRCP 8.5'

#create df for end of century only and clean up
metrics_2070 <- metrics[metrics$period != 2040,]
metrics_2070$model[metrics_2070$model == 'base'] <- 'Baseline'
metrics_2070$rcp[metrics_2070$rcp == 'BIT 4.5'] <- 'rcp45'
metrics_2070$rcp[metrics_2070$rcp == 'BIT 8.5'] <- 'rcp85'

#separate by rcp
metrics_2070_rcp_45 <- metrics_2070[metrics_2070$rcp != 'rcp85',]
metrics_2070_rcp_45$model <- ordered(metrics_2070_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
metrics_2070_rcp_85 <- metrics_2070[metrics_2070$rcp != 'rcp45',]
metrics_2070_rcp_85$model <- ordered(metrics_2070_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
#set group
metrics_2070_rcp_45$group <- 'CEnd of Century \nRCP 4.5'
metrics_2070_rcp_85$group <- 'DEnd of Century \nRCP 8.5'

#create baseline dataframes using average data
base_out_mean <- data.frame(group = c('AMid-Century \nRCP 4.5',
                                     'DEnd of Century \nRCP 8.5'),
                           dur = metrics_2070_rcp_45$mean_dur[metrics_2070_rcp_45$model == 'Baseline'],
                           ord = metrics_2070_rcp_45$mean_ord[metrics_2070_rcp_45$model == 'Baseline'],
                           ss = metrics_2070_rcp_45$mean_rate[metrics_2070_rcp_45$model == 'Baseline'])

base_in_mean <- data.frame(group = c('BMid-Century \nRCP 8.5',
                                    'CEnd of Century \nRCP 4.5'),
                          dur = metrics_2070_rcp_45$mean_dur[metrics_2070_rcp_45$model == 'Baseline'],
                          ord = metrics_2070_rcp_45$mean_ord[metrics_2070_rcp_45$model == 'Baseline'],
                          ss = metrics_2070_rcp_45$mean_rate[metrics_2070_rcp_45$model == 'Baseline'])

base_out_loci <- data.frame(group = c('AMid-Century \nRCP 4.5',
                                      'DEnd of Century \nRCP 8.5'),
                            dur = metrics_2070_rcp_45$low_ci_dur[metrics_2070_rcp_45$model == 'Baseline'],
                            ord = metrics_2070_rcp_45$low_ci_ord[metrics_2070_rcp_45$model == 'Baseline'],
                            ss = metrics_2070_rcp_45$low_ci_rate[metrics_2070_rcp_45$model == 'Baseline'])

base_in_loci <- data.frame(group = c('BMid-Century \nRCP 8.5',
                                     'CEnd of Century \nRCP 4.5'),
                            dur = metrics_2070_rcp_45$low_ci_dur[metrics_2070_rcp_45$model == 'Baseline'],
                            ord = metrics_2070_rcp_45$low_ci_ord[metrics_2070_rcp_45$model == 'Baseline'],
                            ss = metrics_2070_rcp_45$low_ci_rate[metrics_2070_rcp_45$model == 'Baseline'])

base_out_hici <- data.frame(group = c('AMid-Century \nRCP 4.5',
                                      'DEnd of Century \nRCP 8.5'),
                            dur = metrics_2070_rcp_45$hi_ci_dur[metrics_2070_rcp_45$model == 'Baseline'],
                            ord = metrics_2070_rcp_45$hi_ci_ord[metrics_2070_rcp_45$model == 'Baseline'],
                            ss = metrics_2070_rcp_45$hi_ci_rate[metrics_2070_rcp_45$model == 'Baseline'])

base_in_hici <- data.frame(group = c('BMid-Century \nRCP 8.5',
                                     'CEnd of Century \nRCP 4.5'),
                            dur = metrics_2070_rcp_45$hi_ci_dur[metrics_2070_rcp_45$model == 'Baseline'],
                            ord = metrics_2070_rcp_45$hi_ci_ord[metrics_2070_rcp_45$model == 'Baseline'],
                            ss = metrics_2070_rcp_45$hi_ci_rate[metrics_2070_rcp_45$model == 'Baseline'])

#set col scale
myColors_avg <- c('#4477AA', '#000000', '#228833', '#CCBB44', '#EE6677', 
                  '#AA3377')
names(myColors_avg) <- c('CNRM-CM5', 'Baseline',
                         'HadGEM2-ES365', 'inmcm4', 'IPSL-CM5A-MR',
                         'MIROC-ESM-CHEM')
colScale_avg <- scale_fill_manual(name = "Model (BIT Transparent)", values = myColors_avg,
                                  breaks = c('Baseline', 'HadGEM2-ES365', 'inmcm4', 'CNRM-CM5',
                                             'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'),
                                  labels = c('Baseline' = 'Baseline',
                                             'HadGEM2-ES365' = 'Middle',
                                             'inmcm4' = 'Low Temp/Low Precip',
                                             'CNRM-CM5' = 'Low Temp/High Precip', 
                                             'IPSL-CM5A-MR' = 'High Temp/Low Precip',
                                             'MIROC-ESM-CHEM' = 'High Temp/High Precip'))

#set linetype aesthetic
line_types <- c("Mean" = 'dashed',"CI" = 'dotted')

#set universal theme for out
theme_out <- theme(plot.title = element_text(size = 20, hjust = 0.5), axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                   axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                   legend.title = element_text(size = 15),
                   legend.position = '',
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#set alpha aesthetic
alpha_types <- c('1999' = 0.4, '2040' = 1, '2070' = 1)

#put all organized metrics back together
metrics_plot <- rbind(metrics_2040_rcp_45, metrics_2040_rcp_85,
                      metrics_2070_rcp_45, metrics_2070_rcp_85)

#plot bar charts separately
#duration
ggplot(metrics_plot[metrics_plot$model != 'Baseline',], aes(fill = model, x = group, y = mean_dur, alpha = factor(period))) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = low_ci_dur, ymax = hi_ci_dur),
                width = 0.4, position=position_dodge(.9)) +
  xlab('Metric and Emissions Scenario') + 
  theme_bw() + theme_out + colScale_avg +
  scale_y_continuous(name = 'Days', limits = c(0, 20)) +
  scale_x_discrete(labels = c('Mid-Century \nRCP 4.5', 'Mid-Century \nRCP 8.5', 
                              'End of Century \nRCP 4.5', 'End of Century \nRCP 8.5')) +
  ggtitle('') +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                           '2040' = 'Projection')) +
  geom_hline(yintercept = base_out_mean$dur[1], linetype = 'dashed', size = 0.75) +
  geom_hline(yintercept = base_out_loci$dur[1], linetype = 'dotted', size = 0.75) +
  geom_hline(yintercept = base_out_hici$dur[1], linetype = 'dotted', size = 0.75)

ggsave('output/final_plots/duration_bar_ci.png',
       width = 9, height = 3, unit = "in")

#order
ggplot(metrics_plot[metrics_plot$model != 'Baseline',], aes(fill = model, x = group, y = mean_ord, alpha = factor(period))) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = low_ci_ord, ymax = hi_ci_ord),
                width = 0.4, position=position_dodge(.9)) +
  colScale_avg + xlab('Metric and Emissions Scenario') + 
  theme_bw() + theme_out +
  scale_y_continuous(name = 'Correlation', limits = c(0, 0.75)) +
  scale_x_discrete(labels = c('Mid-Century \nRCP 4.5', 'Mid-Century \nRCP 8.5', 
                              'End of Century \nRCP 4.5', 'End of Century \nRCP 8.5')) +
  ggtitle('') +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                           '2040' = 'Projection')) +
  geom_hline(yintercept = base_out_mean$ord[1], linetype = 'dashed', size = 0.75) +
  geom_hline(yintercept = base_out_loci$ord[1], linetype = 'dotted', size = 0.75) +
  geom_hline(yintercept = base_out_hici$ord[1], linetype = 'dotted', size = 0.75)

ggsave('output/final_plots/order_bar_ci.png',
       width = 9, height = 3, unit = "in")

#mean spring scale
ggplot(metrics_plot[metrics_plot$model != 'Baseline',], aes(fill = model, x = group, y = mean_rate, alpha = factor(period))) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = low_ci_rate, ymax = hi_ci_rate),
                width = 0.4, position=position_dodge(.9)) +
  colScale_avg + xlab('Emissions Scenario') + 
  theme_bw() + theme_out +
  scale_y_continuous(name = 'Days', limits = c(0,20)) +
  scale_x_discrete(labels = c('Mid-Century \nRCP 4.5', 'Mid-Century \nRCP 8.5', 
                              'End of Century \nRCP 4.5', 'End of Century \nRCP 8.5')) +
  ggtitle('') +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                           '2040' = 'Projection')) +
  geom_hline(yintercept = base_out_mean$ss[1], linetype = 'dashed', size = 0.75) +
  geom_hline(yintercept = base_out_loci$ss[1], linetype = 'dotted', size = 0.75) +
  geom_hline(yintercept = base_out_hici$ss[1], linetype = 'dotted', size = 0.75)

ggsave('output/final_plots/mean_spring_scale_bar_ci.png',
       width = 9, height = 3, unit = "in")

#################################################
###TABLES AND GRAPHICS FOR MID AND END CENTURY###
#################################################

#create list of herds so can separate values
herds <- c('AtlanticRimNorth',
           'AtlanticRimSouth',
           'Chokecherry',
           'Cody',
           'DEERP',
           'Dubois',
           'Lander',
           'Meeteetse',
           'Pinedale',
           'PlatteValleyNorth',
           'PlatteValleySouth',
           'Superior',
           'SweetWaterGreenMtn',
           'Teton',
           'WyomingRange')

#create df of mean values across herds and scenarios
mig_dat <- data.frame(herd = character(), model = character(), duration = numeric(), order = numeric(), rate = numeric())

#fill df with values from herds
for(i in herds){
  
  #subset and calc the mean dur
  dur_h <- dur[str_detect(dur$global_id, i) == T,] %>% group_by(model) %>% 
    mutate(duration = mean(duration, na.rm = T)) %>% na.omit
  
  #subset and calc the mean order
  order_h <- order[str_detect(order$global_id, i) == T,] %>% group_by(model) %>% 
    mutate(order = mean(order, na.rm = T)) %>% na.omit
  
  #subset and calc the mean rate
  rate_h <- rate[str_detect(rate$global_id, i) == T,] %>% group_by(model) %>% 
    mutate(rate = mean(rate, na.rm = T)) %>% na.omit
  
  #take one row for each unique value
  for(j in unique(dur_h$model)){
    
    #rbind to main df
    mig_dat <- rbind(mig_dat, data.frame(herd = i,
                                         model = j,
                                         duration = dur_h$duration[dur_h$model == j][1],
                                         order = order_h$order[order_h$model == j][1],
                                         rate = rate_h$rate[rate_h$model == j][1]))
  }
}

rm(dur_h, order_h, rate_h)

#separate model parameters into parts
mig_dat[,c('model', 'rcp', 'period')] <- str_split_fixed(mig_dat$model, '_', 3)
mig_dat$model <- str_replace_all(mig_dat$model, '[.]', '-')

#set rcp category for baseline
mig_dat$rcp[mig_dat$model == 'base'] <- 'Baseline'

#set rcp category for BIT (2000-2019)
mig_dat$rcp[mig_dat$period == '1999' & mig_dat$rcp == 'rcp45'] <- 'BIT 4.5'
mig_dat$rcp[mig_dat$period == '1999' & mig_dat$rcp == 'rcp85'] <- 'BIT 8.5'

#create mig_dat df to export as csv table and clean up
mig_out <- mig_dat[, c('herd', 'model', 'rcp', 'period', 'duration', 'order', 'rate')]
mig_out$duration <- round(mig_out$duration, 2)
mig_out$order <- round(mig_out$order, 2)
mig_out$rate <- round(mig_out$rate, 2)
mig_out$model[mig_out$model == 'base'] <- 'Baseline'
mig_out$period[mig_out$period == ''] <- '2000-2019'
mig_out$period[mig_out$period == '1999'] <- '2000-2019'
mig_out$period[mig_out$period == '2040'] <- '2040-2069'
mig_out$period[mig_out$period == '2070'] <- '2070-2099'
mig_out$rcp[mig_out$rcp == 'rcp45'] <- 'RCP 4.5'
mig_out$rcp[mig_out$rcp == 'rcp85'] <- 'RCP 8.5'

#export table
write.csv(mig_out, file = 'output/final_plots/tables/greenscape_metrics_by_herd.csv', row.names = F)

#set color scale for all 6 models
#need to choose a color palette colorblind friendly
#https://personal.sron.nl/~pault/#sec:qualitative
myColors <- c('#BBBBBB', '#4477AA', '#228833', '#CCBB44', '#EE6677', 
              '#AA3377')
names(myColors) <- unique(mig_dat$model)
colScale <- scale_color_manual(name = "Model", values = myColors, labels = c('base' = 'Baseline'))
fillScale <- scale_fill_manual(name = "Model", values = myColors)

#set universal theme
theme_e <- theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                 axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                 legend.title = element_text(size = 15))

#load source function for 3d grids
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

# #loop through each herd to plot individual herd outputs
# for(herd in herds){
#   
#   #plot order vs. duration mid century
#   ggplot(mig_dat[mig_dat$period != 2070 & mig_dat$herd == herd,], 
#          aes(x = order, y = duration, col = model, shape = rcp)) +
#     geom_point(size = 4) + theme_bw() + colScale + theme_e +
#     scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#     xlab('Order') + ylab('Duration') + ggtitle(str_c(herd, ' Mid-Century (2040-2069)')) + 
#     scale_x_continuous(limits = c(-0.8, 1)) + scale_y_continuous(limits = c(-5, 31))
#   
#   ggsave(filename = str_c('output/final_plots/by_herd/avg_order_duration_mid_century_', herd,'.png'))
#   
#   #plot order vs. duration end of century
#   ggplot(mig_dat[mig_dat$period != 2040 & mig_dat$herd == herd,], 
#          aes(x = order, y = duration, col = model, shape = rcp)) +
#     geom_point(size = 4) + theme_bw() + colScale + theme_e +
#     scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#     xlab('Order') + ylab('Duration') + ggtitle(str_c(herd, ' End of Century (2070-2099)')) +
#     scale_x_continuous(limits = c(-0.8, 1)) + scale_y_continuous(limits = c(-5, 31))
#   
#   ggsave(filename = str_c('output/final_plots/by_herd/avg_order_duration_end_century_', herd,'.png'))
#   
#   #plot order vs. rate mid century
#   ggplot(mig_dat[mig_dat$period != 2070 & mig_dat$herd == herd,], 
#          aes(x = order, y = rate, col = model, shape = rcp)) +
#     geom_point(size = 4) + theme_bw() + colScale + theme_e +
#     scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#     xlab('Order') + ylab('Rate') + ggtitle(str_c(herd, ' Mid-Century (2040-2069)')) + 
#     scale_x_continuous(limits = c(-0.8, 1)) + scale_y_continuous(limits = c(10, 22))
#   
#   ggsave(filename = str_c('output/final_plots/by_herd/avg_order_rate_mid_century_', herd,'.png'))
#   
#   #plot order vs. rate end of century
#   ggplot(mig_dat[mig_dat$period != 2040 & mig_dat$herd == herd,], 
#          aes(x = order, y = rate, col = model, shape = rcp)) +
#     geom_point(size = 4) + theme_bw() + colScale + theme_e +
#     scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#     xlab('Order') + ylab('Rate') + ggtitle(str_c(herd, ' End of Century (2070-2099)')) +
#     scale_x_continuous(limits = c(-0.8, 1)) + scale_y_continuous(limits = c(10, 22))
#   
#   ggsave(filename = str_c('output/final_plots/by_herd/avg_order_rate_end_century_', herd,'.png'))
#   
#   #plot rate vs. duration mid century
#   ggplot(mig_dat[mig_dat$period != 2070 & mig_dat$herd == herd,], 
#          aes(x = rate, y = duration, col = model, shape = rcp)) +
#     geom_point(size = 4) + theme_bw() + colScale + theme_e +
#     scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#     xlab('Rate') + ylab('Duration') + ggtitle(str_c(herd, ' Mid-Century (2040-2069)')) + 
#     scale_x_continuous(limits = c(10, 22)) + scale_y_continuous(limits = c(-5, 31))
#   
#   ggsave(filename = str_c('output/final_plots/by_herd/avg_rate_duration_mid_century_', herd,'.png'))
#   
#   #plot rate vs. duration end of century
#   ggplot(mig_dat[mig_dat$period != 2040 & mig_dat$herd == herd,], 
#          aes(x = rate, y = duration, col = model, shape = rcp)) +
#     geom_point(size = 4) + theme_bw() + colScale + theme_e +
#     scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#     xlab('Rate') + ylab('Duration') + ggtitle(str_c(herd, ' End of Century (2070-2099)')) +
#     scale_x_continuous(limits = c(10, 22)) + scale_y_continuous(limits = c(-5, 31))
#   
#   ggsave(filename = str_c('output/final_plots/by_herd/avg_rate_duration_end_century_', herd,'.png'))
#   
# }

# #create plots by model showing various herds
# #plot order vs. duration mid century
# ggplot(mig_dat[mig_dat$period != 2070 & mig_dat$model == 'HadGEM2-ES365',], 
#        aes(x = order, y = duration, col = herd, shape = rcp)) +
#   geom_point(size = 4) + theme_bw() + 
# #  colScale + 
#   theme_e +
#   scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   xlab('Order') + ylab('Duration') + ggtitle('HadGEM2-ES365 Mid-Century (2040-2069)') + 
#   scale_x_continuous(limits = c(-0.8, 1)) + scale_y_continuous(limits = c(-5, 31))

#create mig data average across herds
mig_dat_avg <- mig_dat

#create unique scenarios
mig_dat_avg$scenario <- str_c(mig_dat_avg$model, '_', mig_dat_avg$rcp, '_', mig_dat_avg$period)

#take average order, duration, and rate
mig_dat_avg <- mig_dat_avg %>% group_by(scenario) %>%
  mutate(duration_avg = mean(duration, na.rm = T)) %>% 
  mutate(order_avg = mean(order, na.rm = T)) %>%
  mutate(rate_avg = mean(rate, na.rm = T))

#remove duplicates
mig_dat_avg <- mig_dat_avg[!duplicated(mig_dat_avg$scenario),]

#remove unnecessary columns to limit confusion
mig_dat_avg <- mig_dat_avg[!(colnames(mig_dat_avg) %in% c('herd', 'duration', 'order', 'rate'))]

#create mig_dat_avg table to export
mig_out_avg <- mig_dat_avg[, c('model', 'rcp', 'period', 'duration_avg', 'order_avg', 'rate_avg')]
colnames(mig_out_avg) <- c('model', 'rcp', 'period', 'duration', 'order', 'rate')
mig_out_avg$duration <- round(mig_out_avg$duration, 2)
mig_out_avg$order <- round(mig_out_avg$order, 2)
mig_out_avg$rate <- round(mig_out_avg$rate, 2)
mig_out_avg$model[mig_out_avg$model == 'base'] <- 'Baseline'
mig_out_avg$period[mig_out_avg$period == ''] <- '2000-2019'
mig_out_avg$period[mig_out_avg$period == '1999'] <- '2000-2019'
mig_out_avg$period[mig_out_avg$period == '2040'] <- '2040-2069'
mig_out_avg$period[mig_out_avg$period == '2070'] <- '2070-2099'
mig_out_avg$rcp[mig_out_avg$rcp == 'rcp45'] <- 'RCP 4.5'
mig_out_avg$rcp[mig_out_avg$rcp == 'rcp85'] <- 'RCP 8.5'

#export table
write.csv(mig_out_avg, file = 'output/final_plots/tables/greenscape_metrics_avg.csv', row.names = F)

#set universal theme for out
theme_out <- theme(plot.title = element_text(size = 20, hjust = 0.5), axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                   axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                   legend.position = '')

# #plot avg order vs. avg duration mid century
# ggplot(mig_dat_avg[mig_dat_avg$period != 2070,], 
#        aes(x = order_avg, y = duration_avg, col = model, shape = rcp)) +
#   geom_point(size = 4) + theme_bw() + colScale + theme_out +
#   scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   xlab('Order (Cor)') + ylab('Duration (Days)') + ggtitle('Mid-Century (2040-2069)') +
#   scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 25))
# 
# ggsave(filename = 'output/final_plots/avg_order_duration_mid_century.png',
#        width = 5, height = 5, unit = "in")
# 
# #plot avg order vs. avg duration end of century
# ggplot(mig_dat_avg[mig_dat_avg$period != 2040,], 
#        aes(x = order_avg, y = duration_avg, col = model, shape = rcp)) +
#   geom_point(size = 4) + theme_bw() + colScale + theme_out +
#   scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   xlab('Order (Cor)') + ylab('Duration (Days)') + ggtitle('End of Century (2070-2099)') +
#   scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 25))
# 
# ggsave(filename = 'output/final_plots/avg_order_duration_end_century.png',
#        width = 5, height = 5, unit = "in")
# 
# #plot avg order vs. avg rate mid century
# ggplot(mig_dat_avg[mig_dat_avg$period != 2070,], 
#        aes(x = order_avg, y = rate_avg, col = model, shape = rcp)) +
#   geom_point(size = 4) + theme_bw() + colScale + theme_out +
#   scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   xlab('Order (Cor)') + ylab('Mean Spring Scale (Days)') + ggtitle('Mid-Century (2040-2069)') +
#   scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(10, 20))
# 
# ggsave(filename = 'output/final_plots/avg_order_rate_mid_century.png',
#        width = 5, height = 5, unit = "in")
# 
# #plot avg order vs. avg rate end of century
# ggplot(mig_dat_avg[mig_dat_avg$period != 2040,], 
#        aes(x = order_avg, y = rate_avg, col = model, shape = rcp)) +
#   geom_point(size = 4) + theme_bw() + colScale + theme_out +
#   scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   xlab('Order (Cor)') + ylab('Mean Spring Scale (Days)') + ggtitle('End of Century (2070-2099)') +
#   scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(10, 20))
# 
# ggsave(filename = 'output/final_plots/avg_order_rate_end_century.png',
#        width = 5, height = 5, unit = "in")
# 
# #plot avg rate vs. avg duration mid century
# ggplot(mig_dat_avg[mig_dat_avg$period != 2070,], 
#        aes(x = rate_avg, y = duration_avg, col = model, shape = rcp)) +
#   geom_point(size = 4) + theme_bw() + colScale + theme_out +
#   scale_shape_manual(values = c(16, 17, 15, 2, 12), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   xlab('Mean Spring Scale (Days)') + ylab('Duration (Days)') + ggtitle('Mid-Century (2040-2069)') +
#   scale_x_continuous(limits = c(10, 20)) + scale_y_continuous(limits = c(0, 25))
# 
# ggsave(filename = 'output/final_plots/avg_rate_duration_mid_century.png',
#        width = 5, height = 5, unit = "in")
# 
# #plot avg rate vs. avg duration end of century
# ggplot(mig_dat_avg[mig_dat_avg$period != 2040,], 
#        aes(x = rate_avg, y = duration_avg, col = model, shape = rcp)) +
#   geom_point(size = 4) + theme_bw() + colScale + theme_out +
#   scale_shape_manual(values = c(16, 15, 17, 0, 2), labels = c('rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   xlab('Mean Spring Scale (Days)') + ylab('Duration (Days)') + ggtitle('End of Century (2070-2099)') +
#   scale_x_continuous(limits = c(10, 20)) + scale_y_continuous(limits = c(0, 25))
# 
# ggsave(filename = 'output/final_plots/avg_rate_duration_end_century.png',
#        width = 5, height = 5, unit = "in")

#create df for mid century only and clean up
mig_dat_avg_2040 <- mig_dat_avg[mig_dat_avg$period != 2070,]
mig_dat_avg_2040$model[mig_dat_avg_2040$model == 'base'] <- 'Baseline'
#mig_dat_avg_2040$model[mig_dat_avg_2040$model == 'HadGEM2-ES365' & mig_dat_avg_2040$period == 1999] <- 'Back in Time (HadGEM2)'
mig_dat_avg_2040$rcp[mig_dat_avg_2040$rcp == 'BIT 4.5'] <- 'rcp45'
mig_dat_avg_2040$rcp[mig_dat_avg_2040$rcp == 'BIT 8.5'] <- 'rcp85'

#duration
mig_dat_avg_2040_dur_rcp_45 <- mig_dat_avg_2040[mig_dat_avg_2040$rcp != 'rcp85',]
mig_dat_avg_2040_dur_rcp_45$model <- ordered(mig_dat_avg_2040_dur_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                        'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
mig_dat_avg_2040_dur_rcp_85 <- mig_dat_avg_2040[mig_dat_avg_2040$rcp != 'rcp45',]
mig_dat_avg_2040_dur_rcp_85$model <- ordered(mig_dat_avg_2040_dur_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))

#spring scale
mig_dat_avg_2040_ss_rcp_45 <- mig_dat_avg_2040[mig_dat_avg_2040$rcp != 'rcp85',]
mig_dat_avg_2040_ss_rcp_45$model <- ordered(mig_dat_avg_2040_ss_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
mig_dat_avg_2040_ss_rcp_85 <- mig_dat_avg_2040[mig_dat_avg_2040$rcp != 'rcp45',]
mig_dat_avg_2040_ss_rcp_85$model <- ordered(mig_dat_avg_2040_ss_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                         'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))

#order
mig_dat_avg_2040_ord_rcp_45 <- mig_dat_avg_2040[mig_dat_avg_2040$rcp != 'rcp85',]
mig_dat_avg_2040_ord_rcp_45$model <- ordered(mig_dat_avg_2040_ord_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
mig_dat_avg_2040_ord_rcp_85 <- mig_dat_avg_2040[mig_dat_avg_2040$rcp != 'rcp45',]
mig_dat_avg_2040_ord_rcp_85$model <- ordered(mig_dat_avg_2040_ord_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))

#set group
mig_dat_avg_2040_dur_rcp_45$group <- 'AMid-Century \nRCP 4.5'
mig_dat_avg_2040_dur_rcp_85$group <- 'BMid-Century \nRCP 8.5'
mig_dat_avg_2040_ss_rcp_45$group <- 'AMid-Century \nRCP 4.5'
mig_dat_avg_2040_ss_rcp_85$group <- 'BMid-Century \nRCP 8.5'
mig_dat_avg_2040_ord_rcp_45$group <- 'AMid-Century \nRCP 4.5'
mig_dat_avg_2040_ord_rcp_85$group <- 'BMid-Century \nRCP 8.5'

#create df for end of century only and clean up
mig_dat_avg_2070 <- mig_dat_avg[mig_dat_avg$period != 2040,]
mig_dat_avg_2070$model[mig_dat_avg_2070$model == 'base'] <- 'Baseline'
#mig_dat_avg_2070$model[mig_dat_avg_2070$model == 'HadGEM2-ES365' & mig_dat_avg_2070$period == 1999] <- 'Back in Time (HadGEM2)'
mig_dat_avg_2070$rcp[mig_dat_avg_2070$rcp == 'BIT 4.5'] <- 'rcp45'
mig_dat_avg_2070$rcp[mig_dat_avg_2070$rcp == 'BIT 8.5'] <- 'rcp85'

#duration
mig_dat_avg_2070_dur_rcp_45 <- mig_dat_avg_2070[mig_dat_avg_2070$rcp != 'rcp85',]
mig_dat_avg_2070_dur_rcp_45$model <- ordered(mig_dat_avg_2070_dur_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
mig_dat_avg_2070_dur_rcp_85 <- mig_dat_avg_2070[mig_dat_avg_2070$rcp != 'rcp45',]
mig_dat_avg_2070_dur_rcp_85$model <- ordered(mig_dat_avg_2070_dur_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))

#spring scale
mig_dat_avg_2070_ss_rcp_45 <- mig_dat_avg_2070[mig_dat_avg_2070$rcp != 'rcp85',]
mig_dat_avg_2070_ss_rcp_45$model <- ordered(mig_dat_avg_2070_ss_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
mig_dat_avg_2070_ss_rcp_85 <- mig_dat_avg_2070[mig_dat_avg_2070$rcp != 'rcp45',]
mig_dat_avg_2070_ss_rcp_85$model <- ordered(mig_dat_avg_2070_ss_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))

#order
mig_dat_avg_2070_ord_rcp_45 <- mig_dat_avg_2070[mig_dat_avg_2070$rcp != 'rcp85',]
mig_dat_avg_2070_ord_rcp_45$model <- ordered(mig_dat_avg_2070_ord_rcp_45$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))
mig_dat_avg_2070_ord_rcp_85 <- mig_dat_avg_2070[mig_dat_avg_2070$rcp != 'rcp45',]
mig_dat_avg_2070_ord_rcp_85$model <- ordered(mig_dat_avg_2070_ord_rcp_85$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'Baseline',
                                                                                           'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))

#set group
mig_dat_avg_2070_dur_rcp_45$group <- 'CEnd of Century \nRCP 4.5'
mig_dat_avg_2070_dur_rcp_85$group <- 'DEnd of Century \nRCP 8.5'
mig_dat_avg_2070_ss_rcp_45$group <- 'CEnd of Century \nRCP 4.5'
mig_dat_avg_2070_ss_rcp_85$group <- 'DEnd of Century \nRCP 8.5'
mig_dat_avg_2070_ord_rcp_45$group <- 'CEnd of Century \nRCP 4.5'
mig_dat_avg_2070_ord_rcp_85$group <- 'DEnd of Century \nRCP 8.5'

#create baseline dataframes using migration group avg
baseline_out <- data.frame(group = c('AMid-Century \nRCP 4.5',
                                     'DEnd of Century \nRCP 8.5'),
                           dur = mig_dat_avg_2070_dur_rcp_45$duration_avg[mig_dat_avg_2070_dur_rcp_45$model == 'Baseline'],
                           ord = mig_dat_avg_2070_ord_rcp_45$order_avg[mig_dat_avg_2070_ord_rcp_45$model == 'Baseline'],
                           ss = mig_dat_avg_2070_ss_rcp_45$rate_avg[mig_dat_avg_2070_ss_rcp_45$model == 'Baseline'])
baseline_in <- data.frame(group = c('BMid-Century \nRCP 8.5',
                                     'CEnd of Century \nRCP 4.5'),
                           dur = mig_dat_avg_2070_dur_rcp_45$duration_avg[mig_dat_avg_2070_dur_rcp_45$model == 'Baseline'],
                           ord = mig_dat_avg_2070_ord_rcp_45$order_avg[mig_dat_avg_2070_ord_rcp_45$model == 'Baseline'],
                           ss = mig_dat_avg_2070_ss_rcp_45$rate_avg[mig_dat_avg_2070_ss_rcp_45$model == 'Baseline'])

#set col scale
myColors_avg <- c('#4477AA', '#000000', '#228833', '#CCBB44', '#EE6677', 
              '#AA3377')
names(myColors_avg) <- c('CNRM-CM5', 'Baseline',
                         'HadGEM2-ES365', 'inmcm4', 'IPSL-CM5A-MR',
                         'MIROC-ESM-CHEM')
colScale_avg <- scale_fill_manual(name = "Model (BIT Transparent)", values = myColors_avg,
                                  breaks = c('Baseline', 'HadGEM2-ES365', 'inmcm4', 'CNRM-CM5',
                                             'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'),
                                  labels = c('Baseline' = 'Baseline',
                                             'HadGEM2-ES365' = 'Middle',
                                             'inmcm4' = 'Low Temp/Low Precip',
                                             'CNRM-CM5' = 'Low Temp/High Precip', 
                                             'IPSL-CM5A-MR' = 'High Temp/Low Precip',
                                             'MIROC-ESM-CHEM' = 'High Temp/High Precip'))

#set line type
line_types <- c('Baseline' = 'dashed')

#set universal theme for out
theme_out <- theme(plot.title = element_text(size = 20, hjust = 0.5), axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                   axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                   legend.title = element_text(size = 15),
                   legend.position = '')

#set alpha aesthetic
alpha_types <- c('1999' = 0.4, '2040' = 1, '2070' = 1)

#plot bar charts separately
#duration
ggplot() +
  geom_bar(data = mig_dat_avg_2040_dur_rcp_45[mig_dat_avg_2040_dur_rcp_45$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2040_dur_rcp_85[mig_dat_avg_2040_dur_rcp_85$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2070_dur_rcp_45[mig_dat_avg_2040_dur_rcp_45$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2070_dur_rcp_85[mig_dat_avg_2040_dur_rcp_85$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg, alpha = factor(period))) +
  xlab('Metric and Emissions Scenario') + 
  theme_bw() + theme_out + colScale_avg +
  scale_y_continuous(name = 'Days', limits = c(0, 20)) +
  scale_x_discrete(labels = c('Mid-Century \nRCP 4.5', 'Mid-Century \nRCP 8.5', 
                              'End of Century \nRCP 4.5', 'End of Century \nRCP 8.5')) +
  geom_errorbar(data = baseline_out, width = 1.16, aes(x=group, ymax = dur, ymin = dur, linetype = 'Baseline'), 
                colour="#000000") +
  geom_errorbar(data = baseline_in, width = 1, aes(x=group, ymax = dur, ymin = dur, linetype = 'Baseline'), 
                colour="#000000") +
  scale_linetype_manual(name = '', values = line_types) + 
  ggtitle('') +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                           '2040' = 'Projection'))

ggsave('output/final_plots/duration_bar_legend.png',
       width = 9, height = 3, unit = "in")

#Order
ggplot() +
  geom_bar(data = mig_dat_avg_2040_ord_rcp_45[mig_dat_avg_2040_ord_rcp_45$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2040_ord_rcp_85[mig_dat_avg_2040_ord_rcp_85$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2070_ord_rcp_45[mig_dat_avg_2040_ord_rcp_45$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2070_ord_rcp_85[mig_dat_avg_2040_ord_rcp_85$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg, alpha = factor(period))) +
  colScale_avg + xlab('Metric and Emissions Scenario') + 
  theme_bw() + theme_out +
  scale_y_continuous(name = 'Correlation', limits = c(0, 0.75)) +
  scale_x_discrete(labels = c('Mid-Century \nRCP 4.5', 'Mid-Century \nRCP 8.5', 
                              'End of Century \nRCP 4.5', 'End of Century \nRCP 8.5')) +
  geom_errorbar(data = baseline_out, width = 1.16, aes(x=group, ymax = ord, ymin = ord, linetype = 'Baseline'), 
                colour="#000000") +
  geom_errorbar(data = baseline_in, width = 1, aes(x=group, ymax = ord, ymin = ord, linetype = 'Baseline'), 
                colour="#000000") +
  scale_linetype_manual(name = '', values = line_types) +
  ggtitle('') +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                           '2040' = 'Projection'))

ggsave('output/final_plots/order_bar.png',
       width = 9, height = 3, unit = "in")

#mean spring scale
ggplot() +
  geom_bar(data = mig_dat_avg_2040_ss_rcp_45[mig_dat_avg_2040_ss_rcp_45$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2040_ss_rcp_85[mig_dat_avg_2040_ss_rcp_85$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2070_ss_rcp_45[mig_dat_avg_2040_ss_rcp_45$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg, alpha = factor(period))) +
  geom_bar(data = mig_dat_avg_2070_ss_rcp_85[mig_dat_avg_2040_ss_rcp_85$model != 'Baseline',], 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg, alpha = factor(period))) +
  colScale_avg + xlab('Emissions Scenario') + 
  theme_bw() + theme_out +
  scale_y_continuous(name = 'Days', limits = c(0,20)) +
  scale_x_discrete(labels = c('Mid-Century \nRCP 4.5', 'Mid-Century \nRCP 8.5', 
                              'End of Century \nRCP 4.5', 'End of Century \nRCP 8.5')) +
  geom_errorbar(data = baseline_out, width = 1.16, aes(x=group, ymax = ss, ymin = ss, linetype = 'Baseline'), 
                colour="#000000", linetype = 'dashed') +
  geom_errorbar(data = baseline_in, width = 1, aes(x=group, ymax = ss, ymin = ss, linetype = 'Baseline'), 
                colour="#000000") +
  scale_linetype_manual(name = '', values = line_types) +
  ggtitle('') +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                           '2040' = 'Projection'))

ggsave('output/final_plots/mean_spring_scale_bar.png',
       width = 9, height = 3, unit = "in")

#plot three metrics mid-century
ggplot() +
  geom_bar(data = mig_dat_avg_2040_dur_rcp_45, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg)) +
  geom_bar(data = mig_dat_avg_2040_dur_rcp_85, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg)) +
  geom_bar(data = mig_dat_avg_2040_ss_rcp_45, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg)) +
  geom_bar(data = mig_dat_avg_2040_ss_rcp_85, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg)) +
  geom_bar(data = mig_dat_avg_2040_ord_rcp_45, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg*20)) +
  geom_bar(data = mig_dat_avg_2040_ord_rcp_85, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg*20)) +
  colScale_avg + xlab('Metric and Emissions Scenario') + 
  theme_bw() + theme_out +
  scale_y_continuous(name = 'Days', limits = c(0, 20), sec.axis = sec_axis(~./20, name = 'Correlation')) +
  ggtitle('Mid-Century (2040-2069)')

ggsave('output/final_plots/greenscape_2040.png',
       width = 9, height = 3, unit = "in")

#plot three metrics end of century
ggplot() +
  geom_bar(data = mig_dat_avg_2070_dur_rcp_45, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg)) +
  geom_bar(data = mig_dat_avg_2070_dur_rcp_85, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = duration_avg)) +
  geom_bar(data = mig_dat_avg_2070_ss_rcp_45, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg)) +
  geom_bar(data = mig_dat_avg_2070_ss_rcp_85, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = rate_avg)) +
  geom_bar(data = mig_dat_avg_2070_ord_rcp_45, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg*20)) +
  geom_bar(data = mig_dat_avg_2070_ord_rcp_85, 
           stat = 'identity', position = position_dodge(), 
           aes(fill = model, x = group, y = order_avg*20)) +
  colScale_avg + xlab('Metric and Emissions Scenario') + 
  theme_bw() + theme_out +
  scale_y_continuous(name = 'Days', limits = c(0, 20), sec.axis = sec_axis(~./20, name = 'Correlation')) +
  ggtitle('End of Century (2070-2099)')

ggsave('output/final_plots/greenscape_2070.png',
       width = 9, height = 3, unit = "in")


#modify df for 3d plots
mig_3d_avg <- mig_dat_avg
mig_3d_avg$model <- ordered(mig_3d_avg$model, 
                            levels = c("base", "CNRM-CM5", "HadGEM2-ES365", "inmcm4", "IPSL-CM5A-MR", "MIROC-ESM-CHEM") %>% rev,
                            labels = c("Baseline", "CNRM-CM5", "HadGEM2-ES365", "inmcm4", "IPSL-CM5A-MR", "MIROC-ESM-CHEM") %>% rev)
mig_3d_avg$rcp <- ordered(mig_3d_avg$rcp,
                          levels = c('Baseline', 'BIT 4.5', 'BIT 8.5', 'rcp45', 'rcp85') %>% rev,
                          labels = c('Baseline', 'BIT 4.5', 'BIT 8.5', 'RCP 4.5', 'RCP 8.5') %>% rev)

#set camera
camera <- list(eye = list(x = -.3, y = -2, z = .85))

#3d plot mid century
s3d <- plot_ly(mig_3d_avg[mig_3d_avg$period != 2070,],
               x = ~order_avg, y = ~duration_avg, z = ~rate_avg, color = ~model,
               colors = c('#BBBBBB', '#4477AA', '#228833', '#CCBB44', '#EE6677', '#AA3377') %>% rev,
               symbol = ~rcp,
               symbols = c('circle', 'square', 'diamond', 'square-open', 'diamond-open') %>% rev,
               legendgroup = ~rcp)
s3d <- s3d %>% add_markers()
s3d <- s3d %>% layout(title = 'Average Across Herds Mid-Century (2040-2069)',
                      scene = list(xaxis = list(title = 'Order (cor)',
                                                range = c(0, 0.75),
                                                dtick = 0.25,
                                                tickfont = list(size = 15),
                                                titlefont = list(size = 20),
                                                gridwidth = 2,
                                                gridcolor = toRGB("gray50")),
                                   yaxis = list(title = 'Duration (Days)',
                                                range = c(0, 16),
                                                dtick = 4,
                                                tickfont = list(size = 15),
                                                titlefont = list(size = 20),
                                                gridwidth = 2,
                                                gridcolor = toRGB("gray50")),
                                   zaxis = list(title = 'Spring Scale (Days)',
                                                range = c(12, 18),
                                                dtick = 2,
                                                tickfont = list(size = 15),
                                                titlefont = list(size = 20),
                                                gridwidth = 2,
                                                gridcolor = toRGB("gray50")),
                                   camera = camera),
                      showlegend = FALSE)

s3d
  
#3d plot end century
s3d <- plot_ly(mig_3d_avg[mig_3d_avg$period != 2040,],
               x = ~order_avg, y = ~duration_avg, z = ~rate_avg, color = ~model,
               colors = c('#BBBBBB', '#4477AA', '#228833', '#CCBB44', '#EE6677', '#AA3377') %>% rev,
               symbol = ~rcp,
               symbols = c('circle', 'square', 'diamond', 'square-open', 'diamond-open') %>% rev,
               legendgroup = ~rcp)
s3d <- s3d %>% add_markers()
s3d <- s3d %>% layout(title = 'Average Across Herds End of Century (2070-2099)',
                      scene = list(xaxis = list(title = 'Order (cor)',
                                                range = c(0, 0.75),
                                                dtick = 0.25,
                                                tickfont = list(size = 15),
                                                titlefont = list(size = 20),
                                                gridwidth = 2,
                                                gridcolor = toRGB("gray50")),
                                   yaxis = list(title = 'Duration (Days)',
                                                range = c(0, 16),
                                                dtick = 4,
                                                tickfont = list(size = 15),
                                                titlefont = list(size = 20),
                                                gridwidth = 2,
                                                gridcolor = toRGB("gray50")),
                                   zaxis = list(title = 'Spring Scale (Days)',
                                                range = c(12, 18),
                                                dtick = 2,
                                                tickfont = list(size = 15),
                                                titlefont = list(size = 20),
                                                gridwidth = 2,
                                                gridcolor = toRGB("gray50")),
                                   camera = camera),
                      showlegend = FALSE)
s3d

#########################################
###TABLES AND GRAPHICS FOR ANNUAL DATA###
#########################################

#create df of mean values across herds and scenarios
ann_dat <- data.frame(herd = character(), 
                      model = character(), 
                      rcp = character(), 
                      year = numeric(),
                      duration = numeric(), 
                      order = numeric(),
                      rate = numeric())

#create column for unique model_year
dur_ann$model_year <- str_c(dur_ann$model, '_', dur_ann$year)
order_ann$model_year <- str_c(order_ann$model, '_', order_ann$year)
rate_ann$model_year <- str_c(rate_ann$model, '_', rate_ann$year)

#fill df with values from herds
for(i in herds){
  
  #subset by year and rcp and calc the mean dur
  dur_h_45 <- dur_ann[str_detect(dur_ann$global_id, i) == T & dur_ann$rcp != 'rcp85',] %>% group_by(model_year) %>%
    mutate(duration = mean(duration, na.rm = T))
  
  dur_h_85 <- dur_ann[str_detect(dur_ann$global_id, i) == T & dur_ann$rcp != 'rcp45',] %>% group_by(model_year) %>% 
    mutate(duration = mean(duration, na.rm = T))
  
  #subset by year and calc the mean order
  order_h_45 <- order_ann[str_detect(order_ann$global_id, i) == T & order_ann$rcp != 'rcp85',] %>% group_by(model_year) %>% 
    mutate(order = mean(order, na.rm = T))
  
  order_h_85 <- order_ann[str_detect(order_ann$global_id, i) == T & order_ann$rcp != 'rcp45',] %>% group_by(model_year) %>% 
    mutate(order = mean(order, na.rm = T))
  
  #subset by year and calc the mean rate
  rate_h_45 <- rate_ann[str_detect(rate_ann$global_id, i) == T & rate_ann$rcp != 'rcp85',] %>% group_by(model_year) %>% 
    mutate(rate = mean(rate, na.rm = T))
  
  rate_h_85 <- rate_ann[str_detect(rate_ann$global_id, i) == T & rate_ann$rcp != 'rcp45',] %>% group_by(model_year) %>% 
    mutate(order = mean(rate, na.rm = T))
  
  #remove year NAs for now
  dur_h_45 <- na.omit(dur_h_45)
  dur_h_85 <- na.omit(dur_h_85)
  order_h_45 <- na.omit(order_h_45)
  order_h_85 <- na.omit(order_h_85)
  rate_h_45 <- na.omit(rate_h_45)
  rate_h_85 <- na.omit(rate_h_85)
  
  #take one row for each unique year and rcp combo
  for(j in unique(dur_h_45$year)){
    for(k in unique(dur_h_45$rcp)){
      
      if(NROW(dur_h_45[dur_h_45$year == j & dur_h_45$rcp == k,]) != 0){
        if(k == 'baseline'){
          
          #only rbind baseline once
          ann_dat <- rbind(ann_dat, 
                           data.frame(herd = i,
                                      model = dur_ann$model[dur_h_45$year == j & dur_h_45$rcp == k][1],
                                      rcp = k,
                                      year = j,
                                      duration = dur_h_45$duration[dur_h_45$year == j & dur_h_45$rcp == k][1],
                                      order = order_h_45$order[order_h_45$year == j & order_h_45$rcp == k][1],
                                      rate = rate_h_45$rate[rate_h_45$year == j & rate_h_45$rcp == k][1]))
          
        }else{
          
          #rbind both rcps main df
          ann_dat <- rbind(ann_dat, 
                           data.frame(herd = i,
                                      model = dur_ann$model[dur_h_45$year == j & dur_h_45$rcp == k][1],
                                      rcp = k,
                                      year = j,
                                      duration = dur_h_45$duration[dur_h_45$year == j & dur_h_45$rcp == k][1],
                                      order = order_h_45$order[order_h_45$year == j & order_h_45$rcp == k][1],
                                      rate = rate_h_45$rate[rate_h_45$year == j & rate_h_45$rcp == k][1]),
                           data.frame(herd = i,
                                      model = dur_ann$model[dur_h_85$year == j & dur_h_85$rcp == 'rcp85'][1],
                                      rcp = 'rcp85',
                                      year = j,
                                      duration = dur_h_85$duration[dur_h_85$year == j & dur_h_85$rcp == 'rcp85'][1],
                                      order = order_h_85$order[order_h_85$year == j & order_h_85$rcp == 'rcp85'][1],
                                      rate = rate_h_85$rate[rate_h_85$year == j & rate_h_85$rcp == 'rcp85'][1]))
        }

      }
        }
  }
}

#set BIT rcp names
ann_dat$rcp[ann_dat$rcp == 'rcp45' & ann_dat$year < 2020] <- 'BIT 4.5'
ann_dat$rcp[ann_dat$rcp == 'rcp85' & ann_dat$year < 2020] <- 'BIT 8.5'

#set color scheme for rcp scenarios
ann_col <- c('baseline' = '#228833', 
             'BIT 4.5' = '#CCBB44', 
             'BIT 8.5' = '#AA3377',
             'rcp45' = '#4477AA',
             'rcp85' = '#EE6677')

#c('#BBBBBB', '#4477AA', '#228833', '#CCBB44', '#EE6677', 
#  '#AA3377')

# #line graph of duration
# ggplot(data = ann_dat, aes(x = year, y = duration, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_point() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   geom_smooth(method = 'lm', se = F) + xlab('Year') + ylab('Duration') + ggtitle('Average Duration Middle Scenario (2000-2099)') +
#   geom_vline(xintercept = 2020)
# 
# #line graph of order
# ggplot(data = ann_dat, aes(x = year, y = order, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_point() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   geom_smooth(method = 'lm', se = F) + xlab('Year') + ylab('Order') + ggtitle('Average Order Middle Scenario (2000-2099)') +
#   geom_vline(xintercept = 2020)
# 
# #calculate mean and quartile of baseline duration
# dur_base <- dur_ann[dur_ann$rcp == 'baseline', ]
# 
# dur_base_med <- median(dur_base$duration, na.rm = T)
# dur_base_q3 <- quantile(dur_base$duration, probs = 0.75, na.rm = T)
# dur_base_q1 <- quantile(dur_base$duration, probs = 0.25, na.rm = T)
# 
# #line graph of duration with hline
# ggplot(data = ann_dat, aes(x = year, y = duration, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_point() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   geom_smooth(method = 'lm', se = F, size = 0.8) + xlab('Year') + ylab('Duration') + ggtitle('Average Duration Middle Scenario (2000-2099)') + 
#   geom_hline(yintercept = dur_base_med, colour="#000000", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = dur_base_q3, colour="#BBBBBB", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = dur_base_q1, colour="#BBBBBB", linetype="dashed", size = 1.2)
# 
# #calc number of points each rcp in quartiles
# #all data points
# ann_dat_45 <- ann_dat[ann_dat$rcp == 'rcp45',] %>% na.omit
# ann_dat_85 <- ann_dat[ann_dat$rcp == 'rcp85',] %>% na.omit
# 
# #points > baseline q3
# NROW(ann_dat_45[ann_dat_45$duration > dur_base_q3,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$duration > dur_base_q3,]) / NROW(ann_dat_85) * 100
# 
# #points in box
# NROW(ann_dat_45[ann_dat_45$duration <= dur_base_q3 & ann_dat_45$duration >= dur_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$duration <= dur_base_q3 & ann_dat_85$duration >= dur_base_q1,]) / NROW(ann_dat_85) * 100
# 
# #points < baseline q1
# NROW(ann_dat_45[ann_dat_45$duration < dur_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$duration < dur_base_q1,]) / NROW(ann_dat_85) * 100
# 
# #calculate mean and quartile of baseline order
# order_base <- order_ann[order_ann$rcp == 'baseline', ]
# 
# order_base_med <- median(order_base$order, na.rm = T)
# order_base_q3 <- quantile(order_base$order, probs = 0.75, na.rm = T)
# order_base_q1 <- quantile(order_base$order, probs = 0.25, na.rm = T)
# 
# #line graph of order with hline
# ggplot(data = ann_dat, aes(x = year, y = order, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_point() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   geom_smooth(method = 'lm', se = F, size = 0.8) + xlab('Year') + ylab('Order') + ggtitle('Average Order Middle Scenario (2000-2099)') +
#   geom_hline(yintercept = order_base_med, colour="#000000", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = order_base_q3, colour="#BBBBBB", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = order_base_q1, colour="#BBBBBB", linetype="dashed", size = 1.2)
# 
# #calc number of points each rcp in quartiles
# #points > baseline q3
# NROW(ann_dat_45[ann_dat_45$order > order_base_q3,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$order > order_base_q3,]) / NROW(ann_dat_85) * 100
# 
# #points in box
# NROW(ann_dat_45[ann_dat_45$order <= order_base_q3 & ann_dat_45$order >= order_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$order <= order_base_q3 & ann_dat_85$order >= order_base_q1,]) / NROW(ann_dat_85) * 100
# 
# #points < baseline q1
# NROW(ann_dat_45[ann_dat_45$order < order_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$order < order_base_q1,]) / NROW(ann_dat_85) * 100

#create ann data average across herds
ann_dat_avg <- ann_dat

#create unique scenarios
ann_dat_avg$scenario <- str_c(ann_dat_avg$rcp, '_', ann_dat_avg$year)

#take average order, duration and rate
ann_dat_avg <- ann_dat_avg %>% group_by(scenario) %>%
  mutate(duration_avg = mean(duration, na.rm = T)) %>% 
  mutate(order_avg = mean(order, na.rm = T)) %>%
  mutate(rate_avg = mean(rate, na.rm = T))

#remove duplicates
ann_dat_avg <- ann_dat_avg[!duplicated(ann_dat_avg$scenario),]

#remove unnecessary columns to limit confusion
ann_dat_avg <- ann_dat_avg[!(colnames(ann_dat_avg) %in% c('herd', 'model', 'duration', 'order', 'rate'))]

#bring in drought data to assess relationships between drought and variability
#rbind all data together into large df for 2 RCPS
pdsi_45 <- rbind(read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp45_1999.csv'),
                 read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp45_2020.csv'),
                 read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp45_2040.csv'),
                 read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp45_2070.csv'))

pdsi_85 <- rbind(read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp85_1999.csv'),
                 read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp85_2020.csv'),
                 read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp85_2040.csv'),
                 read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp85_2070.csv'))

#add rcp col
pdsi_45$rcp <- 'rcp45'
pdsi_85$rcp <- 'rcp85'

#take the average annual pdsi value across state
pdsi_45 <- pdsi_45 %>% group_by(year) %>%
  mutate(pdsi_avg = mean(pdsi, na.rm = T))

pdsi_85 <- pdsi_85 %>% group_by(year) %>%
  mutate(pdsi_avg = mean(pdsi, na.rm = T))

#remove duplicates
pdsi_45 <- pdsi_45[!duplicated(pdsi_45$pdsi_avg),]
pdsi_85 <- pdsi_85[!duplicated(pdsi_85$pdsi_avg),]

#bind together
pdsi <- rbind(data.frame(year = pdsi_45$year, rcp = pdsi_45$rcp, pdsi = pdsi_45$pdsi_avg),
              data.frame(year = pdsi_85$year, rcp = pdsi_85$rcp, pdsi = pdsi_85$pdsi_avg))
rm(pdsi_45, pdsi_85)

#remove 1999 and round
pdsi <- pdsi[pdsi$year != 1999,]
pdsi$pdsi <- round(pdsi$pdsi, 2)

#change RCP category for BIT
pdsi$rcp[pdsi$year < 2020 & pdsi$rcp == 'rcp45'] <- 'BIT 4.5'
pdsi$rcp[pdsi$year < 2020 & pdsi$rcp == 'rcp85'] <- 'BIT 8.5'

#add pdsi to ann_dat_avg df
ann_dat_avg$pdsi <- NA
for(i in 1:NROW(pdsi)){
  ann_dat_avg$pdsi[ann_dat_avg$year == pdsi$year[i] & ann_dat_avg$rcp == pdsi$rcp[i]] <- pdsi$pdsi[i]
}

#check correlation between pdsi and greenscape metrics
cor(ann_dat_avg$pdsi[ann_dat_avg$rcp != 'baseline'], ann_dat_avg$duration_avg[ann_dat_avg$rcp != 'baseline'], use = 'complete.obs')
cor(ann_dat_avg$pdsi[ann_dat_avg$rcp != 'baseline'], ann_dat_avg$order_avg[ann_dat_avg$rcp != 'baseline'], use = 'complete.obs')
cor(ann_dat_avg$pdsi[ann_dat_avg$rcp != 'baseline'], ann_dat_avg$rate_avg[ann_dat_avg$rcp != 'baseline'], use = 'complete.obs')

#load contemporary pdsi data
load(file = '/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/dat2_2021_01_27.RData')

#drop everything but year and pdsi
dat2 <- dat2[,colnames(dat2) %in% c('year', 'pdsi_mar_apr_min')]

#create annual averages
dat2 <- dat2 %>% group_by(year) %>%
  mutate(pdsi_avg = round(mean(pdsi_mar_apr_min, na.rm = T), 2))

#remove duplicates
dat2 <- dat2[!duplicated(dat2$pdsi_avg),] 

#add contemporary pdsi to ann_dat_avg df
for(i in 1:NROW(dat2)){
  ann_dat_avg$pdsi[ann_dat_avg$year == dat2$year[i] & ann_dat_avg$rcp == 'baseline'] <- dat2$pdsi_avg[i]
}

#check baseline correlation between greenscape metrics and pdsi
cor(ann_dat_avg$pdsi[ann_dat_avg$rcp == 'baseline'], ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'baseline'], use = 'complete.obs')
cor(ann_dat_avg$pdsi[ann_dat_avg$rcp == 'baseline'], ann_dat_avg$order_avg[ann_dat_avg$rcp == 'baseline'], use = 'complete.obs')
cor(ann_dat_avg$pdsi[ann_dat_avg$rcp == 'baseline'], ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'baseline'], use = 'complete.obs')

#plot average annual drought by year
ggplot(data = ann_dat_avg, aes(x = year, y = pdsi, color = rcp)) +
  geom_point() + theme_bw()

 #average order duration and rate under different drought scenarios
pdsi_df <- data.frame(rcp = character(), scenario = character(), duration = numeric(),
                      order = numeric(), rate = numeric(), n = numeric(), total = numeric(), perc = numeric())

scen <- c('moderate drought', 'normal', 'unusual moist spell')

#loop through rcp and drought scenarios
for(r in unique(ann_dat_avg$rcp)){
  for(s in scen){
    #create subset of rcp we want
    sub <- ann_dat_avg[ann_dat_avg$rcp == r,]
    
    #remove missing rows
    sub <- na.omit(sub)
    
    #calc relavent rows
    if(s == 'moderate drought'){
      sub_c <- sub[sub$pdsi <= -2,]
    }else if(s == 'normal'){
      sub_c <- sub[sub$pdsi > -2 & sub$pdsi < 2,]
    }else{sub_c <- sub[sub$pdsi >= 2,]}
    
    if(NROW(sub_c != 0)){
      pdsi_df <- rbind(pdsi_df,
                       data.frame(
                         rcp = r,
                         scenario = s,
                         duration = mean(sub_c$duration_avg) %>% round(2),
                         order = mean(sub_c$order_avg) %>% round(2),
                         rate = mean(sub_c$rate_avg) %>% round(2),
                         n = NROW(sub_c),
                         total = NROW(sub),
                         perc = round((NROW(sub_c)/NROW(sub))*100, 2)))
    }else{
      pdsi_df <- rbind(pdsi_df,
                       data.frame(
                         rcp = r,
                         scenario = s,
                         duration = NA,
                         order = NA,
                         rate = NA,
                         n = 0,
                         total = NROW(sub),
                         perc = 0))
    }
  }
}

#write csv to disk
write.csv(pdsi_df, file = 'output/final_plots/tables/annual_drought_greenscape.csv', row.names = F)

#set annual theme
theme_ann <- theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                 axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                 legend.title = element_text(size = 15), legend.position = 'none')

#lm of duration, order, ss
lm(duration_avg ~ year, data = subset(ann_dat_avg, rcp == 'rcp45')) %>% summary
lm(duration_avg ~ year, data = subset(ann_dat_avg, rcp == 'rcp85')) %>% summary
lm(order_avg ~ year, data = subset(ann_dat_avg, rcp == 'rcp45')) %>% summary
lm(order_avg ~ year, data = subset(ann_dat_avg, rcp == 'rcp85')) %>% summary
lm(rate_avg ~ year, data = subset(ann_dat_avg, rcp == 'rcp45')) %>% summary
lm(rate_avg ~ year, data = subset(ann_dat_avg, rcp == 'rcp85')) %>% summary

#line graph of duration
ggplot(data = ann_dat_avg, aes(x = year, y = duration_avg, color = rcp)) +
  #  geom_line(aes(color = rcp)) +
  geom_point() +
  theme_bw() + theme_ann +
  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  geom_smooth(data = subset(ann_dat_avg, rcp == 'rcp45' | rcp == 'rcp85'),
              aes(x = year, y = duration_avg, color = rcp, linetype = 'Trend'), method = 'lm', se = F) + 
  xlab('Year') + ylab('Duration (Days)') + ggtitle('Annual Greenscape Metrics Middle Model (2000-2099)') +
  geom_vline(xintercept = 2020) +
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'baseline']),
                   yend = mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'baseline']),
                   linetype = 'Mean',  color = 'baseline'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'BIT 4.5']),
                   yend = mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'BIT 4.5']),
                   linetype = 'Mean', color = 'BIT 4.5'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'BIT 8.5']),
                   yend = mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'BIT 8.5']),
                   linetype = 'Mean',  color = 'BIT 8.5'), lwd = 0.8) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2, keywidth = 1.8)) +
  scale_linetype_manual('Line Type', values = c('Trend' = 1, 'Mean' = 2))

ggsave('output/final_plots/avg_duration_annual_lm.png', width = 8, height = 3.5)

#line graph of order
ggplot(data = ann_dat_avg, aes(x = year, y = order_avg, color = rcp)) +
  #  geom_line(aes(color = rcp)) +
  geom_point() +
  theme_bw() + theme_ann +
  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  geom_smooth(data = subset(ann_dat_avg, rcp == 'rcp45' | rcp == 'rcp85'),
              aes(x = year, y = order_avg, color = rcp, linetype = 'Trend'), method = 'lm', se = F) + 
  xlab('Year') + ylab('Order (Correlation)') + ggtitle('Annual Order Middle Scenario (2000-2099)') +
  geom_vline(xintercept = 2020) +
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'baseline']),
                   yend = mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'baseline']),
                   linetype = 'Mean',  color = 'baseline'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'BIT 4.5']),
                   yend = mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'BIT 4.5']),
                   linetype = 'Mean', color = 'BIT 4.5'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'BIT 8.5']),
                   yend = mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'BIT 8.5']),
                   linetype = 'Mean',  color = 'BIT 8.5'), lwd = 0.8) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2, keywidth = 1.8)) +
  scale_linetype_manual('Line Type', values = c('Trend' = 1, 'Mean' = 2))

ggsave('output/final_plots/avg_order_annual_lm.png', width = 8, height = 3.5)

#line graph of rate
ggplot(data = ann_dat_avg, aes(x = year, y = rate_avg, color = rcp)) +
  #  geom_line(aes(color = rcp)) +
  geom_point() +
  theme_bw() + theme_ann +
  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  geom_smooth(data = subset(ann_dat_avg, rcp == 'rcp45' | rcp == 'rcp85'),
              aes(x = year, y = rate_avg, color = rcp, linetype = 'Trend'), method = 'lm', se = F) + 
  xlab('Year') + ylab('Mean Spring Scale (Days)') + ggtitle('Spring Scale Middle Scenario (2000-2099)') +
  geom_vline(xintercept = 2020) + scale_y_continuous(limits = c(8, 24), breaks = c(8, 12, 16, 20, 24)) +
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'baseline']),
                   yend = mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'baseline']),
                   linetype = 'Mean',  color = 'baseline'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'BIT 4.5']),
                   yend = mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'BIT 4.5']),
                   linetype = 'Mean', color = 'BIT 4.5'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'BIT 8.5']),
                   yend = mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'BIT 8.5']),
                   linetype = 'Mean',  color = 'BIT 8.5'), lwd = 0.8) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2, keywidth = 1.8)) +
  scale_linetype_manual('Line Type', values = c('Trend' = 1, 'Mean' = 2))

ggsave('output/final_plots/avg_rate_annual_lm.png', width = 8, height = 3.5)

#calculate number of years that at least one greenscape metric has "negative" outcome

#allocate output df
ann_results <- data.frame(metric = character(), num_years = numeric(), perc_years = numeric())

#based on baseline mean
#create new df to populate
ann_outcome <- ann_dat_avg

#duration - true if duration less than contemp
ann_outcome$dur_out <- (ann_dat_avg$duration_avg - mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'baseline'])) < 0

#order - true if order less than contemp
ann_outcome$ord_out <- (ann_dat_avg$order_avg - mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'baseline'])) < 0

#mean spring scale - true if mean spring scale greater than contemp
ann_outcome$ss_out <- (ann_dat_avg$rate_avg - mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'baseline'])) > 0

#calculate years that at least one is true
ann_outcome$neg <- (ann_outcome$dur_out + ann_outcome$ord_out + ann_outcome$ss_out) >= 1

#populate output df
ann_results <- rbind(ann_results,
                     data.frame(metric = 'baseline to rcp 4.5',
                                num_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp45']),
                                perc_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp45'])/80),
                     data.frame(metric = 'baseline to rcp 8.5',
                                num_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp85']),
                                perc_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp85'])/80))

#based on BIT 4.5 mean
#create new df to populate
ann_outcome <- ann_dat_avg

#duration - true if duration less than contemp
ann_outcome$dur_out <- (ann_dat_avg$duration_avg - mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'BIT 4.5'])) < 0

#order - true if order less than contemp
ann_outcome$ord_out <- (ann_dat_avg$order_avg - mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'BIT 4.5'])) < 0

#mean spring scale - true if mean spring scale greater than contemp
ann_outcome$ss_out <- (ann_dat_avg$rate_avg - mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'BIT 4.5'])) > 0

#calculate years that at least one is true
ann_outcome$neg <- (ann_outcome$dur_out + ann_outcome$ord_out + ann_outcome$ss_out) >= 1

#populate output df
ann_results <- rbind(ann_results,
                     data.frame(metric = 'bit 4.5 to rcp 4.5',
                                num_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp45']),
                                perc_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp45'])/80))

#based on BIT 8.5 mean
#create new df to populate
ann_outcome <- ann_dat_avg

#duration - true if duration less than contemp
ann_outcome$dur_out <- (ann_dat_avg$duration_avg - mean(ann_dat_avg$duration_avg[ann_dat_avg$rcp == 'BIT 8.5'])) < 0

#order - true if order less than contemp
ann_outcome$ord_out <- (ann_dat_avg$order_avg - mean(ann_dat_avg$order_avg[ann_dat_avg$rcp == 'BIT 8.5'])) < 0

#mean spring scale - true if mean spring scale greater than contemp
ann_outcome$ss_out <- (ann_dat_avg$rate_avg - mean(ann_dat_avg$rate_avg[ann_dat_avg$rcp == 'BIT 8.5'])) > 0

#calculate years that at least one is true
ann_outcome$neg <- (ann_outcome$dur_out + ann_outcome$ord_out + ann_outcome$ss_out) >= 1

#populate output df
ann_results <- rbind(ann_results,
                     data.frame(metric = 'bit 8.5 to rcp 8.5',
                                num_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp85']),
                                perc_years = sum(ann_outcome$neg[ann_outcome$rcp == 'rcp85'])/80))

 ##########################
###LOOK INTO HERD LINES###
##########################

#line graph of duration
ggplot(data = subset(ann_dat, rcp == 'BIT 4.5' | rcp == 'rcp45'), aes(x = year, y = duration, color = herd)) +
  geom_line() +
  theme_bw() + theme_e +
#  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  xlab('Year') + ylab('Duration') + ggtitle('Duration Middle Model RCP 4.5 (2000-2099)') +
  geom_vline(xintercept = 2020)

ggsave('output/final_plots/duration_annual_byherd_45.png')

#8.5
ggplot(data = subset(ann_dat, rcp == 'BIT 8.5' | rcp == 'rcp85'), aes(x = year, y = duration, color = herd)) +
  geom_line() +
  theme_bw() + theme_e +
  #  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  xlab('Year') + ylab('Duration') + ggtitle('Duration Middle Model RCP 8.5 (2000-2099)') +
  geom_vline(xintercept = 2020)

ggsave('output/final_plots/duration_annual_byherd_85.png')

#line graph of order
ggplot(data = subset(ann_dat, rcp == 'BIT 4.5' | rcp == 'rcp45'), aes(x = year, y = order, color = herd)) +
  geom_line() +
  theme_bw() + theme_e +
  #  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  xlab('Year') + ylab('Order') + ggtitle('Order Middle Model RCP 4.5 (2000-2099)') +
  geom_vline(xintercept = 2020)

ggsave('output/final_plots/order_annual_byherd_45.png')

#8.5
ggplot(data = subset(ann_dat, rcp == 'BIT 8.5' | rcp == 'rcp85'), aes(x = year, y = order, color = herd)) +
  geom_line() +
  theme_bw() + theme_e +
  #  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  xlab('Year') + ylab('Order') + ggtitle('Order Middle Model RCP 8.5 (2000-2099)') +
  geom_vline(xintercept = 2020)

ggsave('output/final_plots/order_annual_byherd_85.png')

#line graph of rate
ggplot(data = subset(ann_dat, rcp == 'BIT 4.5' | rcp == 'rcp45'), aes(x = year, y = rate, color = herd)) +
  geom_line() +
  theme_bw() + theme_e +
  #  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  xlab('Year') + ylab('Spring Scale') + ggtitle('Mean Spring Scale Middle Model RCP 4.5 (2000-2099)') +
  geom_vline(xintercept = 2020)

ggsave('output/final_plots/ss_annual_byherd_45.png')

#8.5
ggplot(data = subset(ann_dat, rcp == 'BIT 8.5' | rcp == 'rcp85'), aes(x = year, y = rate, color = herd)) +
  geom_line() +
  theme_bw() + theme_e +
  #  scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
  xlab('Year') + ylab('Spring Scale') + ggtitle('Mean Spring Scale Middle Model RCP 8.5 (2000-2099)') +
  geom_vline(xintercept = 2020)

ggsave('output/final_plots/ss_annual_byherd_85.png')

# #calculate mean and quartile of baseline duration
# dur_base <- dur_ann[dur_ann$rcp == 'baseline', ]
# 
# dur_base_med <- median(dur_base$duration, na.rm = T)
# dur_base_q3 <- quantile(dur_base$duration, probs = 0.75, na.rm = T)
# dur_base_q1 <- quantile(dur_base$duration, probs = 0.25, na.rm = T)
# 
# #line graph of duration with hline
# ggplot(data = ann_dat_avg, aes(x = year, y = duration_avg, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_point() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   geom_smooth(method = 'lm', se = F, size = 0.8) + xlab('Year') + ylab('Duration') + ggtitle('Average Duration Middle Scenario (2000-2099)') + 
#   geom_hline(yintercept = dur_base_med, colour="#000000", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = dur_base_q3, colour="#BBBBBB", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = dur_base_q1, colour="#BBBBBB", linetype="dashed", size = 1.2)
# 
# #calc number of points each rcp in quartiles
# #all data points
# ann_dat_45 <- ann_dat[ann_dat$rcp == 'rcp45',] %>% na.omit
# ann_dat_85 <- ann_dat[ann_dat$rcp == 'rcp85',] %>% na.omit
# 
# #points > baseline q3
# NROW(ann_dat_45[ann_dat_45$duration > dur_base_q3,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$duration > dur_base_q3,]) / NROW(ann_dat_85) * 100
# 
# #points in box
# NROW(ann_dat_45[ann_dat_45$duration <= dur_base_q3 & ann_dat_45$duration >= dur_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$duration <= dur_base_q3 & ann_dat_85$duration >= dur_base_q1,]) / NROW(ann_dat_85) * 100
# 
# #points < baseline q1
# NROW(ann_dat_45[ann_dat_45$duration < dur_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$duration < dur_base_q1,]) / NROW(ann_dat_85) * 100
# 
# #calculate mean and quartile of baseline order
# order_base <- order_ann[order_ann$rcp == 'baseline', ]
# 
# order_base_med <- median(order_base$order, na.rm = T)
# order_base_q3 <- quantile(order_base$order, probs = 0.75, na.rm = T)
# order_base_q1 <- quantile(order_base$order, probs = 0.25, na.rm = T)
# 
# #line graph of order with hline
# ggplot(data = ann_dat_avg, aes(x = year, y = order_avg, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_point() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('baseline' = 'Baseline', 'rcp45' = 'RCP 4.5', 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   geom_smooth(method = 'lm', se = F, size = 0.8) + xlab('Year') + ylab('Order') + ggtitle('Average Order Middle Scenario (2000-2099)') +
#   geom_hline(yintercept = order_base_med, colour="#000000", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = order_base_q3, colour="#BBBBBB", linetype="dashed", size = 1.2) +
#   geom_hline(yintercept = order_base_q1, colour="#BBBBBB", linetype="dashed", size = 1.2)
# 
# #calc number of points each rcp in quartiles
# #points > baseline q3
# NROW(ann_dat_45[ann_dat_45$order > order_base_q3,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$order > order_base_q3,]) / NROW(ann_dat_85) * 100
# 
# #points in box
# NROW(ann_dat_45[ann_dat_45$order <= order_base_q3 & ann_dat_45$order >= order_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$order <= order_base_q3 & ann_dat_85$order >= order_base_q1,]) / NROW(ann_dat_85) * 100
# 
# #points < baseline q1
# NROW(ann_dat_45[ann_dat_45$order < order_base_q1,]) / NROW(ann_dat_45) * 100
# NROW(ann_dat_85[ann_dat_85$order < order_base_q1,]) / NROW(ann_dat_85) * 100
# 
# 

# #graphs of variability across herds
# #line graph of duration rcp 4.5
# ggplot(data = ann_dat[ann_dat$rcp %in% c('BIT 4.5', 'rcp45'),], aes(x = factor(year), y = duration, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_boxplot() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('rcp45' = 'RCP 4.5'), name = 'RCP') +
#   #geom_smooth(method = 'lm', se = F) + 
#   xlab('Year') + ylab('Duration (Days)') + ggtitle('Duration Middle Scenario RCP 4.5 (2000-2099)') +
#   geom_vline(xintercept = factor(2020)) + 
#   scale_x_discrete(breaks = c('2000', '2020', '2040', '2060', '2080', '2099')) +
#   scale_y_continuous(limits = c(-20, 50), breaks = c(-20, -10, 0, 10, 20, 30, 40, 50))
# 
# ggsave('output/final_plots/duration_annual_rcp45.png')
# 
# #line graph of duration rcp 4.5
# ggplot(data = ann_dat[ann_dat$rcp %in% c('BIT 8.5', 'rcp85'),], aes(x = factor(year), y = duration, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_boxplot() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('rcp85' = 'RCP 8.5'), name = 'RCP') +
#   #geom_smooth(method = 'lm', se = F) + 
#   xlab('Year') + ylab('Duration (Days)') + ggtitle('Duration Middle Scenario RCP 8.5 (2000-2099)') +
#   geom_vline(xintercept = factor(2020)) + 
#   scale_x_discrete(breaks = c('2000', '2020', '2040', '2060', '2080', '2099')) +
#   scale_y_continuous(limits = c(-20, 50), breaks = c(-20, -10, 0, 10, 20, 30, 40, 50))
# 
# ggsave('output/final_plots/duration_annual_rcp85.png')
# 
# #line graph of order rcp 4.5
# ggplot(data = ann_dat[ann_dat$rcp %in% c('BIT 4.5', 'rcp45'),], aes(x = factor(year), y = order, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_boxplot() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c( 'rcp45' = 'RCP 4.5'), name = 'RCP') +
#   #geom_smooth(method = 'lm', se = F) + 
#   xlab('Year') + ylab('Order (Cor)') + ggtitle('Order Middle Scenario RCP 4.5 (2000-2099)') +
#   geom_vline(xintercept = factor(2020)) + 
#   scale_x_discrete(breaks = c('2000', '2020', '2040', '2060', '2080', '2099')) +
#   scale_y_continuous(limits = c(-1, 1))
# 
# ggsave('output/final_plots/order_annual_rcp45.png')
# 
# #line graph of order rcp 8.5
# ggplot(data = ann_dat[ann_dat$rcp %in% c('BIT 8.5', 'rcp85'),], aes(x = factor(year), y = order, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_boxplot() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c( 'rcp85' = 'RCP 8.5'), name = 'RCP') +
#   #geom_smooth(method = 'lm', se = F) + 
#   xlab('Year') + ylab('Order (Cor)') + ggtitle('Order Middle Scenario RCP 8.5 (2000-2099)') +
#   geom_vline(xintercept = factor(2020)) + 
#   scale_x_discrete(breaks = c('2000', '2020', '2040', '2060', '2080', '2099')) +
#   scale_y_continuous(limits = c(-1, 1))
# 
# ggsave('output/final_plots/order_annual_rcp85.png')
# 
# #line graph of rate rcp 4.5
# ggplot(data = ann_dat[ann_dat$rcp %in% c('BIT 4.5', 'rcp45'),], aes(x = factor(year), y = rate, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_boxplot() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('rcp45' = 'RCP 4.5'), name = 'RCP') +
#   #geom_smooth(method = 'lm', se = F) + 
#   xlab('Year') + ylab('Spring Scale (Days)') + ggtitle('Rate Middle Scenario RCP 4.5 (2000-2099)') +
#   geom_vline(xintercept = factor(2020)) + 
#   scale_x_discrete(breaks = c('2000', '2020', '2040', '2060', '2080', '2099')) +
#   scale_y_continuous(limits = c(5, 30), breaks = c(5, 10, 15, 20, 25, 30))
# 
# ggsave('output/final_plots/rate_annual_rcp45.png')
# 
# #line graph of rate rcp 8.5
# ggplot(data = ann_dat[ann_dat$rcp %in% c('BIT 8.5', 'rcp85'),], aes(x = factor(year), y = rate, color = rcp)) +
#   #  geom_line(aes(color = rcp)) +
#   geom_boxplot() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, labels = c('rcp85' = 'RCP 8.5'), name = 'RCP') +
#   #geom_smooth(method = 'lm', se = F) + 
#   xlab('Year') + ylab('Spring Scale (Days)') + ggtitle('Rate Middle Scenario RCP 8.5 (2000-2099)') +
#   geom_vline(xintercept = factor(2020)) + 
#   scale_x_discrete(breaks = c('2000', '2020', '2040', '2060', '2080', '2099')) +
#   scale_y_continuous(limits = c(5, 30), breaks = c(5, 10, 15, 20, 25, 30))
# 
# ggsave('output/final_plots/rate_annual_rcp85.png')

##########################################################
###ADD GREENSCAPE METRICS TO INDIVIDUAL MIGRATION LINES###
##########################################################

#load shapefiles of migration lines
files <- list.files('mig_data/processed/shp', pattern = glob2rx('mig_lines*.shp'), full.names = T)

#remove names that have greenscape
files <- files[str_detect(files, 'greenscape') == F]

#loop through herds
for(i in files){
  
  #load herd lines shp
  shp <- readOGR(i)
  
  #population shapefile data with greenscape metrics
  #baseline
  shp$d_base_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'base') %>% select(duration) %>% round %>% .[,1]
  shp$o_base_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'base') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_base_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'base') %>% select(rate) %>% round(2) %>% .[,1]
  
  #CNRM.CM5
  #2000-2019 rcp45
  shp$d_CNR45_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_CNR45_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_CNR45_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2000-2019 rcp85
  shp$d_CNR85_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_CNR85_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_CNR85_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp45
  shp$d_CNR45_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_CNR45_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_CNR45_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp85
  shp$d_CNR85_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_CNR85_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_CNR85_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_2040') %>% select(rate) %>% round(2) %>% .[,1]

  #2070-2099 rcp45
  shp$d_CNR45_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_CNR45_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_CNR45_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp45_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp85
  shp$d_CNR85_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_CNR85_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_CNR85_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'CNRM.CM5_rcp85_2070') %>% select(rate) %>% round(2) %>% .[,1]
 
  #HadGEM2.ES365
  #2000-2019 rcp45
  shp$d_Had45_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_Had45_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_Had45_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2000-2019 rcp85
  shp$d_Had85_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_Had85_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_Had85_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp45
  shp$d_Had45_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_Had45_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_Had45_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp85
  shp$d_Had85_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_Had85_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_Had85_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp45
  shp$d_Had45_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_Had45_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_Had45_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp45_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp85
  shp$d_Had85_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_Had85_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_Had85_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'HadGEM2.ES365_rcp85_2070') %>% select(rate) %>% round(2) %>% .[,1]
   
  #inmcm4
  #2000-2019 rcp45
  shp$d_inm45_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_inm45_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_inm45_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp45
  shp$d_inm45_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_inm45_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_inm45_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp85
  shp$d_inm85_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp85_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_inm85_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp85_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_inm85_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp85_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp45
  shp$d_inm45_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_inm45_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_inm45_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp45_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp85
  shp$d_inm85_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp85_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_inm85_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp85_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_inm85_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'inmcm4_rcp85_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #IPSL.CM5A.MR
  #2000-2019 rcp45
  shp$d_IPS45_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_IPS45_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_IPS45_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2000-2019 rcp85
  shp$d_IPS85_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_IPS85_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_IPS85_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp45
  shp$d_IPS45_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_IPS45_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_IPS45_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp85
  shp$d_IPS85_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_IPS85_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_IPS85_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp45
  shp$d_IPS45_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_IPS45_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_IPS45_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp45_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp85
  shp$d_IPS85_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_IPS85_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_IPS85_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'IPSL.CM5A.MR_rcp85_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #MIROC.ESM.CHEM
  #2000-2019 rcp45
  shp$d_MIR45_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_MIR45_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_MIR45_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2000-2019 rcp85
  shp$d_MIR85_00 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_1999') %>% select(duration) %>% round %>% .[,1]
  shp$o_MIR85_00 <- order %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_1999') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_MIR85_00 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_1999') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp45
  shp$d_MIR45_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_MIR45_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_MIR45_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2040-2069 rcp85
  shp$d_MIR85_40 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_2040') %>% select(duration) %>% round %>% .[,1]
  shp$o_MIR85_40 <- order %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_2040') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_MIR85_40 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_2040') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp45
  shp$d_MIR45_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_MIR45_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_MIR45_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp45_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #2070-2099 rcp85
  shp$d_MIR85_70 <- dur %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_2070') %>% select(duration) %>% round %>% .[,1]
  shp$o_MIR85_70 <- order %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_2070') %>% select(order) %>% round(2) %>% .[,1]
  shp$s_MIR85_70 <- rate %>% filter(global_id %in% shp$GlobalID & model == 'MIROC.ESM.CHEM_rcp85_2070') %>% select(rate) %>% round(2) %>% .[,1]
  
  #export to shp
  exp <- str_replace(i, 'lines', 'lines_greenscape')
  
  raster::shapefile(shp, exp, overwrite = T)
  
}

##########################################
###CALCULATE ANNUAL CORRELATION METRICS###
##########################################

#create dataframes for annual greenscape correlation tables rcp 4.5
ann_cor_dur <- ann_dat %>% filter(rcp == 'BIT 4.5' | rcp == 'rcp45') %>% 
  select(herd, year, duration)

ann_cor_ord <- ann_dat %>% filter(rcp == 'BIT 4.5' | rcp == 'rcp45') %>% 
  select(herd, year, order)

ann_cor_ss <- ann_dat %>% filter(rcp == 'BIT 4.5' | rcp == 'rcp45') %>% 
  select(herd, year, rate)

#cast df by herd
ann_cor_dur <- dcast(ann_cor_dur, year ~ herd)
ann_cor_ord <- dcast(ann_cor_ord, year ~ herd)
ann_cor_ss <- dcast(ann_cor_ss, year ~ herd)

#change row names to year and remove year column
rownames(ann_cor_dur) <- ann_cor_dur$year
ann_cor_dur <- ann_cor_dur %>% select(!year)

rownames(ann_cor_ord) <- ann_cor_ord$year
ann_cor_ord <- ann_cor_ord %>% select(!year)

rownames(ann_cor_ss) <- ann_cor_ss$year
ann_cor_ss <- ann_cor_ss %>% select(!year)

#calculate correlation
ann_cor_dur <- cor(ann_cor_dur, use = 'pairwise.complete.obs') %>% round(2)
ann_cor_ord <- cor(ann_cor_ord, use = 'pairwise.complete.obs') %>% round(2)
ann_cor_ss <- cor(ann_cor_ss, use = 'pairwise.complete.obs') %>% round(2)

#plot correlation output pixel width 800
# corrplot(ann_cor_dur, method = 'number', type = 'upper',
#          title = 'Duration Correlation Across Herds 2000-2099 From Middle Model RCP 4.5',
#          mar=c(0,0,1,0))
# 
# corrplot(ann_cor_ord, method = 'circle', type = 'upper',
#          title = '\nOrder Correlation Across Herds 2000-2099 From Middle Model RCP 4.5')
# 
# corrplot(ann_cor_ss, method = 'circle', type = 'upper',
#          title = '\nMean Spring Scale Correlation Across Herds 2000-2099 From Middle Model RCP 4.5')

ggcorrplot(ann_cor_dur, type = "lower",
           lab = TRUE, title = 'Duration Correlation Across Herds 2000-2099 From Middle Model RCP 4.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_dur_45.png', width = 12, height = 8)

ggcorrplot(ann_cor_ord, type = "lower",
           lab = TRUE, title = 'Order Correlation Across Herds 2000-2099 From Middle Model RCP 4.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_ord_45.png', width = 12, height = 8)

ggcorrplot(ann_cor_ss, type = "lower",
           lab = TRUE, title = 'Mean Spring Scale Correlation Across Herds 2000-2099 From Middle Model RCP 4.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_ss_45.png', width = 12, height = 8)

#create dataframes for annual greenscape correlation tables rcp 4.5
ann_cor_dur <- ann_dat %>% filter(rcp == 'BIT 8.5' | rcp == 'rcp85') %>% 
  select(herd, year, duration)

ann_cor_ord <- ann_dat %>% filter(rcp == 'BIT 8.5' | rcp == 'rcp85') %>% 
  select(herd, year, order)

ann_cor_ss <- ann_dat %>% filter(rcp == 'BIT 8.5' | rcp == 'rcp85') %>% 
  select(herd, year, rate)

#cast df by herd
ann_cor_dur <- dcast(ann_cor_dur, year ~ herd)
ann_cor_ord <- dcast(ann_cor_ord, year ~ herd)
ann_cor_ss <- dcast(ann_cor_ss, year ~ herd)

#change row names to year and remove year column
rownames(ann_cor_dur) <- ann_cor_dur$year
ann_cor_dur <- ann_cor_dur %>% select(!year)

rownames(ann_cor_ord) <- ann_cor_ord$year
ann_cor_ord <- ann_cor_ord %>% select(!year)

rownames(ann_cor_ss) <- ann_cor_ss$year
ann_cor_ss <- ann_cor_ss %>% select(!year)

#calculate correlation
ann_cor_dur <- cor(ann_cor_dur, use = 'pairwise.complete.obs') %>% round(2)
ann_cor_ord <- cor(ann_cor_ord, use = 'pairwise.complete.obs') %>% round(2)
ann_cor_ss <- cor(ann_cor_ss, use = 'pairwise.complete.obs') %>% round(2)

#plot correlation output pixel width 800
# corrplot(ann_cor_dur, method = 'circle', type = 'upper',
#          title = '\nDuration Correlation Across Herds 2000-2099 From Middle Model RCP 8.5')
# 
# corrplot(ann_cor_ord, method = 'circle', type = 'upper',
#          title = '\nOrder Correlation Across Herds 2000-2099 From Middle Model RCP 8.5')
# 
# corrplot(ann_cor_ss, method = 'circle', type = 'upper',
#          title = '\nMean Spring Scale Correlation Across Herds 2000-2099 From Middle Model RCP 8.5')

ggcorrplot(ann_cor_dur, type = "lower",
           lab = TRUE, title = 'Duration Correlation Across Herds 2000-2099 From Middle Model RCP 8.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_dur_85.png', width = 12, height = 8)

ggcorrplot(ann_cor_ord, type = "lower",
           lab = TRUE, title = 'Order Correlation Across Herds 2000-2099 From Middle Model RCP 8.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_ord_85.png', width = 12, height = 8)

ggcorrplot(ann_cor_ss, type = "lower",
           lab = TRUE, title = 'Mean Spring Scale Correlation Across Herds 2000-2099 From Middle Model RCP 8.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_ss_85.png', width = 12, height = 8)

#now look into mean metrics across herds

#create dataframe rcp 4.5
ann_cor <- ann_dat_avg %>% as.data.frame %>% filter(rcp == 'BIT 4.5' | rcp == 'rcp45') %>% 
  select(year, duration_avg, order_avg, rate_avg)

#change row names to year and remove year and scenario column
rownames(ann_cor) <- ann_cor$year
ann_cor <- ann_cor %>% select(!year)

#change col names
colnames(ann_cor) <- c('duration', 'order', 'mean spring scale')

#calculate correlation
ann_cor <- cor(ann_cor, use = 'pairwise.complete.obs') %>% round(2)

#plot correlation output pixel width 800
# corrplot(ann_cor, method = 'circle', type = 'upper',
#          title = '\nCorrelation Across Greenscape Metrics 2000-2099 From Middle Model RCP 4.5')

ggcorrplot(ann_cor, type = "lower",
           lab = TRUE, title = 'Correlation Across Greenscape Metrics \n2000-2099 From Middle Model RCP 4.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_greenscape_45.png', width = 6, height = 4)

#create dataframe rcp 8.5
ann_cor <- ann_dat_avg %>% as.data.frame %>% filter(rcp == 'BIT 8.5' | rcp == 'rcp85') %>% 
  select(year, duration_avg, order_avg, rate_avg)

#change row names to year and remove year and scenario column
rownames(ann_cor) <- ann_cor$year
ann_cor <- ann_cor %>% select(!year)

#change col names
colnames(ann_cor) <- c('duration', 'order', 'mean spring scale')

#calculate correlation
ann_cor <- cor(ann_cor, use = 'pairwise.complete.obs') %>% round(2)

#plot correlation output pixel width 800
# corrplot(ann_cor, method = 'circle', type = 'upper',
#          title = '\nCorrelation Across Greenscape Metrics 2000-2099 From Middle Model RCP 8.5')

ggcorrplot(ann_cor, type = "lower",
           lab = TRUE, title = 'Correlation Across Greenscape Metrics \n2000-2099 From Middle Model RCP 8.5',
           colors = c("#6D9EC1", "white", "#E46726"))

ggsave('output/final_plots/correlation/ann_cor_greenscape_85.png', width = 6, height = 4)

herds <- c('AtlanticRimNorth',
           'AtlanticRimSouth',
           'Chokecherry',
           'Cody',
           'DEERP',
           'Dubois',
           'Lander',
           'Meeteetse',
           'Pinedale',
           'PlatteValleyNorth',
           'PlatteValleySouth',
           'Superior',
           'SweetWaterGreenMtn',
           'Teton',
           'WyomingRange')

