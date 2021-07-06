#this code compares future outputs of pirgd from various scenarios against
#contemporary pirgd from the DLC curve

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(ggthemes)

#set wd
setwd('/Volumes/SSD/climate_effects')

#####################################################
###CREATE LANDCOVER RASTER FOR CLIMATE PROJECTIONS###
#####################################################

#first we need to resample landcover
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

#clean up
rm(clim_ras)

###############################################
###GENERATE AVERAGES FOR CLIMATE PROJECTIONS###
###############################################

#allocate output dataframe to store summary results
dat <- data.frame(model = character(), rcp = character(), period = numeric(),
                  lc = character(),
                  min = numeric(), q1 = numeric(), med = numeric(),
                  mean = numeric(), q3 = numeric(), max = numeric())

#allocate output dataframe to store raw results
dat_raw <- data.frame(model = character(), rcp = character(), period = numeric(),
                      lc = character(), value = numeric())

#allocate output dataframe to store annual spatial average
dat_ann <- data.frame(model = character(), rcp = character(), year = numeric(), period = numeric(),
                      lc = character(), value = numeric())

#create file list of pirgd projection rasters
files <- list.files('wy_projections/pirgd_log_gdd_log_rain', 'pirgd_ann',
                         full.names = T)
names <- list.files('wy_projections/pirgd_log_gdd_log_rain', 'pirgd_ann')

#loop through projection rasters
for(i in 1:length(files)){
  
  #create dataframe output variables
  str <- str_split(names[i], "_")
  model_h <- str[[1]][3]
  rcp_h <- str[[1]][4]
  period_h <- str[[1]][5] %>% str_replace(., '.tif', '') %>% as.numeric
  
  #load projection brick
  proj_ann <- brick(files[i])
  
  #calculate per pixel mean across period
  proj <- proj_ann %>% calc(., function(x) mean(x, na.rm = T)) %>% getValues %>% round
  
  #only move forward if have values
  if(sum(is.na(proj)) != 22967){
  
  #fill summary dataframe for 3 lc types
  #shrub
  dat <- rbind(dat,
               data.frame(model = model_h, rcp = rcp_h, period = period_h, lc = 'shrub',
                          min = min(proj[lc_v == 1], na.rm = T),
                          q1 = quantile(proj[lc_v == 1], probs = 0.25, na.rm = T) %>% as.numeric,
                          med = median(proj[lc_v == 1], na.rm = T),
                          mean = mean(proj[lc_v == 1], na.rm = T) %>% round,
                          q3 = quantile(proj[lc_v == 1], probs = 0.75, na.rm = T) %>% as.numeric,
                          max = max(proj[lc_v == 1], na.rm = T)))
  
  #herb
  dat <- rbind(dat,
               data.frame(model = model_h, rcp = rcp_h, period = period_h, lc = 'herb',
                          min = min(proj[lc_v == 2], na.rm = T),
                          q1 = quantile(proj[lc_v == 2], probs = 0.25, na.rm = T) %>% as.numeric,
                          med = median(proj[lc_v == 2], na.rm = T),
                          mean = mean(proj[lc_v == 2], na.rm = T) %>% round,
                          q3 = quantile(proj[lc_v == 2], probs = 0.75, na.rm = T) %>% as.numeric,
                          max = max(proj[lc_v == 2], na.rm = T)))
  
  #evergreen
  dat <- rbind(dat,
               data.frame(model = model_h, rcp = rcp_h, period = period_h, lc = 'evergreen',
                          min = min(proj[lc_v == 3], na.rm = T),
                          q1 = quantile(proj[lc_v == 3], probs = 0.25, na.rm = T) %>% as.numeric,
                          med = median(proj[lc_v == 3], na.rm = T),
                          mean = mean(proj[lc_v == 3], na.rm = T) %>% round,
                          q3 = quantile(proj[lc_v == 3], probs = 0.75, na.rm = T) %>% as.numeric,
                          max = max(proj[lc_v == 3], na.rm = T)))
  
  #create raw dataframes for 3 landcover types
  #shrub
  shrub <- data.frame(value = proj[lc_v == 1]) %>% na.omit
  shrub$model <- model_h
  shrub$rcp <- rcp_h
  shrub$period <- period_h
  shrub$lc <- 'shrub'
  
  #herb
  herb <- data.frame(value = proj[lc_v == 2]) %>% na.omit
  herb$model <- model_h
  herb$rcp <- rcp_h
  herb$period <- period_h
  herb$lc <- 'herb'
  
  #evergreen
  ever <- data.frame(value = proj[lc_v == 3]) %>% na.omit
  ever$model <- model_h
  ever$rcp <- rcp_h
  ever$period <- period_h
  ever$lc <- 'evergreen'
  
  #rbind raw dataframes
  dat_raw <- rbind(dat_raw, shrub, herb, ever)
  
  #clean up
  rm(shrub, herb, ever)
  
  #calculate annual spatial averages for all years
  #shrub
  shrub <- data.frame(value = proj_ann %>% mask(., mask = lc, inverse = T, maskvalue = 1) %>% 
                        cellStats(., stat = 'mean', na.rm = T) %>% round)
  shrub$model <- model_h
  shrub$rcp <- rcp_h
  shrub$year <- period_h + 1:nlayers(proj_ann) - 1
  shrub$period <- period_h
  shrub$lc <- 'shrub'
  
  #herb
  herb <- data.frame(value = proj_ann %>% mask(., mask = lc, inverse = T, maskvalue = 2) %>% 
                       cellStats(., stat = 'mean', na.rm = T) %>% round)
  herb$model <- model_h
  herb$rcp <- rcp_h
  herb$year <- period_h + 1:nlayers(proj_ann) - 1
  herb$period <- period_h
  herb$lc <- 'herb'
  
  #evergreen
  ever <- data.frame(value = proj_ann %>% mask(., mask = lc, inverse = T, maskvalue = 3) %>% 
                       cellStats(., stat = 'mean', na.rm = T) %>% round)
  ever$model <- model_h
  ever$rcp <- rcp_h
  ever$year <- period_h + 1:nlayers(proj_ann) - 1
  ever$period <- period_h
  ever$lc <- 'evergreen'
    
  #rbind to df
  dat_ann <- rbind(dat_ann, shrub, herb, ever)

  #clean up
  rm(shrub, herb, ever)
  rm(proj, str, model_h, rcp_h, period_h, proj_ann)
  }
  
}

#clean up
rm(files, names, lc_v, lc, i)

#############################################
###GENERATE AVERAGES FOR CONTEMPORARY DATA###
#############################################

#load mean pirgd raster and disturbance mask
pirgd <- raster('dlc/mean_maxIRGdate_wy_laea_2001_2018.tif')
mask <- raster('dlc/maxIRGdate_mask_disturbances.tif')

#mask disturbances
pirgd <- mask(pirgd, mask)
rm(mask)

#get pirgd values
pirgd <- getValues(pirgd)

#load lc for pirgd from dlc
lc <- raster('landcover/five_class_landcover_wy_laea_2016.tif')
lc[lc == 3] <- NA
lc[lc == 5] <- NA
lc[lc == 4] <- 3
lc <- getValues(lc)

#fill dataframe for 3 lc types
#shrub
dat <- rbind(dat,
             data.frame(model = 'dlc', rcp = 'na', period = 2001, lc = 'shrub',
                        min = min(pirgd[lc == 1], na.rm = T),
                        q1 = quantile(pirgd[lc == 1], probs = 0.25, na.rm = T) %>% as.numeric,
                        med = median(pirgd[lc == 1], na.rm = T),
                        mean = mean(pirgd[lc == 1], na.rm = T) %>% round,
                        q3 = quantile(pirgd[lc == 1], probs = 0.75, na.rm = T) %>% as.numeric,
                        max = max(pirgd[lc == 1], na.rm = T)))

#herb
dat <- rbind(dat,
             data.frame(model = 'dlc', rcp = 'na', period = 2001, lc = 'herb',
                        min = min(pirgd[lc == 2], na.rm = T),
                        q1 = quantile(pirgd[lc == 2], probs = 0.25, na.rm = T) %>% as.numeric,
                        med = median(pirgd[lc == 2], na.rm = T),
                        mean = mean(pirgd[lc == 2], na.rm = T) %>% round,
                        q3 = quantile(pirgd[lc == 2], probs = 0.75, na.rm = T) %>% as.numeric,
                        max = max(pirgd[lc == 2], na.rm = T)))

#evergreen
dat <- rbind(dat,
             data.frame(model = 'dlc', rcp = 'na', period = 2001, lc = 'evergreen',
                        min = min(pirgd[lc == 3], na.rm = T),
                        q1 = quantile(pirgd[lc == 3], probs = 0.25, na.rm = T) %>% as.numeric,
                        med = median(pirgd[lc == 3], na.rm = T),
                        mean = mean(pirgd[lc == 3], na.rm = T) %>% round,
                        q3 = quantile(pirgd[lc == 3], probs = 0.75, na.rm = T) %>% as.numeric,
                        max = max(pirgd[lc == 3], na.rm = T)))

#create raw dataframes for 3 landcover types
#name rcp 45 and rcp 85 to simplify figures
#shrub
shrub <- data.frame(value = pirgd[lc == 1]) %>% na.omit
shrub$model <- 'DLC'
shrub$rcp <- 'rcp45'
shrub$period <- 2001
shrub$lc <- 'shrub'

#herb
herb <- data.frame(value = pirgd[lc == 2]) %>% na.omit
herb$model <- 'DLC'
herb$rcp <- 'rcp45'
herb$period <- 2001
herb$lc <- 'herb'

#evergreen
ever <- data.frame(value = pirgd[lc == 3]) %>% na.omit
ever$model <- 'DLC'
ever$rcp <- 'rcp45'
ever$period <- 2001
ever$lc <- 'evergreen'

#rbind raw dataframes
dat_raw <- rbind(dat_raw, shrub, herb, ever)

#rcp85
#shrub
shrub <- data.frame(value = pirgd[lc == 1]) %>% na.omit
shrub$model <- 'DLC'
shrub$rcp <- 'rcp85'
shrub$period <- 2001
shrub$lc <- 'shrub'

#herb
herb <- data.frame(value = pirgd[lc == 2]) %>% na.omit
herb$model <- 'DLC'
herb$rcp <- 'rcp85'
herb$period <- 2001
herb$lc <- 'herb'

#evergreen
ever <- data.frame(value = pirgd[lc == 3]) %>% na.omit
ever$model <- 'DLC'
ever$rcp <- 'rcp85'
ever$period <- 2001
ever$lc <- 'evergreen'

#rbind raw dataframes
dat_raw <- rbind(dat_raw, shrub, herb, ever)

#clean up
rm(pirgd, lc, shrub, herb, ever)

#calculate contemporary values for annual spatial average
pirgd_ann <- brick('dlc/maxIRGdate_wy_laea_2000_2019_masked.tif')
lc <- raster('landcover/five_class_landcover_wy_laea_2016.tif')
lc[lc == 3] <- NA
lc[lc == 5] <- NA
lc[lc == 4] <- 3

#create ann dataframes for 3 landcover types
#name rcp 45 and rcp 85 to simplify figures
#shrub
shrub <- data.frame(value = pirgd_ann %>% mask(., mask = lc, inverse = T, maskvalue = 1) %>% 
                      cellStats(., stat = 'mean', na.rm = T) %>% round)
shrub$model <- 'DLC'
shrub$rcp <- 'rcp45'
shrub$year <- 1999 + 1:nlayers(pirgd_ann)
shrub$period <- 2001
shrub$lc <- 'shrub'

#herb
herb <- data.frame(value = pirgd_ann %>% mask(., mask = lc, inverse = T, maskvalue = 2) %>% 
                     cellStats(., stat = 'mean', na.rm = T) %>% round)
herb$model <- 'DLC'
herb$rcp <- 'rcp45'
herb$year <- 1999 + 1:nlayers(pirgd_ann)
herb$period <- 2001
herb$lc <- 'herb'

#evergreen
ever <- data.frame(value = pirgd_ann %>% mask(., mask = lc, inverse = T, maskvalue = 3) %>% 
                     cellStats(., stat = 'mean', na.rm = T) %>% round)
ever$model <- 'DLC'
ever$rcp <- 'rcp45'
ever$year <- 1999 + 1:nlayers(pirgd_ann)
ever$period <- 2001
ever$lc <- 'evergreen'

#rbind raw dataframes
dat_ann <- rbind(dat_ann, shrub, herb, ever)

#rcp85
#shrub
shrub$rcp <- 'rcp85'

#herb
herb$rcp <- 'rcp85'

#evergreen
ever$rcp <- 'rcp85'

#rbind raw dataframes
dat_ann <- rbind(dat_ann, shrub, herb, ever)

#clean up
rm(pirgd_ann, lc, shrub, herb, ever)

####################
###CREATE FIGURES###
####################

#check erroneous values
#dat_raw3 <- dat_raw[dat_raw$value > 365 | dat_raw$value < 1,]

#clean data for erroneous values
#reset to 1 and 365
dat_raw2 <- dat_raw
dat_raw2$value[dat_raw2$value <= 0] <- 1
dat_raw2$value[dat_raw2$value >= 366] <- 365

#change BIT data to different model name
#dat_raw2$model[dat_raw2$period == 1999] <- 'Back in Time (HADGEM2)'

#set models to factor levels
dat_raw2$model <- as.factor(dat_raw2$model)

#set color scale for all 6 models
#need to choose a color palette colorblind friendly
#https://personal.sron.nl/~pault/#sec:qualitative
myColors <- c('#4477AA', '#000000', '#228833', '#CCBB44', '#EE6677', 
              '#AA3377')
names(myColors) <- c("CNRM-CM5", "DLC", "HadGEM2-ES365", "inmcm4", "IPSL-CM5A-MR", "MIROC-ESM-CHEM")
colScale <- scale_color_manual(name = "Model", values = myColors,
                               breaks = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'DLC',
                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'),
                               labels = c('HadGEM2-ES365' = 'HadGEM2-ES365',
                                          'inmcm4' = 'inmcm4',
                                          'CNRM-CM5' = 'CNRM-CM5',
                                          'DLC' = 'Baseline',
                                          'IPSL-CM5A-MR' = 'IPSL-CM5A-MR',
                                          'MIROC-ESM-CHEM' = 'MIROC-ESM-CHEM'))

fillScale <- scale_fill_manual(name = "Model (back in time light colors)", values = myColors,
                               breaks = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'DLC',
                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'),
                               labels = c('HadGEM2-ES365' = 'Middle',
                                          'inmcm4' = 'Low Temp/Low Precip',
                                          'CNRM-CM5' = 'Low Temp/High Precip',
                                          'DLC' = 'Baseline',
                                          'IPSL-CM5A-MR' = 'High Temp/Low Precip',
                                          'MIROC-ESM-CHEM' = 'High Temp/High Precip'))

fillScale2 <- scale_fill_manual(name = "Model (back in time light colors)", values = myColors,
                               breaks = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'DLC',
                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'),
                               labels = c('HadGEM2-ES365' = 'Mid',
                                          'inmcm4' = 'LT/LP',
                                          'CNRM-CM5' = 'LT/HP',
                                          'DLC' = 'Baseline',
                                          'IPSL-CM5A-MR' = 'HT/LP',
                                          'MIROC-ESM-CHEM' = 'HT/HP'))

myColors2 <- c('#4477AAFF', '#4477AA4D', '#000000FF', '#0000004D', 
               '#228833FF', '#2288334D', '#CCBB44FF', '#CCBB444D',
               '#EE6677FF', '#EE66774D', '#AA3377FF', '#AA33774D')

names(myColors2) <- c("CNRM-CM5", "CNRM-CM5 BIT", "DLC", "DLC BIT",
                      "HadGEM2-ES365", "HadGEM2-ES365 BIT",
                      "inmcm4", "inmcm4 BIT",
                      "IPSL-CM5A-MR", "IPSL-CM5A-MR BIT",
                      "MIROC-ESM-CHEM", "MIROC-ESM-CHEM BIT")

fillScale_bit <- scale_fill_manual(name = "Model", values = myColors2,
                                breaks = c('HadGEM2-ES365', "HadGEM2-ES365 BIT", 
                                           'inmcm4', 'inmcm4 BIT',
                                           'CNRM-CM5', 'CNRM-CM5 BIT',
                                           'DLC', 'DLC BIT',
                                           'IPSL-CM5A-MR', 'IPSL-CM5A-MR BIT',
                                           'MIROC-ESM-CHEM', 'MIROC-ESM-CHEM BIT'),
                                labels = c('HadGEM2-ES365' = 'Mid',
                                           'HadGEM2-ES365 BIT' = 'Mid BIT',
                                           'inmcm4' = 'LT/LP',
                                           'inmcm4 BIT' = 'LT/LP BIT',
                                           'CNRM-CM5' = 'LT/HP',
                                           'CNRM-CM5 BIT' = 'LT/HP BIT',
                                           'DLC' = 'Baseline',
                                           'DLC BIT' = 'Baseline BIT',
                                           'IPSL-CM5A-MR' = 'HT/LP',
                                           'IPSL-CM5A-MR BIT' = 'HT/LP BIT',
                                           'MIROC-ESM-CHEM' = 'HT/HP',
                                           'MIROC-ESM-CHEM BIT' = 'HT/HP BIT'))

#set universal theme
theme_e <- theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                 axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

#reorder model factor order
dat_raw2$model <- ordered(dat_raw2$model, levels = c('HadGEM2-ES365', 'inmcm4', 'CNRM-CM5', 'DLC',
                                                     'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'))

#PLOT OVERALL IN STATE
#organize rcp45 and rcp85 data separately based on median value. Include BIT data.
#2040
dat_rcp45_2040 <- dat_raw2[dat_raw2$rcp == 'rcp45' & (dat_raw2$period %in% c(2040, 2001, 1999)),]
#dat_rcp45_2040$model <- fct_reorder(dat_rcp45_2040$model, dat_rcp45_2040$value, .fun = median)
dat_rcp85_2040 <- dat_raw2[dat_raw2$rcp == 'rcp85' & (dat_raw2$period %in% c(2040, 2001, 1999)),]
#dat_rcp85_2040$model <- fct_reorder(dat_rcp85_2040$model, dat_rcp85_2040$value, .fun = median)

#2070
dat_rcp45_2070 <- dat_raw2[dat_raw2$rcp == 'rcp45' & (dat_raw2$period %in% c(2070, 2001, 1999)),]
#dat_rcp45_2070$model <- fct_reorder(dat_rcp45_2070$model, dat_rcp45_2070$value, .fun = median)
dat_rcp85_2070 <- dat_raw2[dat_raw2$rcp == 'rcp85' & (dat_raw2$period %in% c(2070, 2001, 1999)),]
#dat_rcp85_2070$model <- fct_reorder(dat_rcp85_2070$model, dat_rcp85_2070$value, .fun = median)

#PLOT BY LC TYPE
#organize data by rcp and lc type
#2040
dat_rcp45_2040_shrub <- dat_rcp45_2040[dat_rcp45_2040$lc == 'shrub',]
#dat_rcp45_2040_shrub$model <- fct_reorder(dat_rcp45_2040_shrub$model, dat_rcp45_2040_shrub$value, .fun = median)
dat_rcp85_2040_shrub <- dat_rcp85_2040[dat_rcp85_2040$lc == 'shrub',]
#dat_rcp85_2040_shrub$model <- fct_reorder(dat_rcp85_2040_shrub$model, dat_rcp85_2040_shrub$value, .fun = median)

dat_rcp45_2040_herb <- dat_rcp45_2040[dat_rcp45_2040$lc == 'herb',]
#dat_rcp45_2040_herb$model <- fct_reorder(dat_rcp45_2040_herb$model, dat_rcp45_2040_herb$value, .fun = median)
dat_rcp85_2040_herb <- dat_rcp85_2040[dat_rcp85_2040$lc == 'herb',]
#dat_rcp85_2040_herb$model <- fct_reorder(dat_rcp85_2040_herb$model, dat_rcp85_2040_herb$value, .fun = median)

dat_rcp45_2040_ever <- dat_rcp45_2040[dat_rcp45_2040$lc == 'evergreen',]
#dat_rcp45_2040_ever$model <- fct_reorder(dat_rcp45_2040_ever$model, dat_rcp45_2040_ever$value, .fun = median)
dat_rcp85_2040_ever <- dat_rcp85_2040[dat_rcp85_2040$lc == 'evergreen',]
#dat_rcp85_2040_ever$model <- fct_reorder(dat_rcp85_2040_ever$model, dat_rcp85_2040_ever$value, .fun = median)

#2070
dat_rcp45_2070_shrub <- dat_rcp45_2070[dat_rcp45_2070$lc == 'shrub',]
#dat_rcp45_2070_shrub$model <- fct_reorder(dat_rcp45_2070_shrub$model, dat_rcp45_2070_shrub$value, .fun = median)
dat_rcp85_2070_shrub <- dat_rcp85_2070[dat_rcp85_2070$lc == 'shrub',]
#dat_rcp85_2070_shrub$model <- fct_reorder(dat_rcp85_2070_shrub$model, dat_rcp85_2070_shrub$value, .fun = median)

dat_rcp45_2070_herb <- dat_rcp45_2070[dat_rcp45_2070$lc == 'herb',]
#dat_rcp45_2070_herb$model <- fct_reorder(dat_rcp45_2070_herb$model, dat_rcp45_2070_herb$value, .fun = median)
dat_rcp85_2070_herb <- dat_rcp85_2070[dat_rcp85_2070$lc == 'herb',]
#dat_rcp85_2070_herb$model <- fct_reorder(dat_rcp85_2070_herb$model, dat_rcp85_2070_herb$value, .fun = median)

dat_rcp45_2070_ever <- dat_rcp45_2070[dat_rcp45_2070$lc == 'evergreen',]
#dat_rcp45_2070_ever$model <- fct_reorder(dat_rcp45_2070_ever$model, dat_rcp45_2070_ever$value, .fun = median)
dat_rcp85_2070_ever <- dat_rcp85_2070[dat_rcp85_2070$lc == 'evergreen',]
#dat_rcp85_2070_ever$model <- fct_reorder(dat_rcp85_2070_ever$model, dat_rcp85_2070_ever$value, .fun = median)


#create rcp and landcover variable so can plot all together
#2040
dat_rcp45_2040$group <- 'DAOverall_RCP45'
dat_rcp85_2040$group <- 'DBOverall_RCP85'
dat_rcp45_2040_shrub$group <- 'AAShrub_RCP45'
dat_rcp85_2040_shrub$group <- 'ABShrub_RCP85'
dat_rcp45_2040_herb$group <- 'BAHerb_RCP45'
dat_rcp85_2040_herb$group <- 'BBHerb_RCP85'
dat_rcp45_2040_ever$group <- 'CAEvergreen_RCP45'
dat_rcp85_2040_ever$group <- 'CBEvergreen_RCP85'

#2070
dat_rcp45_2070$group <- 'DAOverall_RCP45'
dat_rcp85_2070$group <- 'DBOverall_RCP85'
dat_rcp45_2070_shrub$group <- 'AAShrub_RCP45'
dat_rcp85_2070_shrub$group <- 'ABShrub_RCP85'
dat_rcp45_2070_herb$group <- 'BAHerb_RCP45'
dat_rcp85_2070_herb$group <- 'BBHerb_RCP85'
dat_rcp45_2070_ever$group <- 'CAEvergreen_RCP45'
dat_rcp85_2070_ever$group <- 'CBEvergreen_RCP85'

#create new df with median values and quartiles of DLC to plot as error bar
dat_baseline <- data.frame(group = c('DAOverall_RCP45', 'DBOverall_RCP85',
                                     'AAShrub_RCP45', 'ABShrub_RCP85',
                                     'BAHerb_RCP45', 'BBHerb_RCP85',
                                     'CAEvergreen_RCP45', 'CBEvergreen_RCP85'),
                           med_line = c(median(dat_rcp45_2040$value[dat_rcp45_2040$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040$value[dat_rcp85_2040$model == 'DLC'], na.rm = T),
                                        median(dat_rcp45_2040_shrub$value[dat_rcp45_2040_shrub$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040_shrub$value[dat_rcp85_2040_shrub$model == 'DLC'], na.rm = T),
                                        median(dat_rcp45_2040_herb$value[dat_rcp45_2040_herb$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040_herb$value[dat_rcp85_2040_herb$model == 'DLC'], na.rm = T),
                                        median(dat_rcp45_2040_ever$value[dat_rcp45_2040_ever$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040_ever$value[dat_rcp85_2040_ever$model == 'DLC'], na.rm = T)),
                           q1 = c(quantile(dat_rcp45_2040$value[dat_rcp45_2040$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040$value[dat_rcp85_2040$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp45_2040_shrub$value[dat_rcp45_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040_shrub$value[dat_rcp85_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp45_2040_herb$value[dat_rcp45_2040_herb$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040_herb$value[dat_rcp85_2040_herb$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp45_2040_ever$value[dat_rcp45_2040_ever$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040_ever$value[dat_rcp85_2040_ever$model == 'DLC'], na.rm = T, probs = 0.25)),
                           q3 = c(quantile(dat_rcp45_2040$value[dat_rcp45_2040$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040$value[dat_rcp85_2040$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp45_2040_shrub$value[dat_rcp45_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040_shrub$value[dat_rcp85_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp45_2040_herb$value[dat_rcp45_2040_herb$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040_herb$value[dat_rcp85_2040_herb$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp45_2040_ever$value[dat_rcp45_2040_ever$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040_ever$value[dat_rcp85_2040_ever$model == 'DLC'], na.rm = T, probs = 0.75)))

#separate first and last category because if issues with line size
dat_baseline2 <- dat_baseline[dat_baseline$group %in% c('AAShrub_RCP45', 'DBOverall_RCP85'),]
dat_baseline <- dat_baseline[!(dat_baseline$group %in% c('AAShrub_RCP45', 'DBOverall_RCP85')),]

#set linetype aesthetic
line_types <- c("Median" = 'dashed',"Q1/Q3" = 'dotted')

# #calculate number of outliers
# #pre allocate outlier size
# out <- 0
# 
# #loop through all values
# for(g in unique(dat_rcp45_2040_herb$model[dat_rcp45_2040_herb$model != 'DLC'])){
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp45_2040_herb$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp85_2040_herb$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp45_2040_shrub$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp85_2040_shrub$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp45_2040_ever$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp85_2040_ever$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp45_2040$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2040_herb$value[dat_rcp85_2040$model == g])$out)
# }
# 
# #calc total number of vals
# tot <- NROW(dat_rcp45_2040_herb[dat_rcp45_2040_herb$model != 'DLC',]) +
#   NROW(dat_rcp45_2040_herb[dat_rcp85_2040_herb$model != 'DLC',]) +
#   NROW(dat_rcp45_2040_herb[dat_rcp45_2040_shrub$model != 'DLC',]) +
#   NROW(dat_rcp45_2040_herb[dat_rcp85_2040_shrub$model != 'DLC',]) +
#   NROW(dat_rcp45_2040_herb[dat_rcp45_2040_ever$model != 'DLC',]) +
#   NROW(dat_rcp45_2040_herb[dat_rcp85_2040_ever$model != 'DLC',]) +
#   NROW(dat_rcp45_2040_herb[dat_rcp45_2040$model != 'DLC',]) +
#   NROW(dat_rcp45_2040_herb[dat_rcp85_2040$model != 'DLC',])
# 
# #calc percent of outliers
# out / tot * 100

#set universal theme
theme_out <- theme(plot.title = element_text(size = 20, hjust = 0.5), axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                   axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                   legend.position = '',
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#calculate boxplot statistics and remove duplicates
dat_rcp45_2040_herb_stats <- dat_rcp45_2040_herb %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2040_herb_stats <- dat_rcp85_2040_herb %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp45_2040_shrub_stats <- dat_rcp45_2040_shrub %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2040_shrub_stats <- dat_rcp85_2040_shrub %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp45_2040_ever_stats <- dat_rcp45_2040_ever %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2040_ever_stats <- dat_rcp85_2040_ever %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp45_2040_stats <- dat_rcp45_2040 %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2040_stats <- dat_rcp85_2040 %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

#set alpha aesthetic
alpha_types <- c('1999' = 0.4, '2040' = 1, '2070' = 1)

#create boxplot with hlines
#2040
ggplot() +
  geom_boxplot(data = dat_rcp45_2040_herb_stats[dat_rcp45_2040_herb_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_herb_stats[dat_rcp85_2040_herb_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2040_shrub_stats[dat_rcp45_2040_shrub_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_shrub_stats[dat_rcp85_2040_shrub_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2040_ever_stats[dat_rcp45_2040_ever_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_ever_stats[dat_rcp85_2040_ever_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2040_stats[dat_rcp45_2040_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_stats[dat_rcp85_2040_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  fillScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (Date)') +
  ggtitle('Mid-Century (2040-2069)') + theme_bw() + theme_out +
  scale_y_continuous(breaks = c(91, 121, 152), labels = c('Apr 1', 'May 1', 'Jun 1'), limits = c(85, 166)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5) +
  geom_errorbar(data = dat_baseline, width = 1, size = 0.75, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, size = 0.75, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, size = 0.75, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, size = 0.75, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, size = 0.75, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, size = 0.75, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  scale_linetype_manual(name = 'Baseline Values', values = line_types, labels = c('Q1/Q3' = '1st/3rd Quartile')) +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                      '2040' = 'Projection')) +
  guides(alpha = guide_legend(order = 1))

ggsave('output/final_plots/pirgd_comparisons_2040_ci_bit.png',
       width = 9, height = 5, unit = "in")

#create boxplot with hlines
#2040
ggplot() +
  geom_boxplot(data = dat_rcp45_2040_herb_stats[dat_rcp45_2040_herb_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_herb_stats[dat_rcp85_2040_herb_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2040_shrub_stats[dat_rcp45_2040_shrub_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_shrub_stats[dat_rcp85_2040_shrub_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2040_ever_stats[dat_rcp45_2040_ever_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_ever_stats[dat_rcp85_2040_ever_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2040_stats[dat_rcp45_2040_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2040_stats[dat_rcp85_2040_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, color = factor(period)), stat = 'identity') +
  fillScale2 + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (Date)') +
  ggtitle('Mid-Century (2040-2069)') + theme_bw() + theme_e +
  scale_y_continuous(breaks = c(91, 121, 152), labels = c('Apr 1', 'May 1', 'Jun 1'), limits = c(85, 166)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5) +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  scale_linetype_manual(name = 'Baseline Values', values = line_types)

# #calculate number of outliers
# #pre allocate outlier size
# out <- 0
# 
# #loop through all values
# for(g in unique(dat_rcp45_2070_herb$model[dat_rcp45_2070_herb$model != 'DLC'])){
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp45_2070_herb$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp85_2070_herb$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp45_2070_shrub$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp85_2070_shrub$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp45_2070_ever$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp85_2070_ever$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp45_2070$model == g])$out)
#   out <- out + length(boxplot.stats(dat_rcp45_2070_herb$value[dat_rcp85_2070$model == g])$out)
# }
# 
# #calc total number of vals
# tot <- NROW(dat_rcp45_2070_herb[dat_rcp45_2070_herb$model != 'DLC',]) +
#   NROW(dat_rcp45_2070_herb[dat_rcp85_2070_herb$model != 'DLC',]) +
#   NROW(dat_rcp45_2070_herb[dat_rcp45_2070_shrub$model != 'DLC',]) +
#   NROW(dat_rcp45_2070_herb[dat_rcp85_2070_shrub$model != 'DLC',]) +
#   NROW(dat_rcp45_2070_herb[dat_rcp45_2070_ever$model != 'DLC',]) +
#   NROW(dat_rcp45_2070_herb[dat_rcp85_2070_ever$model != 'DLC',]) +
#   NROW(dat_rcp45_2070_herb[dat_rcp45_2070$model != 'DLC',]) +
#   NROW(dat_rcp45_2070_herb[dat_rcp85_2070$model != 'DLC',])
# 
# #calc percent of outliers
# out / tot * 100

#calculate boxplot statistics and remove duplicates
dat_rcp45_2070_herb_stats <- dat_rcp45_2070_herb %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2070_herb_stats <- dat_rcp85_2070_herb %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp45_2070_shrub_stats <- dat_rcp45_2070_shrub %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2070_shrub_stats <- dat_rcp85_2070_shrub %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp45_2070_ever_stats <- dat_rcp45_2070_ever %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2070_ever_stats <- dat_rcp85_2070_ever %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp45_2070_stats <- dat_rcp45_2070 %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

dat_rcp85_2070_stats <- dat_rcp85_2070 %>% group_by(model, period) %>%
  mutate(y05 = quantile(value, 0.05), y25 = quantile(value, 0.25),
         y50 = median(value), y75 = quantile(value, 0.75),
         y95 = quantile(value, 0.95), model_period = str_c(model, period)) %>%
  .[!duplicated(.$model_period),] %>%
  select(-'value')

#2070
ggplot() +
  geom_boxplot(data = dat_rcp45_2070_herb_stats[dat_rcp45_2070_herb_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2070_herb_stats[dat_rcp85_2070_herb_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2070_shrub_stats[dat_rcp45_2070_shrub_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2070_shrub_stats[dat_rcp85_2070_shrub_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2070_ever_stats[dat_rcp45_2070_ever_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2070_ever_stats[dat_rcp85_2070_ever_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp45_2070_stats[dat_rcp45_2070_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  geom_boxplot(data = dat_rcp85_2070_stats[dat_rcp85_2070_stats$model != 'DLC',], 
               aes(x = group, fill = model, ymin = y05, lower = y25,
                   middle = y50, upper = y75, ymax = y95, alpha = factor(period)), stat = 'identity') +
  fillScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (Date)') +
  ggtitle('End of Century (2070-2099)') + theme_bw() + theme_out +
  scale_y_continuous(breaks = c(91, 121, 152), labels = c('Apr 1', 'May 1', 'Jun 1'), limits = c(85, 166)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5) +
  geom_errorbar(data = dat_baseline, width = 1, size = 0.75, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, size = 0.75, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, size = 0.75, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, size = 0.75, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, size = 0.75, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, size = 0.75, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  scale_linetype_manual(name = 'Baseline Values', values = line_types) +
  scale_alpha_manual(name = 'Comparison', values = alpha_types, labels = c('1999' = 'Back in Time (transparent)',
                                                                           '2070' = 'Projection'))

ggsave('output/final_plots/pirgd_comparisons_2070_ci_bit.png',
       width = 9, height = 5, unit = "in")

#create boxplot with baseline boxes

#2040
ggplot() +
  geom_boxplot(data = dat_rcp45_2040_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  colScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (DOY)') +
  ggtitle('PIRGd Projections for Mid-Century in WY (2040-2069)') + theme_bw() + theme_e +
  scale_y_continuous(limits = c(0, 365)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5)

ggsave('output/final_plots/pirgd_comparisons_2040_basebox.png',
       width = 12, height = 8, unit = "in")

#2070
ggplot() +
  geom_boxplot(data = dat_rcp45_2070_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  colScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (DOY)') +
  ggtitle('PIRGd Projections for Mid-Century in WY (2070-2099)') + theme_bw() + theme_e +
  scale_y_continuous(limits = c(0, 365)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5)

ggsave('output/final_plots/pirgd_comparisons_2070_basebox.png',
       width = 12, height = 8, unit = "in")

#################################################
###COMPARE ANNUAL PIRGd ACROSS MIDDLE SCENARIO###
#################################################

#load mean pirgd raster and disturbance mask
pirgd_ann <- stack('dlc/maxIRGdate_wy_laea_2000_2019_masked.tif')

#create data frame for annual values
dat_sum <- data.frame(Year = 2000:2019, RCP = 'Baseline',
                      PIRGd = round(cellStats(pirgd_ann, stat = 'mean', na.rm = T)))

#clean up
rm(pirgd_ann)

#load future data
pirgd_45_1999 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_1999.tif')
pirgd_85_1999 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_1999.tif')
pirgd_45_2020 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_2020.tif')
pirgd_85_2020 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_2020.tif')
pirgd_45_2040 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_2040.tif')
pirgd_85_2040 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_2040.tif')
pirgd_45_2070 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_2070.tif')
pirgd_85_2070 <- stack('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_2070.tif')

#set min and max values 1 and 365
pirgd_45_1999[pirgd_45_1999 < 1] <- 1
pirgd_45_1999[pirgd_45_1999 > 365] <- 365
pirgd_85_1999[pirgd_85_1999 < 1] <- 1
pirgd_85_1999[pirgd_85_1999 > 365] <- 365
pirgd_45_2020[pirgd_45_2020 < 1] <- 1
pirgd_45_2020[pirgd_45_2020 > 365] <- 365
pirgd_85_2020[pirgd_85_2020 < 1] <- 1
pirgd_85_2020[pirgd_85_2020 > 365] <- 365
pirgd_45_2040[pirgd_45_2040 < 1] <- 1
pirgd_45_2040[pirgd_45_2040 > 365] <- 365
pirgd_85_2040[pirgd_85_2040 < 1] <- 1
pirgd_85_2040[pirgd_85_2040 > 365] <- 365
pirgd_45_2070[pirgd_45_2070 < 1] <- 1
pirgd_45_2070[pirgd_45_2070 > 365] <- 365
pirgd_85_2070[pirgd_85_2070 < 1] <- 1
pirgd_85_2070[pirgd_85_2070 > 365] <- 365

#rbind to df
dat_sum <- rbind(dat_sum,
                 data.frame(Year = 2000:2019, RCP = 'BIT RCP 4.5',
                            PIRGd = round(cellStats(pirgd_45_1999[[2:21]], stat = 'mean', na.rm = T))),
                 data.frame(Year = 2000:2019, RCP = 'BIT RCP 8.5',
                            PIRGd = round(cellStats(pirgd_85_1999[[2:21]], stat = 'mean', na.rm = T))),
                 data.frame(Year = 2020:2039, RCP = 'RCP 4.5',
                            PIRGd = round(cellStats(pirgd_45_2020, stat = 'mean', na.rm = T))),
                 data.frame(Year = 2020:2039, RCP = 'RCP 8.5',
                            PIRGd = round(cellStats(pirgd_85_2020, stat = 'mean', na.rm = T))),
                 data.frame(Year = 2040:2069, RCP = 'RCP 4.5',
                            PIRGd = round(cellStats(pirgd_45_2040, stat = 'mean', na.rm = T))),
                 data.frame(Year = 2040:2069, RCP = 'RCP 8.5',
                            PIRGd = round(cellStats(pirgd_85_2040, stat = 'mean', na.rm = T))),
                 data.frame(Year = 2070:2099, RCP = 'RCP 4.5',
                            PIRGd = round(cellStats(pirgd_45_2070, stat = 'mean', na.rm = T))),
                 data.frame(Year = 2070:2099, RCP = 'RCP 8.5',
                            PIRGd = round(cellStats(pirgd_85_2070, stat = 'mean', na.rm = T))))

#set plot color scheme
ann_col <- c('Baseline' = '#228833', 
             'BIT RCP 4.5' = '#CCBB44', 
             'BIT RCP 8.5' = '#AA3377',
             'RCP 4.5' = '#4477AA',
             'RCP 8.5' = '#EE6677')

# #plot pirgd time series
# ggplot(data = dat_sum, aes(x = Year, y = PIRGd, color = RCP)) +
#   geom_point() +
#   theme_bw() + theme_e +
#   scale_color_manual(values = ann_col, name = 'RCP') + 
#   xlab('Year') + ylab('PIRGd') + ggtitle('Average PIRGd Middle Scenario (2000-2099)') +
#   geom_vline(xintercept = 2020)
# 
# ggsave('output/final_plots/pirgd_ann_middle_scenario.png',
#        width = 12, height = 8, unit = "in")

#plot ss time series
ggplot(data = dat_sum, aes(x = Year, y = PIRGd, color = RCP)) +
  #  geom_line(aes(color = rcp)) +
  geom_point() +
  theme_bw() + theme_e +
  scale_color_manual(values = ann_col,  name = 'RCP', labels = c('BIT RCP 4.5' = 'BIT 4.5', 'BIT RCP 8.5' = 'BIT 8.5')) +
  geom_smooth(data = subset(dat_sum, RCP == 'RCP 4.5' | RCP == 'RCP 8.5'),
              aes(x = Year, y = PIRGd, color = RCP, linetype = 'Trend'), method = 'lm', se = F) + 
  xlab('Year') + ylab('PIRGd (Days)') + ggtitle('Average PIRGd Middle Model (2000-2099)') +
  geom_vline(xintercept = 2020) +
  scale_y_continuous(breaks = c(105, 121, 135), labels = c('Apr 15', 'May 1', 'May 15'), limits = c(98, 145)) +
  geom_segment(aes(x = 2000, xend = 2020, y = mean(dat_sum$PIRGd[dat_sum$RCP == 'Baseline']),
                   yend = mean(dat_sum$PIRGd[dat_sum$RCP == 'Baseline']),
                   linetype = 'Mean',  color = 'Baseline'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(dat_sum$PIRGd[dat_sum$RCP == 'BIT RCP 4.5']),
                   yend = mean(dat_sum$PIRGd[dat_sum$RCP == 'BIT RCP 4.5']),
                   linetype = 'Mean', color = 'BIT RCP 4.5'), lwd = 0.8) + 
  geom_segment(aes(x = 2000, xend = 2020, y = mean(dat_sum$PIRGd[dat_sum$RCP == 'BIT RCP 8.5']),
                   yend = mean(dat_sum$PIRGd[dat_sum$RCP == 'BIT RCP 8.5']),
                   linetype = 'Mean',  color = 'BIT RCP 8.5'), lwd = 0.8) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2, keywidth = 1.8)) +
  scale_linetype_manual('Line Type', values = c('Trend' = 1, 'Mean' = 2))

ggsave('output/final_plots/pirgd_ann_middle_scenario.png', width = 8, height = 4.5)

#########################################
###CREATE FIGURES WITH SPATIAL AVERAGE###
#########################################

#check erroneous values
#dat_raw3 <- dat_raw[dat_raw$value > 365 | dat_raw$value < 1,]

#clean data for erroneous values
#reset to 1 and 365
dat_raw2 <- dat_ann
dat_raw2$value[dat_raw2$value <= 0] <- 1
dat_raw2$value[dat_raw2$value >= 366] <- 365

#change BIT data to different model name
dat_raw2$model[dat_raw2$period == 1999] <- 'Back in Time (HADGEM2)'

#set models to factor levels
dat_raw2$model <- as.factor(dat_raw2$model)

#set color scale for all 6 models
#need to choose a color palette colorblind friendly
#https://personal.sron.nl/~pault/#sec:qualitative
myColors <- c('#BBBBBB', '#4477AA', '#000000', '#228833', '#CCBB44', '#EE6677', 
              '#AA3377')
names(myColors) <- levels(dat_raw2$model)
colScale <- scale_color_manual(name = "Model", values = myColors,
                               breaks = c('Back in Time (HADGEM2)', 'DLC', 'CNRM-CM5', 'HadGEM2-ES365', 'inmcm4',
                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'),
                               labels = c('Back in Time (HADGEM2)' = 'Back in Time (HADGEM2)',
                                          'DLC' = 'Baseline',
                                          'CNRM-CM5' = 'CNRM-CM5',
                                          'HadGEM2-ES365' = 'HadGEM2-ES365',
                                          'inmcm4' = 'inmcm4',
                                          'IPSL-CM5A-MR' = 'IPSL-CM5A-MR',
                                          'MIROC-ESM-CHEM' = 'MIROC-ESM-CHEM'))
fillScale <- scale_fill_manual(name = "Model", values = myColors)

#set universal theme
theme_e <- theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                 axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

#PLOT OVERALL IN STATE
#organize rcp45 and rcp85 data separately based on median value. Include BIT data.
#2040
dat_rcp45_2040 <- dat_raw2[dat_raw2$rcp == 'rcp45' & (dat_raw2$period %in% c(2040, 2001, 1999)),]
dat_rcp45_2040$model <- fct_reorder(dat_rcp45_2040$model, dat_rcp45_2040$value, .fun = median, na.rm = TRUE)
dat_rcp85_2040 <- dat_raw2[dat_raw2$rcp == 'rcp85' & (dat_raw2$period %in% c(2040, 2001, 1999)),]
dat_rcp85_2040$model <- fct_reorder(dat_rcp85_2040$model, dat_rcp85_2040$value, .fun = median, na.rm = TRUE)

#2070
dat_rcp45_2070 <- dat_raw2[dat_raw2$rcp == 'rcp45' & (dat_raw2$period %in% c(2070, 2001, 1999)),]
dat_rcp45_2070$model <- fct_reorder(dat_rcp45_2070$model, dat_rcp45_2070$value, .fun = median, na.rm = TRUE)
dat_rcp85_2070 <- dat_raw2[dat_raw2$rcp == 'rcp85' & (dat_raw2$period %in% c(2070, 2001, 1999)),]
dat_rcp85_2070$model <- fct_reorder(dat_rcp85_2070$model, dat_rcp85_2070$value, .fun = median, na.rm = TRUE)

#PLOT BY LC TYPE
#organize data by rcp and lc type
#2040
dat_rcp45_2040_shrub <- dat_rcp45_2040[dat_rcp45_2040$lc == 'shrub',]
dat_rcp45_2040_shrub$model <- fct_reorder(dat_rcp45_2040_shrub$model, dat_rcp45_2040_shrub$value, .fun = median, na.rm = TRUE)
dat_rcp85_2040_shrub <- dat_rcp85_2040[dat_rcp85_2040$lc == 'shrub',]
dat_rcp85_2040_shrub$model <- fct_reorder(dat_rcp85_2040_shrub$model, dat_rcp85_2040_shrub$value, .fun = median, na.rm = TRUE)

dat_rcp45_2040_herb <- dat_rcp45_2040[dat_rcp45_2040$lc == 'herb',]
dat_rcp45_2040_herb$model <- fct_reorder(dat_rcp45_2040_herb$model, dat_rcp45_2040_herb$value, .fun = median, na.rm = TRUE)
dat_rcp85_2040_herb <- dat_rcp85_2040[dat_rcp85_2040$lc == 'herb',]
dat_rcp85_2040_herb$model <- fct_reorder(dat_rcp85_2040_herb$model, dat_rcp85_2040_herb$value, .fun = median, na.rm = TRUE)

dat_rcp45_2040_ever <- dat_rcp45_2040[dat_rcp45_2040$lc == 'evergreen',]
dat_rcp45_2040_ever$model <- fct_reorder(dat_rcp45_2040_ever$model, dat_rcp45_2040_ever$value, .fun = median, na.rm = TRUE)
dat_rcp85_2040_ever <- dat_rcp85_2040[dat_rcp85_2040$lc == 'evergreen',]
dat_rcp85_2040_ever$model <- fct_reorder(dat_rcp85_2040_ever$model, dat_rcp85_2040_ever$value, .fun = median, na.rm = TRUE)

#2070
dat_rcp45_2070_shrub <- dat_rcp45_2070[dat_rcp45_2070$lc == 'shrub',]
dat_rcp45_2070_shrub$model <- fct_reorder(dat_rcp45_2070_shrub$model, dat_rcp45_2070_shrub$value, .fun = median, na.rm = TRUE)
dat_rcp85_2070_shrub <- dat_rcp85_2070[dat_rcp85_2070$lc == 'shrub',]
dat_rcp85_2070_shrub$model <- fct_reorder(dat_rcp85_2070_shrub$model, dat_rcp85_2070_shrub$value, .fun = median, na.rm = TRUE)

dat_rcp45_2070_herb <- dat_rcp45_2070[dat_rcp45_2070$lc == 'herb',]
dat_rcp45_2070_herb$model <- fct_reorder(dat_rcp45_2070_herb$model, dat_rcp45_2070_herb$value, .fun = median, na.rm = TRUE)
dat_rcp85_2070_herb <- dat_rcp85_2070[dat_rcp85_2070$lc == 'herb',]
dat_rcp85_2070_herb$model <- fct_reorder(dat_rcp85_2070_herb$model, dat_rcp85_2070_herb$value, .fun = median, na.rm = TRUE)

dat_rcp45_2070_ever <- dat_rcp45_2070[dat_rcp45_2070$lc == 'evergreen',]
dat_rcp45_2070_ever$model <- fct_reorder(dat_rcp45_2070_ever$model, dat_rcp45_2070_ever$value, .fun = median, na.rm = TRUE)
dat_rcp85_2070_ever <- dat_rcp85_2070[dat_rcp85_2070$lc == 'evergreen',]
dat_rcp85_2070_ever$model <- fct_reorder(dat_rcp85_2070_ever$model, dat_rcp85_2070_ever$value, .fun = median, na.rm = TRUE)


#create rcp and landcover variable so can plot all together
#2040
dat_rcp45_2040$group <- 'DAOverall_RCP45'
dat_rcp85_2040$group <- 'DBOverall_RCP85'
dat_rcp45_2040_shrub$group <- 'AAShrub_RCP45'
dat_rcp85_2040_shrub$group <- 'ABShrub_RCP85'
dat_rcp45_2040_herb$group <- 'BAHerb_RCP45'
dat_rcp85_2040_herb$group <- 'BBHerb_RCP85'
dat_rcp45_2040_ever$group <- 'CAEvergreen_RCP45'
dat_rcp85_2040_ever$group <- 'CBEvergreen_RCP85'

#2070
dat_rcp45_2070$group <- 'DAOverall_RCP45'
dat_rcp85_2070$group <- 'DBOverall_RCP85'
dat_rcp45_2070_shrub$group <- 'AAShrub_RCP45'
dat_rcp85_2070_shrub$group <- 'ABShrub_RCP85'
dat_rcp45_2070_herb$group <- 'BAHerb_RCP45'
dat_rcp85_2070_herb$group <- 'BBHerb_RCP85'
dat_rcp45_2070_ever$group <- 'CAEvergreen_RCP45'
dat_rcp85_2070_ever$group <- 'CBEvergreen_RCP85'

#create new df with median values and quartiles of DLC to plot as error bar
dat_baseline <- data.frame(group = c('DAOverall_RCP45', 'DBOverall_RCP85',
                                     'AAShrub_RCP45', 'ABShrub_RCP85',
                                     'BAHerb_RCP45', 'BBHerb_RCP85',
                                     'CAEvergreen_RCP45', 'CBEvergreen_RCP85'),
                           med_line = c(median(dat_rcp45_2040$value[dat_rcp45_2040$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040$value[dat_rcp85_2040$model == 'DLC'], na.rm = T),
                                        median(dat_rcp45_2040_shrub$value[dat_rcp45_2040_shrub$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040_shrub$value[dat_rcp85_2040_shrub$model == 'DLC'], na.rm = T),
                                        median(dat_rcp45_2040_herb$value[dat_rcp45_2040_herb$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040_herb$value[dat_rcp85_2040_herb$model == 'DLC'], na.rm = T),
                                        median(dat_rcp45_2040_ever$value[dat_rcp45_2040_ever$model == 'DLC'], na.rm = T),
                                        median(dat_rcp85_2040_ever$value[dat_rcp85_2040_ever$model == 'DLC'], na.rm = T)),
                           q1 = c(quantile(dat_rcp45_2040$value[dat_rcp45_2040$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040$value[dat_rcp85_2040$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp45_2040_shrub$value[dat_rcp45_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040_shrub$value[dat_rcp85_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp45_2040_herb$value[dat_rcp45_2040_herb$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040_herb$value[dat_rcp85_2040_herb$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp45_2040_ever$value[dat_rcp45_2040_ever$model == 'DLC'], na.rm = T, probs = 0.25),
                                  quantile(dat_rcp85_2040_ever$value[dat_rcp85_2040_ever$model == 'DLC'], na.rm = T, probs = 0.25)),
                           q3 = c(quantile(dat_rcp45_2040$value[dat_rcp45_2040$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040$value[dat_rcp85_2040$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp45_2040_shrub$value[dat_rcp45_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040_shrub$value[dat_rcp85_2040_shrub$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp45_2040_herb$value[dat_rcp45_2040_herb$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040_herb$value[dat_rcp85_2040_herb$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp45_2040_ever$value[dat_rcp45_2040_ever$model == 'DLC'], na.rm = T, probs = 0.75),
                                  quantile(dat_rcp85_2040_ever$value[dat_rcp85_2040_ever$model == 'DLC'], na.rm = T, probs = 0.75)))

#separate first and last category because if issues with line size
dat_baseline2 <- dat_baseline[dat_baseline$group %in% c('AAShrub_RCP45', 'DBOverall_RCP85'),]
dat_baseline <- dat_baseline[!(dat_baseline$group %in% c('AAShrub_RCP45', 'DBOverall_RCP85')),]

#set linetype aesthetic
line_types <- c("Median" = 'dashed',"Q1/Q3" = 'dotted')

#create boxplot with hlines
#2040
ggplot() +
  geom_boxplot(data = dat_rcp45_2040_herb[dat_rcp45_2040_herb$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_herb[dat_rcp85_2040_herb$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040_shrub[dat_rcp45_2040_shrub$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_shrub[dat_rcp85_2040_shrub$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040_ever[dat_rcp45_2040_ever$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_ever[dat_rcp85_2040_ever$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040[dat_rcp45_2040$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040[dat_rcp85_2040$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  colScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (DOY)') +
  ggtitle('PIRGd Projections for Mid-Century in WY (2040-2069)') + theme_bw() + theme_e +
  scale_y_continuous(limits = c(0, 365)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5) +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  scale_linetype_manual(name = 'Baseline Values', values = line_types)

ggsave('output/final_plots/pirgd_comparisons_2040_spatial_avg.png',
       width = 12, height = 8, unit = "in")

#2070
ggplot() +
  geom_boxplot(data = dat_rcp45_2070_herb[dat_rcp45_2070_herb$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_herb[dat_rcp85_2070_herb$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070_shrub[dat_rcp45_2070_shrub$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_shrub[dat_rcp85_2070_shrub$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070_ever[dat_rcp45_2070_ever$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_ever[dat_rcp85_2070_ever$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070[dat_rcp45_2070$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070[dat_rcp85_2070$model != 'DLC',], 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  colScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (DOY)') +
  ggtitle('PIRGd Projections for End of Century in WY (2070-2099)') + theme_bw() + theme_e +
  scale_y_continuous(limits = c(0, 365)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5) +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = med_line, ymin = med_line, linetype = 'Median'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q1, ymin = q1, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q3, ymin = q3, linetype = 'Q1/Q3'), 
                colour="#000000") +
  scale_linetype_manual(name = 'Baseline Values', values = line_types)

ggsave('output/final_plots/pirgd_comparisons_2070_spatial_avg.png',
       width = 12, height = 8, unit = "in")

#create boxplot with baseline boxes

#2040
ggplot() +
  geom_boxplot(data = dat_rcp45_2040_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2040, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2040, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  colScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (DOY)') +
  ggtitle('PIRGd Projections for Mid-Century in WY (2040-2069)') + theme_bw() + theme_e +
  scale_y_continuous(limits = c(0, 365)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5)

ggsave('output/final_plots/pirgd_comparisons_2040_basebox_spatial_avg.png',
       width = 12, height = 8, unit = "in")

#2070
ggplot() +
  geom_boxplot(data = dat_rcp45_2070_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_2070, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_2070, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  colScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (DOY)') +
  ggtitle('PIRGd Projections for Mid-Century in WY (2070-2099)') + theme_bw() + theme_e +
  scale_y_continuous(limits = c(0, 365)) + 
  scale_x_discrete(labels = c('BAHerb_RCP45' = 'Herb \nRCP 4.5', 'BBHerb_RCP85' = 'Herb \nRCP 8.5',
                              'AAShrub_RCP45' = 'Shrub \nRCP 4.5', 'ABShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5)

ggsave('output/final_plots/pirgd_comparisons_2070_basebox_spatial_avg.png',
       width = 12, height = 8, unit = "in")

###########################################################################
###COMPARE ANNUAL PIRGd ACROSS MIDDLE SCENARIO WITH CONFIDENCE INTERVALS###
###########################################################################

#load mean pirgd raster and disturbance mask
pirgd_ann <- brick('dlc/maxIRGdate_wy_laea_2000_2019_masked.tif')

#allocate output dataframe to store annual summary results
dat_ann_sum <- data.frame(rcp = character(), year = numeric(),
                  min = numeric(), q1 = numeric(), med = numeric(),
                  mean = numeric(), q3 = numeric(), max = numeric())

#fill summary dataframe not doing by LC
dat_ann_sum <- rbind(dat_ann_sum,
             data.frame(rcp = 'Baseline', year = 2000:2019,
                        min = cellStats(pirgd_ann, stat = 'min', na.rm = T),
                        q1 = quantile(pirgd_ann, probs = 0.25, na.rm = T) %>% round,
                        med = cellStats(pirgd_ann, stat = 'median', na.rm = T) %>% round,
                        mean = cellStats(pirgd_ann, stat = 'mean', na.rm = T) %>% round,
                        q3 = quantile(pirgd_ann, probs = 0.75, na.rm = T) %>% round,
                        max = cellStats(pirgd_ann, stat = 'max', na.rm = T)))

#clean up
rm(pirgd_ann)

#future data already loaded
#rbind to df
dat_ann_sum <- rbind(dat_ann_sum,
                     data.frame(rcp = 'BIT RCP 4.5', year = 2000:2019,
                                min = cellStats(pirgd_45_1999[[2:21]], stat = 'min', na.rm = T),
                                q1 = quantile(pirgd_45_1999[[2:21]], probs = 0.25, na.rm = T) %>% round,
                                med = cellStats(pirgd_45_1999[[2:21]], stat = 'median', na.rm = T) %>% round,
                                mean = cellStats(pirgd_45_1999[[2:21]], stat = 'mean', na.rm = T) %>% round,
                                q3 = quantile(pirgd_45_1999[[2:21]], probs = 0.75, na.rm = T) %>% round,
                                max = cellStats(pirgd_45_1999[[2:21]], stat = 'max', na.rm = T)),
                     data.frame(rcp = 'RCP 4.5', year = 2020:2039,
                                min = cellStats(pirgd_45_2020, stat = 'min', na.rm = T),
                                q1 = quantile(pirgd_45_2020, probs = 0.25, na.rm = T) %>% round,
                                med = cellStats(pirgd_45_2020, stat = 'median', na.rm = T) %>% round,
                                mean = cellStats(pirgd_45_2020, stat = 'mean', na.rm = T) %>% round,
                                q3 = quantile(pirgd_45_2020, probs = 0.75, na.rm = T) %>% round,
                                max = cellStats(pirgd_45_2020, stat = 'max', na.rm = T)),
                     data.frame(rcp = 'RCP 4.5', year = 2040:2069,
                                min = cellStats(pirgd_45_2040, stat = 'min', na.rm = T),
                                q1 = quantile(pirgd_45_2040, probs = 0.25, na.rm = T) %>% round,
                                med = cellStats(pirgd_45_2040, stat = 'median', na.rm = T) %>% round,
                                mean = cellStats(pirgd_45_2040, stat = 'mean', na.rm = T) %>% round,
                                q3 = quantile(pirgd_45_2040, probs = 0.75, na.rm = T) %>% round,
                                max = cellStats(pirgd_45_2040, stat = 'max', na.rm = T)),
                     data.frame(rcp = 'RCP 4.5', year = 2070:2099,
                                min = cellStats(pirgd_45_2070, stat = 'min', na.rm = T),
                                q1 = quantile(pirgd_45_2070, probs = 0.25, na.rm = T) %>% round,
                                med = cellStats(pirgd_45_2070, stat = 'median', na.rm = T) %>% round,
                                mean = cellStats(pirgd_45_2070, stat = 'mean', na.rm = T) %>% round,
                                q3 = quantile(pirgd_45_2070, probs = 0.75, na.rm = T) %>% round,
                                max = cellStats(pirgd_45_2070, stat = 'max', na.rm = T)))

#clean up
rm(pirgd_45_1999, pirgd_45_2020,
   pirgd_45_2040, pirgd_45_2070,
   pirgd_85_1999, pirgd_85_2020,
   pirgd_85_2040, pirgd_85_2070)

#change col names
colnames(dat_ann_sum) <- c('rcp', 'year', 'min', 'q1',
                           'med', 'mean', 'q3', 'max')

#set plot color scheme
ann_col <- c('Baseline' = '#228833', 
             'BIT RCP 4.5' = '#CCBB44', 
             'BIT RCP 8.5' = '#AA3377',
             'RCP 4.5' = '#4477AA',
             'RCP 8.5' = '#EE6677')

#plot pirgd time series
ggplot(data = dat_ann_sum) +
  geom_boxplot(aes(x = factor(year), color = rcp,
    lower = q1,
    upper = q3,
    middle = med,
    ymin = q1,
    ymax = q3), stat = 'identity') +
  theme_bw() + theme_e +
  scale_color_manual(values = ann_col, name = 'RCP') +
  scale_x_discrete(breaks = c(2000, 2025, 2050, 2075, 2099)) +
  xlab('Year') + ylab('PIRGd') + ggtitle('Average PIRGd Middle Scenario (2000-2099)') +
  geom_vline(xintercept = as.factor(2020))

ggsave('output/final_plots/pirgd_ann_middle_scenario_spatial_variability.png',
       width = 12, height = 8, unit = "in")
