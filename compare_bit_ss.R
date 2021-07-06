#this code compares BIT outputs of pirgd from various scenarios to see if all BIT values
#show similar patterns

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

#create file list of ss projection rasters
files <- list.files('wy_projections/ss', 'ss_ann',
                    full.names = T)
names <- list.files('wy_projections/ss', 'ss_ann')

#loop through projection rasters
for(i in 1:length(files)){
  
  #create dataframe output variables
  str <- str_split(names[i], "_")
  model_h <- str[[1]][3]
  rcp_h <- str[[1]][4]
  period_h <- str[[1]][5] %>% str_replace(., '.tif', '') %>% as.numeric
  
  #load projection brick
  proj_ann <- brick(files[i])
  
  #load projection raster
  proj <- proj_ann %>% calc(., function(x) mean(x, na.rm = T)) %>% getValues %>% round(2)
  
  #only move forward if have values
  if(sum(is.na(proj)) != 22967){
    
    #fill summary dataframe for 3 lc types
    #shrub
    dat <- rbind(dat,
                 data.frame(model = model_h, rcp = rcp_h, period = period_h, lc = 'shrub',
                            min = min(proj[lc_v == 1], na.rm = T),
                            q1 = quantile(proj[lc_v == 1], probs = 0.25, na.rm = T) %>% as.numeric,
                            med = median(proj[lc_v == 1], na.rm = T),
                            mean = mean(proj[lc_v == 1], na.rm = T) %>% round(2),
                            q3 = quantile(proj[lc_v == 1], probs = 0.75, na.rm = T) %>% as.numeric,
                            max = max(proj[lc_v == 1], na.rm = T)))
    
    #herb
    dat <- rbind(dat,
                 data.frame(model = model_h, rcp = rcp_h, period = period_h, lc = 'herb',
                            min = min(proj[lc_v == 2], na.rm = T),
                            q1 = quantile(proj[lc_v == 2], probs = 0.25, na.rm = T) %>% as.numeric,
                            med = median(proj[lc_v == 2], na.rm = T),
                            mean = mean(proj[lc_v == 2], na.rm = T) %>% round(2),
                            q3 = quantile(proj[lc_v == 2], probs = 0.75, na.rm = T) %>% as.numeric,
                            max = max(proj[lc_v == 2], na.rm = T)))
    
    #evergreen
    dat <- rbind(dat,
                 data.frame(model = model_h, rcp = rcp_h, period = period_h, lc = 'evergreen',
                            min = min(proj[lc_v == 3], na.rm = T),
                            q1 = quantile(proj[lc_v == 3], probs = 0.25, na.rm = T) %>% as.numeric,
                            med = median(proj[lc_v == 3], na.rm = T),
                            mean = mean(proj[lc_v == 3], na.rm = T) %>% round(2),
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
                          cellStats(., stat = 'mean', na.rm = T) %>% round(2))
    shrub$model <- model_h
    shrub$rcp <- rcp_h
    shrub$year <- period_h + 1:nlayers(proj_ann) - 1
    shrub$period <- period_h
    shrub$lc <- 'shrub'
    
    #herb
    herb <- data.frame(value = proj_ann %>% mask(., mask = lc, inverse = T, maskvalue = 2) %>% 
                         cellStats(., stat = 'mean', na.rm = T) %>% round(2))
    herb$model <- model_h
    herb$rcp <- rcp_h
    herb$year <- period_h + 1:nlayers(proj_ann) - 1
    herb$period <- period_h
    herb$lc <- 'herb'
    
    #evergreen
    ever <- data.frame(value = proj_ann %>% mask(., mask = lc, inverse = T, maskvalue = 3) %>% 
                         cellStats(., stat = 'mean', na.rm = T) %>% round(2))
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
ss <- raster('dlc/mean_springScale_wy_laea_2001_2018.tif')
mask <- raster('dlc/maxIRGdate_mask_disturbances.tif')

#mask disturbances
ss <- mask(ss, mask)
rm(mask)

#get pirgd values
ss <- getValues(ss)

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
                        min = min(ss[lc == 1], na.rm = T),
                        q1 = quantile(ss[lc == 1], probs = 0.25, na.rm = T) %>% as.numeric,
                        med = median(ss[lc == 1], na.rm = T),
                        mean = mean(ss[lc == 1], na.rm = T) %>% round(2),
                        q3 = quantile(ss[lc == 1], probs = 0.75, na.rm = T) %>% as.numeric,
                        max = max(ss[lc == 1], na.rm = T)))

#herb
dat <- rbind(dat,
             data.frame(model = 'dlc', rcp = 'na', period = 2001, lc = 'herb',
                        min = min(ss[lc == 2], na.rm = T),
                        q1 = quantile(ss[lc == 2], probs = 0.25, na.rm = T) %>% as.numeric,
                        med = median(ss[lc == 2], na.rm = T),
                        mean = mean(ss[lc == 2], na.rm = T) %>% round(2),
                        q3 = quantile(ss[lc == 2], probs = 0.75, na.rm = T) %>% as.numeric,
                        max = max(ss[lc == 2], na.rm = T)))

#evergreen
dat <- rbind(dat,
             data.frame(model = 'dlc', rcp = 'na', period = 2001, lc = 'evergreen',
                        min = min(ss[lc == 3], na.rm = T),
                        q1 = quantile(ss[lc == 3], probs = 0.25, na.rm = T) %>% as.numeric,
                        med = median(ss[lc == 3], na.rm = T),
                        mean = mean(ss[lc == 3], na.rm = T) %>% round(2),
                        q3 = quantile(ss[lc == 3], probs = 0.75, na.rm = T) %>% as.numeric,
                        max = max(ss[lc == 3], na.rm = T)))

#create raw dataframes for 3 landcover types
#name rcp 45 and rcp 85 to simplify figures
#shrub
shrub <- data.frame(value = ss[lc == 1]) %>% na.omit
shrub$model <- 'DLC'
shrub$rcp <- 'rcp45'
shrub$period <- 2001
shrub$lc <- 'shrub'

#herb
herb <- data.frame(value = ss[lc == 2]) %>% na.omit
herb$model <- 'DLC'
herb$rcp <- 'rcp45'
herb$period <- 2001
herb$lc <- 'herb'

#evergreen
ever <- data.frame(value = ss[lc == 3]) %>% na.omit
ever$model <- 'DLC'
ever$rcp <- 'rcp45'
ever$period <- 2001
ever$lc <- 'evergreen'

#rbind raw dataframes
dat_raw <- rbind(dat_raw, shrub, herb, ever)

#rcp85
#shrub
shrub <- data.frame(value = ss[lc == 1]) %>% na.omit
shrub$model <- 'DLC'
shrub$rcp <- 'rcp85'
shrub$period <- 2001
shrub$lc <- 'shrub'

#herb
herb <- data.frame(value = ss[lc == 2]) %>% na.omit
herb$model <- 'DLC'
herb$rcp <- 'rcp85'
herb$period <- 2001
herb$lc <- 'herb'

#evergreen
ever <- data.frame(value = ss[lc == 3]) %>% na.omit
ever$model <- 'DLC'
ever$rcp <- 'rcp85'
ever$period <- 2001
ever$lc <- 'evergreen'

#rbind raw dataframes
dat_raw <- rbind(dat_raw, shrub, herb, ever)

#clean up
rm(ss, lc, shrub, herb, ever)

#calculate contemporary values for annual spatial average
ss_ann <- brick('dlc/springScale_wy_laea_2000_2019_masked.tif')
lc <- raster('landcover/five_class_landcover_wy_laea_2016.tif')
lc[lc == 3] <- NA
lc[lc == 5] <- NA
lc[lc == 4] <- 3

#create ann dataframes for 3 landcover types
#name rcp 45 and rcp 85 to simplify figures
#shrub
shrub <- data.frame(value = ss_ann %>% mask(., mask = lc, inverse = T, maskvalue = 1) %>% 
                      cellStats(., stat = 'mean', na.rm = T) %>% round)
shrub$model <- 'DLC'
shrub$rcp <- 'rcp45'
shrub$year <- 1999 + 1:nlayers(ss_ann)
shrub$period <- 2001
shrub$lc <- 'shrub'

#herb
herb <- data.frame(value = ss_ann %>% mask(., mask = lc, inverse = T, maskvalue = 2) %>% 
                     cellStats(., stat = 'mean', na.rm = T) %>% round)
herb$model <- 'DLC'
herb$rcp <- 'rcp45'
herb$year <- 1999 + 1:nlayers(ss_ann)
herb$period <- 2001
herb$lc <- 'herb'

#evergreen
ever <- data.frame(value = ss_ann %>% mask(., mask = lc, inverse = T, maskvalue = 3) %>% 
                     cellStats(., stat = 'mean', na.rm = T) %>% round)
ever$model <- 'DLC'
ever$rcp <- 'rcp45'
ever$year <- 1999 + 1:nlayers(ss_ann)
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
rm(ss_ann, lc, shrub, herb, ever)

####################
###CREATE FIGURES###
####################

#clean data for erroneous values
dat_raw2 <- dat_raw
dat_raw2$value[dat_raw2$value > 50] <- 50

#only keep BIT data
dat_raw2 <- dat_raw2[dat_raw2$period == 1999,]

#set models to factor levels
dat_raw2$model <- as.factor(dat_raw2$model)

#set color scale for all 6 models
#need to choose a color palette colorblind friendly
#https://personal.sron.nl/~pault/#sec:qualitative
myColors <- c('#4477AA', '#228833', '#CCBB44', '#EE6677', 
              '#AA3377')
names(myColors) <- levels(dat_raw2$model)
colScale <- scale_color_manual(name = "Model", values = myColors,
                               breaks = c('CNRM-CM5', 'DLC', 'HadGEM2-ES365', 'inmcm4',
                                          'IPSL-CM5A-MR', 'MIROC-ESM-CHEM'),
                               labels = c('CNRM-CM5' = 'CNRM-CM5',
                                          'DLC' = 'Baseline',
                                          'HadGEM2-ES365' = 'HadGEM2-ES365',
                                          'inmcm4' = 'inmcm4',
                                          'IPSL-CM5A-MR' = 'IPSL-CM5A-MR',
                                          'MIROC-ESM-CHEM' = 'MIROC-ESM-CHEM'))
fillScale <- scale_fill_manual(name = "Model", values = myColors)

#set universal theme
theme_e <- theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                 axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

#plot rcp 45
ggplot(data = dat_raw2[dat_raw2$rcp == 'rcp45',], aes(x = lc, color = model, y = value)) +
  geom_boxplot() + 
  colScale + xlab('Landcover') + ylab('Spring Scale (Days)') +
  ggtitle('BIT Comparison RCP 4.5') + theme_bw() + theme_e +
  scale_y_continuous(breaks = c(10, 15, 20, 25), labels = c('      10', '      15', '      20', '      25'))

ggsave('output/final_plots/bit_comparisons_rcp45_ss.png',
       width = 9, height = 5, unit = "in")

#plot rcp 85
ggplot(data = dat_raw2[dat_raw2$rcp == 'rcp85',], aes(x = lc, color = model, y = value)) +
  geom_boxplot() + 
  colScale + xlab('Landcover') + ylab('Spring Scale (Days)') +
  ggtitle('BIT Comparison RCP 8.5') + theme_bw() + theme_e +
  scale_y_continuous(breaks = c(10, 15, 20, 25), labels = c('      10', '      15', '      20', '      25'))

ggsave('output/final_plots/bit_comparisons_rcp85_ss.png',
       width = 9, height = 5, unit = "in")


