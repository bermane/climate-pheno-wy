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

#create file list of pirgd projection rasters
files <- list.files('wy_projections/pirgd', 'pirgd_ann',
                         full.names = T)
names <- list.files('wy_projections/pirgd', 'pirgd_ann')

#loop through projection rasters
for(i in 1:length(files)){
  
  #load projection raster
  proj <- brick(files[i]) %>% calc(., function(x) mean(x, na.rm = T)) %>% getValues %>% round
  
  #create dataframe output variables
  str <- str_split(names[i], "_")
  model_h <- str[[1]][3]
  rcp_h <- str[[1]][4]
  period_h <- str[[1]][5] %>% str_replace(., '.tif', '') %>% as.numeric
  
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
  rm(proj, str, model_h, rcp_h, period_h, shrub, herb, ever)
  
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

####################
###CREATE FIGURES###
####################

#clean data for erroneous values
dat_raw2 <- dat_raw[dat_raw$value < 366 & dat_raw$value > 0,]
dat_raw2$model <- as.factor(dat_raw2$model)

#set color scale for all 6 models
#need to choose a color palette colorblind friendly
#https://personal.sron.nl/~pault/#sec:qualitative
myColors <- c('#4477AA', '#BBBBBB', '#228833', '#CCBB44', '#EE6677', 
              '#AA3377')
names(myColors) <- levels(dat_raw2$model)
colScale <- scale_color_manual(name = "Model", values = myColors, labels = c('DLC' = 'Baseline (2001-2018)'))
fillScale <- scale_fill_manual(name = "Model", values = myColors)

#set universal theme
theme_e <- theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                 axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

#PLOT OVERALL IN STATE

#organize rcp45 and rcp85 data separately based on median value
#2040
dat_rcp45_2040 <- dat_raw2[dat_raw2$rcp == 'rcp45' & (dat_raw2$period %in% c(2040, 2001)),]
dat_rcp45_2040$model <- fct_reorder(dat_rcp45_2040$model, dat_rcp45_2040$value, .fun = median)
dat_rcp85_2040 <- dat_raw2[dat_raw2$rcp == 'rcp85' & (dat_raw2$period %in% c(2040, 2001)),]
dat_rcp85_2040$model <- fct_reorder(dat_rcp85_2040$model, dat_rcp85_2040$value, .fun = median)

#2070
dat_rcp45_2070 <- dat_raw2[dat_raw2$rcp == 'rcp45' & (dat_raw2$period %in% c(2070, 2001)),]
dat_rcp45_2070$model <- fct_reorder(dat_rcp45_2070$model, dat_rcp45_2070$value, .fun = median)
dat_rcp85_2070 <- dat_raw2[dat_raw2$rcp == 'rcp85' & (dat_raw2$period %in% c(2070, 2001)),]
dat_rcp85_2070$model <- fct_reorder(dat_rcp85_2070$model, dat_rcp85_2070$value, .fun = median)

#PLOT BY LC TYPE
#organize data by rcp and lc type
#2040
dat_rcp45_2040_shrub <- dat_rcp45_2040[dat_rcp45_2040$lc == 'shrub',]
dat_rcp45_2040_shrub$model <- fct_reorder(dat_rcp45_2040_shrub$model, dat_rcp45_2040_shrub$value, .fun = median)
dat_rcp85_2040_shrub <- dat_rcp85_2040[dat_rcp85_2040$lc == 'shrub',]
dat_rcp85_2040_shrub$model <- fct_reorder(dat_rcp85_2040_shrub$model, dat_rcp85_2040_shrub$value, .fun = median)

dat_rcp45_2040_herb <- dat_rcp45_2040[dat_rcp45_2040$lc == 'herb',]
dat_rcp45_2040_herb$model <- fct_reorder(dat_rcp45_2040_herb$model, dat_rcp45_2040_herb$value, .fun = median)
dat_rcp85_2040_herb <- dat_rcp85_2040[dat_rcp85_2040$lc == 'herb',]
dat_rcp85_2040_herb$model <- fct_reorder(dat_rcp85_2040_herb$model, dat_rcp85_2040_herb$value, .fun = median)

dat_rcp45_2040_ever <- dat_rcp45_2040[dat_rcp45_2040$lc == 'evergreen',]
dat_rcp45_2040_ever$model <- fct_reorder(dat_rcp45_2040_ever$model, dat_rcp45_2040_ever$value, .fun = median)
dat_rcp85_2040_ever <- dat_rcp85_2040[dat_rcp85_2040$lc == 'evergreen',]
dat_rcp85_2040_ever$model <- fct_reorder(dat_rcp85_2040_ever$model, dat_rcp85_2040_ever$value, .fun = median)

#2070
dat_rcp45_2070_shrub <- dat_rcp45_2070[dat_rcp45_2070$lc == 'shrub',]
dat_rcp45_2070_shrub$model <- fct_reorder(dat_rcp45_2070_shrub$model, dat_rcp45_2070_shrub$value, .fun = median)
dat_rcp85_2070_shrub <- dat_rcp85_2070[dat_rcp85_2070$lc == 'shrub',]
dat_rcp85_2070_shrub$model <- fct_reorder(dat_rcp85_2070_shrub$model, dat_rcp85_2070_shrub$value, .fun = median)

dat_rcp45_2070_herb <- dat_rcp45_2070[dat_rcp45_2070$lc == 'herb',]
dat_rcp45_2070_herb$model <- fct_reorder(dat_rcp45_2070_herb$model, dat_rcp45_2070_herb$value, .fun = median)
dat_rcp85_2070_herb <- dat_rcp85_2070[dat_rcp85_2070$lc == 'herb',]
dat_rcp85_2070_herb$model <- fct_reorder(dat_rcp85_2070_herb$model, dat_rcp85_2070_herb$value, .fun = median)

dat_rcp45_2070_ever <- dat_rcp45_2070[dat_rcp45_2070$lc == 'evergreen',]
dat_rcp45_2070_ever$model <- fct_reorder(dat_rcp45_2070_ever$model, dat_rcp45_2070_ever$value, .fun = median)
dat_rcp85_2070_ever <- dat_rcp85_2070[dat_rcp85_2070$lc == 'evergreen',]
dat_rcp85_2070_ever$model <- fct_reorder(dat_rcp85_2070_ever$model, dat_rcp85_2070_ever$value, .fun = median)


#create rcp and landcover variable so can plot all together
#2040
dat_rcp45_2040$group <- 'DAOverall_RCP45'
dat_rcp85_2040$group <- 'DBOverall_RCP85'
dat_rcp45_2040_shrub$group <- 'BAShrub_RCP45'
dat_rcp85_2040_shrub$group <- 'BBShrub_RCP85'
dat_rcp45_2040_herb$group <- 'AAHerb_RCP45'
dat_rcp85_2040_herb$group <- 'ABHerb_RCP85'
dat_rcp45_2040_ever$group <- 'CAEvergreen_RCP45'
dat_rcp85_2040_ever$group <- 'CBEvergreen_RCP85'

#2070
dat_rcp45_2070$group <- 'DAOverall_RCP45'
dat_rcp85_2070$group <- 'DBOverall_RCP85'
dat_rcp45_2070_shrub$group <- 'BAShrub_RCP45'
dat_rcp85_2070_shrub$group <- 'BBShrub_RCP85'
dat_rcp45_2070_herb$group <- 'AAHerb_RCP45'
dat_rcp85_2070_herb$group <- 'ABHerb_RCP85'
dat_rcp45_2070_ever$group <- 'CAEvergreen_RCP45'
dat_rcp85_2070_ever$group <- 'CBEvergreen_RCP85'

#create new df with median values and quartiles of DLC to plot as error bar
dat_baseline <- data.frame(group = c('DAOverall_RCP45', 'DBOverall_RCP85',
                                     'BAShrub_RCP45', 'BBShrub_RCP85',
                                     'AAHerb_RCP45', 'ABHerb_RCP85',
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
dat_baseline2 <- dat_baseline[dat_baseline$group %in% c('AAHerb_RCP45', 'DBOverall_RCP85'),]
dat_baseline <- dat_baseline[!(dat_baseline$group %in% c('AAHerb_RCP45', 'DBOverall_RCP85')),]

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
  scale_x_discrete(labels = c('AAHerb_RCP45' = 'Herb \nRCP 4.5', 'ABHerb_RCP85' = 'Herb \nRCP 8.5',
                              'BAShrub_RCP45' = 'Shrub \nRCP 4.5', 'BBShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5) +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = med_line, ymin = med_line), 
                colour="#000000", linetype="dashed") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = med_line, ymin = med_line), 
                colour="#000000", linetype="dashed") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q1, ymin = q1), 
                colour="#BBBBBB", linetype="dashed") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q1, ymin = q1), 
                colour="#BBBBBB", linetype="dashed") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q3, ymin = q3), 
                colour="#BBBBBB", linetype="dashed") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q3, ymin = q3), 
                colour="#BBBBBB", linetype="dashed")

ggsave('output/final_plots/pirgd_proj_combine_2040.png',
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
  scale_x_discrete(labels = c('AAHerb_RCP45' = 'Herb \nRCP 4.5', 'ABHerb_RCP85' = 'Herb \nRCP 8.5',
                              'BAShrub_RCP45' = 'Shrub \nRCP 4.5', 'BBShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5) +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = med_line, ymin = med_line), 
                colour="#000000", linetype="dashed") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = med_line, ymin = med_line), 
                colour="#000000", linetype="dashed") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q1, ymin = q1), 
                colour="#BBBBBB", linetype="dashed") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q1, ymin = q1), 
                colour="#BBBBBB", linetype="dashed") +
  geom_errorbar(data = dat_baseline, width = 1, aes(x=group, ymax = q3, ymin = q3), 
                colour="#BBBBBB", linetype="dashed") +
  geom_errorbar(data = dat_baseline2, width = 1.2, aes(x=group, ymax = q3, ymin = q3), 
                colour="#BBBBBB", linetype="dashed")

ggsave('output/final_plots/pirgd_proj_combine_2070.png',
       width = 12, height = 8, unit = "in")

#create boxplot with baseline boxes

ggplot() +
  geom_boxplot(data = dat_rcp45_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_herb, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_shrub, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85_ever, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp45, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  geom_boxplot(data = dat_rcp85, 
               aes(x = group, y = value, color = model), outlier.size=0.5) +
  colScale + xlab('Landcover and Emissions Scenario') + ylab('PIRGd (DOY)') +
  ggtitle('PIRGd Projections for Mid-Century in WY (2040-2069)') + theme_bw() + theme_e +
  scale_y_continuous(limits = c(0, 365)) + 
  scale_x_discrete(labels = c('AAHerb_RCP45' = 'Herb \nRCP 4.5', 'ABHerb_RCP85' = 'Herb \nRCP 8.5',
                              'BAShrub_RCP45' = 'Shrub \nRCP 4.5', 'BBShrub_RCP85' = 'Shrub \nRCP 8.5',
                              'CAEvergreen_RCP45' = 'Evergreen \nRCP 4.5', 'CBEvergreen_RCP85' = 'Evergreen \nRCP 8.5',
                              'DAOverall_RCP45' = 'Overall \nRCP 4.5', 'DBOverall_RCP85' = 'Overall \nRCP 8.5')) + 
  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = 'black', size = 0.5)

ggsave('output/final_plots/pirgd_proj_combine_2040_basebox.png',
       width = 12, height = 8, unit = "in")

#RE RUN FOR 2070 DATA!

