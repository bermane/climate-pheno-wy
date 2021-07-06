#this code compares future outputs of pirgd and ss to their contemporary
#counterparts to see if the patterns are similar

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

###############################
###LOAD CONTEMPORARY DATASET###
###############################

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

#set zero and decimal gdd, rain and vpd values to 1 from variables we will log transform
dat2$gdd_jan_apr[dat2$gdd_jan_apr < 1] <- 1
dat2$rain_mar_may[dat2$rain_mar_may < 1] <- 1
dat2$vpd_tavg_mean_jan_apr[dat2$vpd_tavg_mean_jan_apr < 1] <- 1

#any other missing points?
#may need to remove the variables that include snowmelt and pirgd...
#16% of samples have PIRGd before snowmelt...
sum(is.na(dat2))

#id col and lc as factors
dat2$id <- as.factor(dat2$id)
dat2$lc <- as.factor(dat2$lc)
levels(dat2$lc) <- c('shrub', 'herb', 'evergreen')


##############################
###LOAD FUTURE PIRGd AND SS###
##############################

#load pirgd and ss files. only middle scenario
#make sure the scenarios match the index locations
pirgd_files <- list.files('wy_projections/pirgd_log_gdd_log_rain', glob2rx('pirgd_mean*HadGEM2*.tif'), full.names = T)
ss_files <- list.files('wy_projections/ss', glob2rx('ss_mean*HadGEM2*.tif'), full.names = T)

#load different scenarios into lists so can check different ones
rcp45_99 <- data.frame(pirgd = getValues(raster(pirgd_files[1])),
                       ss = getValues(raster(ss_files[1])))
rcp45_40 <- data.frame(pirgd = getValues(raster(pirgd_files[3])),
                       ss = getValues(raster(ss_files[3])))
rcp45_70 <- data.frame(pirgd = getValues(raster(pirgd_files[4])),
                       ss = getValues(raster(ss_files[4])))
rcp85_99 <- data.frame(pirgd = getValues(raster(pirgd_files[5])),
                       ss = getValues(raster(ss_files[5])))
rcp85_40 <- data.frame(pirgd = getValues(raster(pirgd_files[7])),
                       ss = getValues(raster(ss_files[7])))
rcp85_70 <- data.frame(pirgd = getValues(raster(pirgd_files[8])),
                       ss = getValues(raster(ss_files[8])))
######################
###PLOT PIRGd VS SS###
######################

#set universal theme
theme_e <- theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15), legend.text=element_text(size = 15),
                 axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                 legend.title = element_text(size = 15))

#plot contemporary data
ggplot(dat2, aes(x = pirgd, y = ss)) +
  geom_point() + theme_bw() + theme_e +
  xlab('PIRGd') + ylab('Spring Scale') + ggtitle('Spring Scale vs. PIRGd \nContemporary (2001-2019)') +
  scale_x_continuous(limits = c(1, 300)) + scale_y_continuous(limits = c(0, 35))

ggsave('output/final_plots/log_model/pirgd_ss/pirgd_ss_dlc.png')

 #plot climate scenario data
ggplot(rcp45_99, aes(x = pirgd, y = ss)) +
  geom_point() + theme_bw() + theme_e +
  xlab('PIRGd') + ylab('Spring Scale') + ggtitle('Spring Scale vs. PIRGd \nMiddle Scenario RCP 4.5 (2000-2019)') +
  scale_x_continuous(limits = c(1, 300)) + scale_y_continuous(limits = c(0, 35))

ggsave('output/final_plots/log_model/pirgd_ss/pirgd_ss_hadgem_45_2000_2019.png')

ggplot(rcp85_99, aes(x = pirgd, y = ss)) +
  geom_point() + theme_bw() + theme_e +
  xlab('PIRGd') + ylab('Spring Scale') + ggtitle('Spring Scale vs. PIRGd \nMiddle Scenario RCP 8.5 (2000-2019)') +
  scale_x_continuous(limits = c(1, 300)) + scale_y_continuous(limits = c(0, 35))

ggsave('output/final_plots/log_model/pirgd_ss/pirgd_ss_hadgem_85_2000_2019.png')

ggplot(rcp45_40, aes(x = pirgd, y = ss)) +
  geom_point() + theme_bw() + theme_e +
  xlab('PIRGd') + ylab('Spring Scale') + ggtitle('Spring Scale vs. PIRGd \nMiddle Scenario RCP 4.5 (2040-2069)') +
  scale_x_continuous(limits = c(1, 300)) + scale_y_continuous(limits = c(0, 35))

ggsave('output/final_plots/log_model/pirgd_ss/pirgd_ss_hadgem_45_2040_2069.png')

ggplot(rcp85_40, aes(x = pirgd, y = ss)) +
  geom_point() + theme_bw() + theme_e +
  xlab('PIRGd') + ylab('Spring Scale') + ggtitle('Spring Scale vs. PIRGd \nMiddle Scenario RCP 8.5 (2040-2069)') +
  scale_x_continuous(limits = c(1, 300)) + scale_y_continuous(limits = c(0, 35))

ggsave('output/final_plots/log_model/pirgd_ss/pirgd_ss_hadgem_85_2040_2069.png')

ggplot(rcp45_70, aes(x = pirgd, y = ss)) +
  geom_point() + theme_bw() + theme_e +
  xlab('PIRGd') + ylab('Spring Scale') + ggtitle('Spring Scale vs. PIRGd \nMiddle Scenario RCP 4.5 (2070-2099)') +
  scale_x_continuous(limits = c(1, 300)) + scale_y_continuous(limits = c(0, 35))

ggsave('output/final_plots/log_model/pirgd_ss/pirgd_ss_hadgem_45_2070_2099.png')

ggplot(rcp85_70, aes(x = pirgd, y = ss)) +
  geom_point() + theme_bw() + theme_e +
  xlab('PIRGd') + ylab('Spring Scale') + ggtitle('Spring Scale vs. PIRGd \nMiddle Scenario RCP 8.5 (2070-2099)') +
  scale_x_continuous(limits = c(1, 300)) + scale_y_continuous(limits = c(0, 35))

ggsave('output/final_plots/log_model/pirgd_ss/pirgd_ss_hadgem_85_2070_2099.png')

