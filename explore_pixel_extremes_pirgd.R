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

#load middle scenario BIT values
gem_1999_rcp45 <- read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp45_1999.csv') 
gem_1999_rcp85 <- read.csv('wy_projections/dat/dat_HadGEM2-ES365_rcp85_1999.csv')

#load vpd values
gem_1999_rcp45$vpd <- read.csv('wy_projections/dat/dat_vpd_HadGEM2-ES365_rcp45_1999.csv')$vpd 
gem_1999_rcp85$vpd <- read.csv('wy_projections/dat/dat_vpd_HadGEM2-ES365_rcp85_1999.csv')$vpd 

#only keep years 2005 and 2015 to keep as examples
gem_2005_rcp45 <- gem_1999_rcp45[gem_1999_rcp45$year == 2005,]
gem_2005_rcp85 <- gem_1999_rcp85[gem_1999_rcp85$year == 2005,]
gem_2015_rcp45 <- gem_1999_rcp45[gem_1999_rcp45$year == 2015,]
gem_2015_rcp85 <- gem_1999_rcp85[gem_1999_rcp85$year == 2015,]

#clean up
rm(gem_1999_rcp45, gem_1999_rcp85)

#add landcover to dataframes
gem_2005_rcp45$lc <- lc_v
gem_2005_rcp85$lc <- lc_v
gem_2015_rcp45$lc <- lc_v
gem_2015_rcp85$lc <- lc_v

#load pirgd values for same years/models. these will be base of raster stacks
ras_gem_2005_rcp45 <- raster('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_1999.tif', band = 7)
ras_gem_2005_rcp85 <- raster('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_1999.tif', band = 7)
ras_gem_2015_rcp45 <- raster('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp45_1999.tif', band = 17)
ras_gem_2015_rcp85 <- raster('wy_projections/pirgd_log_gdd_log_rain/pirgd_ann_HadGEM2-ES365_rcp85_1999.tif', band = 17)

#set ss values to min 1 and max 365
ras_gem_2005_rcp45[ras_gem_2005_rcp45 > 365] <- 365
ras_gem_2005_rcp45[ras_gem_2005_rcp45 < 1] <- 1
ras_gem_2005_rcp85[ras_gem_2005_rcp85 > 365] <- 365
ras_gem_2005_rcp85[ras_gem_2005_rcp85 < 1] <- 1
ras_gem_2015_rcp45[ras_gem_2015_rcp45 > 365] <- 365
ras_gem_2015_rcp45[ras_gem_2015_rcp45 < 1] <- 1
ras_gem_2015_rcp85[ras_gem_2015_rcp85 > 365] <- 365
ras_gem_2015_rcp85[ras_gem_2015_rcp85 < 1] <- 1

#ensure rasters have same number of cells as data frames
ncell(ras_gem_2005_rcp45)
nrow(gem_2005_rcp45)

#create raster stack of covariates on top of pirgd
covar <- c('rain', 'snow', 'gdd', 'pdsi', 'lc')

#loop through co-variates
for(vari in covar){
  
  #load raster template
  ras <- clim_ras
  
  #gem_2005_rcp45
  ras[] <- gem_2005_rcp45[,vari] %>% round(2)
  ras_gem_2005_rcp45 <- addLayer(ras_gem_2005_rcp45, ras)
  
  #gem_2005_rcp85
  ras[] <- gem_2005_rcp85[,vari] %>% round(2)
  ras_gem_2005_rcp85 <- addLayer(ras_gem_2005_rcp85, ras)
  
  #gem_2015_rcp45
  ras[] <- gem_2015_rcp45[,vari] %>% round(2)
  ras_gem_2015_rcp45 <- addLayer(ras_gem_2015_rcp45, ras)
  
  #gem_2015_rcp85
  ras[] <- gem_2015_rcp85[,vari] %>% round(2)
  ras_gem_2015_rcp85 <- addLayer(ras_gem_2015_rcp85, ras)
  
}

#clean up
rm(ras)

#change layer names
names(ras_gem_2005_rcp45) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')
names(ras_gem_2005_rcp85) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')
names(ras_gem_2015_rcp45) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')
names(ras_gem_2015_rcp85) <- c('pirgd', 'rain', 'snow', 'gdd', 'pdsi', 'lc')

############################
###LOAD CONTEMPORARY DATA###
############################
#load workspace
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_10_18.RData')

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

#only keep columns we want
rm(dat)
dat2 <- dat2[, colnames(dat2) %in% c('year', 'pirgd', 'lc', 'pdsi_mar_apr_min', 'rain_mar_may', 'snow_oct_apr', 'gdd_jan_apr')]
colnames(dat2) <- c('year', 'pirgd_dlc', 'lc_dlc', 'pdsi_dlc', 'snow_dlc', 'gdd_dlc', 'rain_dlc')
dat2 <- dat2[, c('year', 'pirgd_dlc', 'rain_dlc', 'snow_dlc', 'gdd_dlc', 'pdsi_dlc', 'lc_dlc')]

#subset 2005 and 2015
dat_2005 <- dat2[dat2$year == 2005,]
dat_2015 <- dat2[dat2$year == 2015,]

#############################
###SAMPLE BIT CLIMATE DATA###
#############################

#load sampling grid
grid <- readOGR('reference/sampling_points_three_equal_classes.shp')

#change sampling grid to climate proj
grid <- spTransform(grid, crs(ras_gem_2005_rcp45))

#sample future points on grid
dat_gem_2005_rcp45 <- raster::extract(ras_gem_2005_rcp45, grid, df = T)
dat_gem_2005_rcp85 <- raster::extract(ras_gem_2005_rcp85, grid, df = T)
dat_gem_2015_rcp45 <- raster::extract(ras_gem_2015_rcp45, grid, df = T)
dat_gem_2015_rcp85 <- raster::extract(ras_gem_2015_rcp85, grid, df = T)

#cbind comtemporary data
dat_gem_2005_rcp45 <- cbind(dat_gem_2005_rcp45, dat_2005[, c('pirgd_dlc', 'rain_dlc', 'snow_dlc', 'gdd_dlc', 'pdsi_dlc', 'lc_dlc')])
dat_gem_2005_rcp85 <- cbind(dat_gem_2005_rcp85, dat_2005[, c('pirgd_dlc', 'rain_dlc', 'snow_dlc', 'gdd_dlc', 'pdsi_dlc', 'lc_dlc')])
dat_gem_2015_rcp45 <- cbind(dat_gem_2015_rcp45, dat_2015[, c('pirgd_dlc', 'rain_dlc', 'snow_dlc', 'gdd_dlc', 'pdsi_dlc', 'lc_dlc')])
dat_gem_2015_rcp85 <- cbind(dat_gem_2015_rcp85, dat_2015[, c('pirgd_dlc', 'rain_dlc', 'snow_dlc', 'gdd_dlc', 'pdsi_dlc', 'lc_dlc')])

#sample specific rows
dat_samples <- rbind(slice_sample(dat_gem_2005_rcp45, n = 100),
                     slice_sample(dat_gem_2005_rcp85, n = 100),
                     slice_sample(dat_gem_2015_rcp45, n = 100),
                     slice_sample(dat_gem_2015_rcp85, n = 100))

dat_samples <- cbind(data.frame(scenario = c(rep('2005_rcp45', 100),
                                             rep('2005_rcp85', 100),
                                             rep('2015_rcp45', 100),
                                             rep('2015_rcp85', 100))),
                     dat_samples)

#change lc col to match types
dat_samples$lc <- as.factor(dat_samples$lc)
levels(dat_samples$lc) <- c('shrub', 'herb', 'evergreen')

#round
dat_samples$rain_dlc <- round(dat_samples$rain_dlc, 2)
dat_samples$snow_dlc <- round(dat_samples$snow_dlc, 2)
dat_samples$gdd_dlc <- round(dat_samples$gdd_dlc, 2)

#remove id col
#dat_samples <- dat_samples[, !(colnames(dat_samples) %in% 'ID')]

#write to csv
write.csv(dat_samples, 'output/final/log_gdd_log_rain/bit_pirgd_comparison.csv', row.names = F)

#plot
ggplot(dat_samples) +
  geom_point(aes(x = rain, y = pirgd, col = 'red')) +
  geom_point(aes(x = rain_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Rain Mar-May (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Rain') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_rain.png', width = 12,
       height = 8)

ggplot(dat_samples) +
  geom_point(aes(x = snow, y = pirgd, col = 'red')) +
  geom_point(aes(x = snow_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Snow Oct-Apr (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Snow') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_snow.png', width = 12,
       height = 8)

 ggplot(dat_samples) +
  geom_point(aes(x = gdd, y = pirgd, col = 'red')) +
  geom_point(aes(x = gdd_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('GDD Jan-Apr') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: GDD') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))
 
 ggsave('output/final_plots/pirgd_ss/bit_pirgd_gdd.png', width = 12,
        height = 8)

ggplot(dat_samples) +
  geom_point(aes(x = pdsi, y = pirgd, col = 'red')) +
  geom_point(aes(x = pdsi_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Min PDSI Mar-Apr') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: VPD') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_pdsi.png', width = 12,
       height = 8)

#plot by landcover
#plot
ggplot(dat_samples[dat_samples$lc == 'evergreen',]) +
  geom_point(aes(x = rain, y = pirgd, col = 'red')) +
  geom_point(aes(x = rain_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Rain Mar-May (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Rain in Evergreen') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_evergreen_rain.png', width = 12,
       height = 8)

ggplot(dat_samples[dat_samples$lc == 'shrub',]) +
  geom_point(aes(x = rain, y = pirgd, col = 'red')) +
  geom_point(aes(x = rain_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Rain Mar-May (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Rain in Shrub') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_shrub_rain.png', width = 12,
       height = 8)

ggplot(dat_samples[dat_samples$lc == 'herb',]) +
  geom_point(aes(x = rain, y = pirgd, col = 'red')) +
  geom_point(aes(x = rain_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Rain Mar-May (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Rain in Herb') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_herb_rain.png', width = 12,
       height = 8)

ggplot(dat_samples[dat_samples$lc == 'evergreen',]) +
  geom_point(aes(x = snow, y = pirgd, col = 'red')) +
  geom_point(aes(x = snow_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Snow Oct-Apr (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Snow in Evergreen') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_evergreen_snow.png', width = 12,
       height = 8)

ggplot(dat_samples[dat_samples$lc == 'shrub',]) +
  geom_point(aes(x = snow, y = pirgd, col = 'red')) +
  geom_point(aes(x = snow_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Snow Oct-Apr (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Snow in Shrub') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_shrub_snow.png', width = 12,
       height = 8)


ggplot(dat_samples[dat_samples$lc == 'herb',]) +
  geom_point(aes(x = snow, y = pirgd, col = 'red')) +
  geom_point(aes(x = snow_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('Snow Oct-Apr (mm)') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: Snow in Herb') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_herb_snow.png', width = 12,
       height = 8)


ggplot(dat_samples[dat_samples$lc == 'evergreen',]) +
  geom_point(aes(x = gdd, y = pirgd, col = 'red')) +
  geom_point(aes(x = gdd_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('GDD Jan-Apr') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: GDD in Evergreen') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_evergreen_gdd.png', width = 12,
       height = 8)

ggplot(dat_samples[dat_samples$lc == 'shrub',]) +
  geom_point(aes(x = gdd, y = pirgd, col = 'red')) +
  geom_point(aes(x = gdd_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('GDD Jan-Apr') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: GDD in Shrub') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_shrub_gdd.png', width = 12,
       height = 8)

ggplot(dat_samples[dat_samples$lc == 'herb',]) +
  geom_point(aes(x = gdd, y = pirgd, col = 'red')) +
  geom_point(aes(x = gdd_dlc, y = pirgd_dlc, col = 'blue')) +
  theme_bw() + xlab('GDD Jan-Apr') + ylab('PIRGd (DOY)') +
  ggtitle('BIT vs. DLC Comparisons 2005 & 2015: GDD in Herb') + 
  scale_color_discrete(name = 'Scenario', labels = c('Observed', 'Projected'))

ggsave('output/final_plots/pirgd_ss/bit_pirgd_herb_gdd.png', width = 12,
       height = 8)

#look at differences in co variates
dat_gem_2005_rcp45$name <- '2005_rcp45'
dat_gem_2005_rcp85$name <- '2005_rcp85'
dat_gem_2015_rcp45$name <- '2015_rcp45'
dat_gem_2015_rcp85$name <- '2015_rcp85'

#create melted df of all data
dat_melt <- rbind(dat_gem_2005_rcp45[,c('ss', 'rain', 'vpd', 'pdsi', 'lc', 'name')],
                  dat_gem_2005_rcp85[,c('ss', 'rain', 'vpd', 'pdsi', 'lc', 'name')],
                  dat_gem_2015_rcp45[,c('ss', 'rain', 'vpd', 'pdsi', 'lc', 'name')],
                  dat_gem_2015_rcp85[,c('ss', 'rain', 'vpd', 'pdsi', 'lc', 'name')],
                  data.frame(ss = dat_gem_2005_rcp45$ss_dlc,
                             rain = dat_gem_2005_rcp45$rain_dlc,
                             vpd = dat_gem_2005_rcp45$vpd_dlc,
                             pdsi = dat_gem_2005_rcp45$pdsi_dlc,
                             lc = dat_gem_2005_rcp45$lc_dlc,
                             name = 'baseline_2005'),
                  data.frame(ss = dat_gem_2015_rcp45$ss_dlc,
                             rain = dat_gem_2015_rcp45$rain_dlc,
                             vpd = dat_gem_2015_rcp45$vpd_dlc,
                             pdsi = dat_gem_2015_rcp45$pdsi_dlc,
                             lc = dat_gem_2015_rcp45$lc_dlc,
                             name = 'baseline_2015'))

#plot histograms of values
#ss
ggplot(dat_melt, aes(x = factor(name), y = ss)) +
  geom_boxplot() + theme_bw() +
  xlab('Scenario') + ylab('Spring Scale (Days)') +
  ggtitle('Comparison of SS Covariates: SS')

ggsave('output/final_plots/pirgd_ss/bit_ss_box.png', width = 12,
       height = 8)


#rain
ggplot(dat_melt, aes(x = factor(name), y = rain)) +
  geom_boxplot() + theme_bw() +
  xlab('Scenario') + ylab('Rain Mar-May (mm)') +
  ggtitle('Comparison of SS Covariates: Rain')

ggsave('output/final_plots/pirgd_ss/bit_rain_box.png', width = 12,
       height = 8)

#vpd
ggplot(dat_melt, aes(x = factor(name), y = vpd)) +
  geom_boxplot() + theme_bw() +
  xlab('Scenario') + ylab('Mean VPD Jan-Apr') +
  ggtitle('Comparison of SS Covariates: VPD')

ggsave('output/final_plots/pirgd_ss/bit_vpd_box.png', width = 12,
       height = 8)

#pdsi
ggplot(dat_melt, aes(x = factor(name), y = pdsi)) +
  geom_boxplot() + theme_bw() +
  xlab('Scenario') + ylab('Min PDSI Mar-Apr') +
  ggtitle('Comparison of SS Covariates: VPD')

ggsave('output/final_plots/pirgd_ss/bit_pdsi_box.png', width = 12,
       height = 8)

#LOOK AT DIFFERENCES IN VALUES

# #create output df of values from different pixel locations
# pix <- data.frame(model = character(), location = character(), ss = numeric(), 
#                   rain = numeric(), vpd = numeric(), pdsi = numeric(),
#                   landcover = numeric())
# 
# #click gem_2005_rcp45
# plot(ras_gem_2005_rcp45[[1]])
# vals <- click(ras_gem_2005_rcp45)
# 
# #rbind to df
# pix <- rbind(pix,
#              data.frame(model = 'HadGEM 2005 RCP 4.5',
#                         location = c('SE', 'SW', 
#                                      'NE', 'NW'),
#                         ss = vals[,'ss'],
#                         rain = vals[,'rain'],
#                         vpd = vals[,'vpd'],
#                         pdsi = vals[,'pdsi'],
#                         landcover = vals[,'lc']))
# 
# #click gem_2005_rcp85
# plot(ras_gem_2005_rcp85[[1]])
# vals <- click(ras_gem_2005_rcp85)
# 
# #rbind to df
# pix <- rbind(pix,
#              data.frame(model = 'HadGEM 2005 RCP 8.5',
#                         location = c('SE', 'SW', 
#                                      'NE', 'NW'),
#                         ss = vals[,'ss'],
#                         rain = vals[,'rain'],
#                         vpd = vals[,'vpd'],
#                         pdsi = vals[,'pdsi'],
#                         landcover = vals[,'lc']))
# 
# #click gem_2015_rcp45
# plot(ras_gem_2015_rcp45[[1]])
# vals <- click(ras_gem_2015_rcp45)
# 
# #rbind to df
# pix <- rbind(pix,
#              data.frame(model = 'HadGEM 2015 RCP 4.5',
#                         location = c('SE', 'SW', 
#                                      'NE', 'NW'),
#                         ss = vals[,'ss'],
#                         rain = vals[,'rain'],
#                         vpd = vals[,'vpd'],
#                         pdsi = vals[,'pdsi'],
#                         landcover = vals[,'lc']))
# 
# #click inmcm_2071
# plot(ras_gem_2015_rcp85[[1]])
# vals <- click(ras_gem_2015_rcp85)
# 
# #rbind to df
# pix <- rbind(pix,
#              data.frame(model = 'HadGEM 2015 RCP 8.5',
#                         location = c('SE', 'SW', 
#                                      'NE', 'NW'),
#                         ss = vals[,'ss'],
#                         rain = vals[,'rain'],
#                         vpd = vals[,'vpd'],
#                         pdsi = vals[,'pdsi'],
#                         landcover = vals[,'lc']))
# 
# # export pix as csv
# # write.csv(pix, file = 'output/pixel_extremes_2020_12_29.csv',
# #          row.names = F, col.names = F)
# 
# ##############################
# ###CALCULATE CUMULATIVE GDD###
# ##############################
#   
# #gem_2041
# #load daily mean temp
# tavg <- brick(str_c('wy_projections/temp/tavg_HadGEM2-ES365_rcp85_2040.tif'))
# 
# #create vector of dates to use to extract data
# #create time series of years, months, days
# ts <- seq(mdy('01-01-2040'), mdy('12-31-2069'), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# #create index of which bands to extract
# y <- 2041
# index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))
# 
# #subset raster dataset
# tavg <- tavg[[index]]
# 
# #plot pirgd map
# plot(ras_gem_2041[[1]])
# 
# #load gdd values from click
# vals <- click(tavg)
# 
# #from tavg variables, GDD equals tavg - 273.15 - 5
# #set negatives to zero
# vals <- round(vals - 273.15 - 5, 2)
# vals[vals < 0] <- 0
# 
# #calc cumsum
# vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t
# 
# #set plotting
# par(mfrow = c(2, 2))
# 
# #calculate cumsum
# plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))
# 
# #reset par
# dev.off()
# 
# #gem_2071
# #load daily mean temp
# tavg <- brick(str_c('wy_projections/temp/tavg_HadGEM2-ES365_rcp85_2070.tif'))
# 
# #create vector of dates to use to extract data
# #create time series of years, months, days
# ts <- seq(mdy('01-01-2070'), mdy('12-31-2099'), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# #create index of which bands to extract
# y <- 2071
# index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))
# 
# #subset raster dataset
# tavg <- tavg[[index]]
# 
# #plot pirgd map
# plot(ras_gem_2071[[1]])
# 
# #load gdd values from click
# vals <- click(tavg)
# 
# #from tavg variables, GDD equals tavg - 273.15 - 5
# #set negatives to zero
# vals <- round(vals - 273.15 - 5, 2)
# vals[vals < 0] <- 0
# 
# #calc cumsum
# vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t
# 
# #set plotting
# par(mfrow = c(2, 2))
# 
# #calculate cumsum
# plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))
# 
# #reset par
# dev.off()
# 
# #inmcm_2041
# #load daily mean temp
# tavg <- brick(str_c('wy_projections/temp/tavg_inmcm4_rcp85_2040.tif'))
# 
# #create vector of dates to use to extract data
# #create time series of years, months, days
# ts <- seq(mdy('01-01-2040'), mdy('12-31-2069'), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# #create index of which bands to extract
# y <- 2041
# index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))
# 
# #subset raster dataset
# tavg <- tavg[[index]]
# 
# #plot pirgd map
# plot(ras_inmcm_2041[[1]])
# 
# #load gdd values from click
# vals <- click(tavg)
# 
# #from tavg variables, GDD equals tavg - 273.15 - 5
# #set negatives to zero
# vals <- round(vals - 273.15 - 5, 2)
# vals[vals < 0] <- 0
# 
# #calc cumsum
# vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t
# 
# #set plotting
# par(mfrow = c(2, 2))
# 
# #calculate cumsum
# plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))
# 
# #reset par
# dev.off()
# 
# #inmcm_2071
# #load daily mean temp
# tavg <- brick(str_c('wy_projections/temp/tavg_inmcm4_rcp85_2070.tif'))
# 
# #create vector of dates to use to extract data
# #create time series of years, months, days
# ts <- seq(mdy('01-01-2070'), mdy('12-31-2099'), by = "day") %>%
#   format(., format = '%Y%m%d')
# 
# #create index of which bands to extract
# y <- 2071
# index <- which(str_detect(ts, str_c(y, "01|", y, "02|", y, "03|", y, "04")))
# 
# #subset raster dataset
# tavg <- tavg[[index]]
# 
# #plot pirgd map
# plot(ras_inmcm_2071[[1]])
# 
# #load gdd values from click
# vals <- click(tavg)
# 
# #from tavg variables, GDD equals tavg - 273.15 - 5
# #set negatives to zero
# vals <- round(vals - 273.15 - 5, 2)
# vals[vals < 0] <- 0
# 
# #calc cumsum
# vals_csum <- apply(vals, 1, FUN = function(x) cumsum(x)) %>% t
# 
# #set plotting
# par(mfrow = c(2, 2))
# 
# #calculate cumsum
# plot(vals_csum[1,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SE (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[2,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'SW (low PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[3,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (high PIRGd)', ylim = c(0, max(vals_csum)))
# plot(vals_csum[4,], xlab = 'DOY', ylab = 'Cumulative GDD',
#      main = 'NW (mid PIRGd)', ylim = c(0, max(vals_csum)))
# 
# #reset par
# dev.off()
# 
# ################################
# ###LOOK AT PIRGd BY ELEVATION###
# ################################
# 
# #load pirgd and dem
# pirgd <- raster('dlc/mean_maxIRGdate_wy_laea_2001_2018.tif')
# dem <- raster('dem/srtm/dem_wy_laea.tif')   
# 
# #create dem categores
# dem_lo <- dem
# dem_lo[dem_lo >= 1645] <- NA
# 
# dem_mi <- dem
# dem_mi[dem_mi < 1645 | dem_mi >= 2298] <- NA
# 
# dem_hi <- dem
# dem_hi[dem_hi < 2298] <- NA
# 
# #load WY shapefile
# wy <- readOGR('reference/wyoming.shp')
# 
# #reproject to PIRGd
# wy <- spTransform(wy, crs(pirgd))
# 
# #create pirgd layers by elevation
# pirgd %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
#                main = "PIRGd at all elevations") 
# plot(wy, add = T)
# 
# pirgd %>% mask(., dem_lo) %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
#                                    main = "PIRGd below 1645 m")
# plot(wy, add = T)
# 
# pirgd %>% mask(., dem_mi) %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
#                                    main = 'PIRGd between 1645 and 2298 m')
# plot(wy, add = T)
# 
# pirgd %>% mask(., dem_hi) %>% plot(col = terrain.colors(length(seq(80, 220, by = 20))), breaks= seq(80, 220, by = 20),
#                                    main = "PIRGd above 2298 m")
# plot(wy, add = T)
# 
# #plot pirgd by elevation
# plot(pirgd[dem < 1645])
# 
# 
# 
# 
