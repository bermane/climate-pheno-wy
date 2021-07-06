#this code compares fitted model outputs for SS against BIT data from the middle model
#to assess the spatial locations where the model is performing better/worse

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

########################
###RUN FINAL SS MODEL###
########################

#run model without rescaled values
m_ss <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + log(rain_mar_may) + 
              log(vpd_tavg_mean_jan_apr) + poly(pdsi_mar_apr_min, 2, raw = TRUE):log(rain_mar_may) +
              log(vpd_tavg_mean_jan_apr):log(rain_mar_may), 
            random = ~ 1 | id, data = dat2)

#create new df with id, year, and fitted values
vals <- data.frame(id = dat2$id, year = dat2$year, fitted = fitted(m_ss) %>% round(2))

#########################################
###ADD BIT VALUES FROM MIDDLE SCENARIO###
#########################################

#load BIT values
pred_rcp45 <- brick('wy_projections/ss/ss_ann_HadGEM2-ES365_rcp45_1999.tif') %>%
  .[[3:20]]

pred_rcp85 <- brick('wy_projections/ss/ss_ann_HadGEM2-ES365_rcp85_1999.tif') %>%
  .[[3:20]]

#change crs of sampling grid
grid_maca <- spTransform(grid, crs(pred_rcp45))

#extract pirgd vals
pred_rcp45 <- raster::extract(pred_rcp45, grid_maca)
pred_rcp85 <- raster::extract(pred_rcp85, grid_maca)

#allocate pred values
vals$pred_rcp45 <- as.numeric(NA)
vals$pred_rcp85 <- as.numeric(NA)

#loop through years
for(i in 1:length(unique(vals$year))){
  
  #load df of vals from year
  id_hold <- vals$id[vals$year == unique(vals$year)[i]] %>% as.numeric
  
  #loop through pixel ids
  for(j in 1:NROW(grid)){
    
    #if id exists in vals
    if((j %in% id_hold) == T){
      
      #load pred value
      vals$pred_rcp45[vals$id == j & vals$year == unique(vals$year)[i]] <- pred_rcp45[j, i]
      vals$pred_rcp85[vals$id == j & vals$year == unique(vals$year)[i]] <- pred_rcp85[j, i]
      
    }
  }
}

#calc difference between fitted and  BIT rcp 4.5
vals$resid <- vals$pred_rcp45 - vals$fitted

#calc average difference across years
vals_avg <- vals %>% group_by(id) %>% mutate(avg_resid = mean(resid, na.rm = T) %>% round(2)) %>%
  .[!duplicated(.$id),] %>% .[order(.$id),]

#set NaN values to 0
vals_avg$avg_resid[is.na(vals_avg$avg_resid)] <- 0

##########
###PLOT###
##########

#add mean residuals to grid of points
grid@data$resid <- vals_avg$avg_resid

#create positive and negative residual values
grid_pos <- grid[grid@data$resid >= 0,]
grid_neg <- grid[grid@data$resid < 0,]

#load wy shapefile
wy <- readOGR('/Volumes/SSD/climate_effects/reference/wyoming.shp')

#load lc
lc <- raster('/Volumes/SSD/climate_effects/landcover/five_class_landcover_wy_laea_2016.tif')
lc[lc == 3] <- NA
lc[lc == 5] <- NA
lc[lc == 4] <- 3

#set up output plot
jpeg('/Volumes/SSD/climate_effects/output/final_plots/resid_size_fitted_bit_ss.jpeg', width = 715*2, height = 495*2)

#plot grid points
plot(lc, legend=FALSE, col=c("coral3", "papayawhip", "forestgreen"), xaxt='n', yaxt='n',
     main = "Average SS Residual Size on Sampling Grid Fitted vs. BIT \n Black Pos Blue Neg") 
par(xpd=TRUE)
legend("bottom", legend=c("Shrub", "Herb", "Evergreen"),
       fill=c("coral3","papayawhip", "forestgreen"),horiz = TRUE, inset=-0.175)
plot(wy, add = T)
plot(grid, pch = 1, add = T, cex = .2)
plot(grid_pos, pch = 19, cex = abs(grid@data$resid), add = T, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.75))
plot(grid_neg, pch = 19, cex = abs(grid@data$resid), add = T, 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.75))

#dev.off
dev.off()

