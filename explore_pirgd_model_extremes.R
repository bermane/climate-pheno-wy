#this code generates various summary statistics to explore issues
#with extreme values and how they are projected into the future.
#it is exploratory in nature and will direct the process used 
#to build the final model

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(ggplot2)
library(qwraps2)

#set wd
setwd('/Volumes/SSD/climate_effects')

options(qwraps2_markup = "markdown")

#######################################################
###SUMMARY STATISTICS CONTEMPORARY MODEL CO-VARIATES###
#######################################################

#load contemporary pirgd model data
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_10_18.RData')

#remove missing and unwanted values

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

#any other missing points?
#may need to remove the variables that include snowmelt and pirgd...
#16% of samples have PIRGd before snowmelt...
sum(is.na(dat2))

#keep only the co-variate columns we want
dat2 <- dat2[,colnames(dat2) %in% c('id', 'year', 'pirgd', 'lc', 'pdsi_mar_apr_min',
                                    'snow_oct_apr', 'gdd_jan_apr', 'rain_mar_may')]

#create summary statistics table of 4 co-variates
#create output list qwraps
summary_contemp <-
  list("Rain Mar-May (mm)" =
         list("min"       = ~ min(rain_mar_may),
              "max"       = ~ max(rain_mar_may),
              "mean (sd)" = ~ qwraps2::mean_sd(rain_mar_may),
              "mean (95% ci)" = ~ qwraps2::mean_ci(rain_mar_may),
              "median (iqr)"    = ~ qwraps2::median_iqr(rain_mar_may)),
       "Snow Oct-Apr (mm)" =
         list("min"       = ~ min(snow_oct_apr),
              "max"       = ~ max(snow_oct_apr),
              "mean (sd)" = ~ qwraps2::mean_sd(snow_oct_apr),
              "mean (95% ci)" = ~ qwraps2::mean_ci(snow_oct_apr),
              "median (iqr)"    = ~ qwraps2::median_iqr(snow_oct_apr)),
       "GDD Jan-Apr (deg C)" =
         list("min"       = ~ min(gdd_jan_apr),
              "max"       = ~ max(gdd_jan_apr),
              "mean (sd)" = ~ qwraps2::mean_sd(gdd_jan_apr),
              "mean (95% ci)" = ~ qwraps2::mean_ci(gdd_jan_apr),
              "median (iqr)"    = ~ qwraps2::median_iqr(gdd_jan_apr)),
       "Minimum scPDSI Mar-Apr" =
         list("min"       = ~ min(pdsi_mar_apr_min),
              "max"       = ~ max(pdsi_mar_apr_min),
              "mean (sd)" = ~ qwraps2::mean_sd(pdsi_mar_apr_min),
              "mean (95% ci)" = ~ qwraps2::mean_ci(pdsi_mar_apr_min),
              "median (iqr)"    = ~ qwraps2::median_iqr(pdsi_mar_apr_min))
  )

#create summary table
summary <- summary_table(dat2, summary_contemp)
summary
