#code to extract sampling grid for study
#based on a 16 km grid and moving points to avoid missing data

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(seegSDM)

#set wd
setwd('/Volumes/SSD/climate_effects')

#set rasteroptions to use a bit more memory
rasterOptions()
rasterOptions(maxmemory = 1e+10)

#load dlc data
pirgd <- stack('/Volumes/SSD/climate_effects/dlc/maxIRGdate_wy_laea_2000_2019.tif')

#we only want a single layer to use as a mask for sampling selection
mask <- pirgd[[1]]

#load treatments throughout WY
t1 <- readOGR('/Volumes/SSD/climate_effects/treatments/Completed_Treatments_December2019_FINAL.shp')
t2 <- readOGR('/Volumes/SSD/climate_effects/treatments/LTDL_Treatment_Polygons_wy.shp')
t3 <- readOGR('/Volumes/SSD/climate_effects/treatments/SWyRange_CompletedTreatments_thru2019.shp')

#transform to CRS of prigd
t1 <- t1 %>% spTransform(., crs(pirgd))
t2 <- t2 %>% spTransform(., crs(pirgd))
t3 <- t3 %>% spTransform(., crs(pirgd))

#mask pirgd based on treatments
mask <- mask %>% mask(., mask = t1, inverse = T) %>% 
  mask(., mask = t2, inverse = T) %>% 
  mask(., mask = t3, inverse = T)

#clean up
rm(t1, t2, t3)

#load fire data
fire <- raster('fire/burn_wy_laea_1984_2017.tif')

#change fire values to use as mask
fire[fire > 1980] <- 1

#apply to mask
mask <- mask %>% mask(., mask = fire, maskvalue = 1)
rm(fire)

#load landcover change (1 = no change)
lc_ch <- raster('landcover/landcover_change_wy_laea_2016.tif')

#change values to use as mask
lc_ch[lc_ch > 1] <- 2

#apply to mask
mask <- mask %>% mask(., mask = lc_ch, maskvalue = 2)
rm(lc_ch)

#load landcover
lc <- raster('landcover/five_class_landcover_wy_laea_2016.tif')

#apply to mask and keep lc for sample stratification
mask <- mask %>% mask(., mask = lc)

#load wy and plot to share
wy <- readOGR('reference/wyoming.shp') %>% spTransform(., crs(mask))
plot(mask, main = "Example PIRGd (DOY) throughout study area with disturbances masked")
plot(wy, add = T)

#load sampling grid at ~16 km
grid <- readOGR('reference/sampling_points.shp')
grid <- spTransform(grid, crs(mask))

#load additional disturbance layers and mask study area
#load geomac fire data
geomac <- readOGR('fire/US_HIST_FIRE_PERIMTRS_DD83/US_HIST_FIRE_PERIMTRS_DD83.shp')

#transform to CRS of prigd
geomac <- geomac %>% spTransform(., crs(mask))

#mask pirgd based on treatments
mask <- mask %>% mask(., mask = geomac, inverse = T) 

#clean up
rm(geomac)

#load tree loss layer
loss <- raster('forest/tree_loss_year_wy_laea_2000_2018.tif')

#change values to use as mask. zero is no loss
loss[loss > 0] <- 1

#apply to mask
mask <- mask %>% mask(., mask = loss, maskvalue = 1)

rm(loss)

#load snowmelt timing stack
snowmelt <- brick('/Volumes/SSD/climate_effects/snowmelt/snowmelt_timing_wy_laea_2001_2018.tif')

#calculate raster of number of years snowmelt is NA (which is snowmelt == 0)
snow_zero <- calc(snowmelt, fun = function(x) sum(x == 0, na.rm = T))

#remove points from sampling if more than half of years have missing snowmelt data
snow_zero[snow_zero <= 3] <- 1
snow_zero[snow_zero > 3] <- 0

#apply to mask
mask <- mask %>% mask(., mask = snow_zero, maskvalue = 0)

rm(snowmelt, snow_zero)

#save image
save.image('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/extract_sampling_grid.RData')

#load image
#load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/extract_sampling_grid.RData')

#calculate coordinates of closest non NA points
coords <- nearestLand(data.frame(x = grid@coords[,1], y = grid@coords[,2]), mask, max_distance = 30000)

#create new grid without NA values
grid2 <- grid
grid2@coords <- coords

#write mask and new points to file
writeRaster(mask, 'dlc/maxIRGdate_mask_disturbances.tif', format = "GTiff", overwrite = T)
writeOGR(grid2, 'reference', 'sampling_points_no_na', driver = 'ESRI Shapefile', overwrite = T)

#plot with points
plot(mask, main = 'Example PIRGd (DOY) throughout study area with disturbances 
     masked and adjusted sampling grid')
plot(wy, add = T)
plot(grid2, add = T)

#now we need to randomly sample the shrub and grasslands pixels to
#match the number of evergreen pixels (which has the least amount)

#lc codes
#lc[lc %in% shrub] <- 1
#lc[lc %in% herb] <- 2
#lc[lc %in% deci] <- 3
#lc[lc %in% ever] <- 4
#lc[lc %in% wetl] <- 5

#extract lc values at data points and reenter into grid2 data
dd <- raster::extract(lc, grid2)
grid2@data <- cbind(grid2@data, lc = dd)

#check number of shrub, herb, and evergreen points
#shrub
nrow(grid2@data[grid2@data$lc == 1,])

#herb
nrow(grid2@data[grid2@data$lc == 2,])

#evergreen
n_samples <- nrow(grid2@data[grid2@data$lc == 4,])

#check number of deciduous points, will keep separate
nrow(grid2@data[grid2@data$lc == 3,])

#create separate spdf's for each lc type we want
#sample number of points to match evergreen
set.seed(2)

#shrub
grid_shrub <- grid2[grid2@data$lc == 1, ]
grid_shrub <- grid_shrub[sample(1:length(grid_shrub), n_samples),]

#herb
grid_herb <- grid2[grid2@data$lc == 2, ]
grid_herb <- grid_herb[sample(1:length(grid_herb), n_samples),]

#evergreen
grid_ever <- grid2[grid2@data$lc == 4, ]

#union grid of sampling points
grid3 <- union(grid_shrub, grid_herb)
grid3 <- union(grid3, grid_ever)

#check that we have correct number of points for each class
check <- raster::extract(lc, grid3)
length(check[check == 1])
length(check[check == 2])
length(check[check == 4])

#write out grid of sampling points
writeOGR(grid3, 'reference', 'sampling_points_three_equal_classes', driver = 'ESRI Shapefile', overwrite = T)

#plot with points
plot(mask, main = 'Example PIRGd (DOY) throughout study area with disturbances 
     masked and adjusted sampling grid')
plot(wy, add = T)
plot(grid3, add = T)

