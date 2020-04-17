#code to organize, crop, and output finalized covariates for analysis

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(dynatopmodel)
library(foreach)
library(spatial.tools)

#set wd
setwd('/Volumes/SSD/climate_effects')

#set rasteroptions to use a bit more memory
rasterOptions()
rasterOptions(maxmemory = 1e+10)

#load dlc data, reproject to eMODIS grid since span UTM zones, output
#pirgd <- stack('dlc/maxIRGdate.grd')

#emodis_grid <- raster('emodis/pirgd/pirgd_mean_wy_2002_2018.tif')

#pirgd <- pirgd %>% projectRaster(., emodis_grid) %>% round
#writeRaster(pirgd, filename = 'dlc/maxIRGdate_wy_laea_2000_2019.tif', format = "GTiff")

#rm(emodis_grid)

pirgd <- stack('/Volumes/SSD/climate_effects/dlc/maxIRGdate_wy_laea_2000_2019.tif')

#we only want 2001-2018
pirgd <- pirgd[[2:19]]

#create single buffered layer that we can reproject and use to pre-crop before processing
pirgd_extent <- pirgd[[1]]
extent(pirgd_extent) <- extent(pirgd_extent) + 5000

#####################
###DEM AND TERRAIN###
#####################

#list dem files
dem_files <- list.files('dem/srtm/raw', full.names = T)

#mosaic rasters
dem_list <- lapply(1:length(dem_files),function(x) raster(dem_files[x]))
dem_list$fun <- mean

dem_mosaic <- do.call(mosaic, dem_list)

#aggregate, reproject, crop, write
dem_ag <- dem_mosaic %>% aggregate(fact = 8)
dem_ag <- dem_ag %>% projectRaster(., pirgd) %>% crop(., pirgd)
writeRaster(dem_ag, filename = 'dem/srtm/dem_wy_laea.tif', format = "GTiff")

#create twi, have to reproj to UTM and back
twi <- upslope.area(dem_ag, log = TRUE, atb = TRUE, deg = 0.1, fill.sinks = TRUE)
writeRaster(twi[[2]], filename = "dem/srtm/twi_wy_laea.tif", format = "GTiff")

rm(dem_list, dem_ag, dem_mosaic, twi)

################
###SNOW COVER###
################

#list snowmelt timing files
smt_files <- Sys.glob('snowmelt/Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_20*.tif', dirmark = T)

#load snowmelt timing maps into raster stack
smt <- stack(smt_files)

#disag, reproj, crop, write
smt <- smt %>% raster::disaggregate(., fact = 2)
smt <- smt %>% projectRaster(., pirgd, method = "bilinear") %>% crop(., pirgd)
smt <- smt %>% round
writeRaster(smt, filename = 'snowmelt/snowmelt_timing_wy_laea_2001_2018.tif', format = "GTiff", overwrite = T)

rm(smt, smt_files)

###############
###LANDCOVER###
###############

#load 2016 land cover map
lc <- raster('landcover/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img')

#project pirgd_extent so we can used it to crop
pirgd_extent_lc <- pirgd_extent %>% projectRaster(., crs = crs(lc))

#crop lc
lc <- lc %>% crop(., pirgd_extent_lc)
rm(pirgd_extent_lc)

#aggregate lc using modal
lc2 <- lc %>% aggregate(., fact = 8, fun = modal)

#reproject and crop land cover map, write
lc2 <- lc2 %>% projectRaster(., pirgd, method = "ngb") %>% crop(.,pirgd)
writeRaster(lc2, filename = 'landcover/landcover_wy_laea_2016.tif', format = "GTiff")

rm(lc, lc2)

#reload landcover to set classes
lc <- raster('/Volumes/SSD/climate_effects/landcover/landcover_wy_laea_2016.tif')

#segment values we want to keep
lc_keep <- c(41, 42, 43, #deci, evergreen, mixed forests
             52, #shrub
             71, #herbaceous
             90, 95) #woody wetlands, herbaceous wetlands

shrub <- 52
herb <- 71
deci <- c(41, 43)
ever <- 42
wetl <- c(90, 95)

#remove other values
lc[!(lc %in% lc_keep)] <- NA

#reset lc values
lc[lc %in% shrub] <- 1
lc[lc %in% herb] <- 2
lc[lc %in% deci] <- 3
lc[lc %in% ever] <- 4
lc[lc %in% wetl] <- 5

#write
writeRaster(lc, filename = 'landcover/five_class_landcover_wy_laea_2016.tif', format = "GTiff")

rm(lc)

#load 2016 land cover change map
lc <- raster('landcover/NLCD_Land_Cover_Change_Index_L48_20190424/NLCD_Land_Cover_Change_Index_L48_20190424.img')

#project pirgd_extent so we can used it to crop
pirgd_extent_lc <- pirgd_extent %>% projectRaster(., crs = crs(lc))

#crop lc
lc <- lc %>% crop(., pirgd_extent_lc)
rm(pirgd_extent_lc)

#aggregate lc using modal
lc2 <- lc %>% aggregate(., fact = 8, fun = modal)

#reproject and crop land cover map, write
lc2 <- lc2 %>% projectRaster(., pirgd, method = "ngb") %>% crop(.,pirgd)
writeRaster(lc2, filename = 'landcover/landcover_change_wy_laea_2016.tif', format = "GTiff")

rm(lc, lc2)

#############
###DAYMET####
#############

#code to mosaic together daymet tiles and crop to study area
###RUN FOR TEMP AND PRECIP AND ALL YEARS
#set initial variables
#clim <- c("tmax", "prcp", "tmin")
clim <- c("tmin")
years <- 2009:2019
#years <- 2019

#project pirgd_extent so we can used it to crop
pirgd_extent_dm <- pirgd_extent %>% projectRaster(., crs = CRS("+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +lat_1=25 +ellps=WGS84 +lat_2=45"))

#run for all variables
for(vari in clim){
  
  #set vari wd
  setwd(paste0('/Volumes/SSD/climate_effects/',vari,'/raw'))
  
  #load all raw tile file names
  vari_f <- list.files()
  
  #check for any duplicate files
  if(anyDuplicated(vari_f) != 0){print(paste('Duplicate files!', vari, 'not processed!!!'))
    } else{print(paste('No duplicate files for', vari, '!!! Processing!!!'))
    
    #run for all years  
    for(yr in years){
      
      #test for first year
      #load annual files
      vari_yr <- grep(yr, vari_f, value = T)
      
      #check for all 12 tiles
      if(length(vari_yr) != 12){print(paste('Not all tiles present!', vari, 'for', yr, 'not processed!!!'))}
      else{
        print(paste('All tiles present!', vari, 'for', yr, 'processing!!!'))
        
        #stack rasters in list
        ras_list <- lapply(1:length(vari_yr), function(x){stack(vari_yr[x])})
        
        #set parameters for mosaic function
        names(ras_list)[1:2] <- c('x', 'y')
        ras_list$fun <- mean
        #ras_list$na.rm <- TRUE
        
        #mosaic raster stacks
        mos <- do.call(mosaic, ras_list)
        
        #disag, reproject, crop
        mos <- mos %>% crop(., pirgd_extent_dm) %>% raster::disaggregate(., fact = 4) %>% 
          projectRaster(., pirgd, method = "bilinear") %>% crop(., pirgd)
        
        #create integer values with two decimals to save writing space
        mos <- mos*100 %>% round
        
        #write to disk
        writeRaster(mos, paste0('/Volumes/SSD/climate_effects/', vari ,'/mosaic/', vari, '_wy_laea', yr, '.tif'), format = "GTiff",
                    datatype = "INT2S", NAflag = -32767)
        
        rm(mos, ras_list, vari_yr)
        removeTmpFiles(h = 0.000000001)
      }
    }
  }
  rm(vari_f) 
}

##########################
###SHRUBLAND COMPONENTS###
##########################

#load sagebrush layer
sage <- raster('shrublands/sagebrush/raw/NLCD_2016_Sagebrush_Shrubland_Fractional_Component_20191021.img')

#project pirgd_extent so we can used it to crop
pirgd_ex_sage <- pirgd_extent %>% projectRaster(., crs = crs(sage))

#initial crop and aggregate
sage <- sage %>% crop(., pirgd_ex_sage) %>% aggregate(., fact = 8, fun = modal)

#convert na data to NA
sage[sage == 101] <- NA
sage[sage == 102] <- NA

#reproject, crop
sage <- sage %>% projectRaster(., pirgd, method = "bilinear") %>% round %>% crop(., pirgd)

#write to disk
writeRaster(sage, filename = 'shrublands/sagebrush/sagebrush_wy_laea_2016.tif', format = "GTiff")
rm(sage)

#do two more times for annual herbaceous and shrub
#load annual herb layer
sage <- raster('shrublands/herbaceous/raw/NLCD_2016_Annual_Herb_Shrubland_Fractional_Component_20191021.img')

#initial crop and aggregate
sage <- sage %>% crop(., pirgd_ex_sage) %>% aggregate(., fact = 8, fun = modal)

#convert na data to NA
sage[sage == 101] <- NA
sage[sage == 102] <- NA

#reproject, crop
sage <- sage %>% projectRaster(., pirgd, method = "bilinear") %>% round %>% crop(., pirgd)

#write to disk
writeRaster(sage, filename = 'shrublands/herbaceous/ann_herbaceous_wy_laea_2016.tif', format = "GTiff")
rm(sage)

#load shrub layer
sage <- raster('shrublands/shrub/raw/NLCD_2016_Shrub_Shrubland_Fractional_Component_20191021.img')

#initial crop and aggregate
sage <- sage %>% crop(., pirgd_ex_sage) %>% aggregate(., fact = 8, fun = modal)

#convert na data to NA
sage[sage == 101] <- NA
sage[sage == 102] <- NA

#reproject, crop
sage <- sage %>% projectRaster(., pirgd, method = "bilinear") %>% round %>% crop(., pirgd)

#write to disk
writeRaster(sage, filename = 'shrublands/shrub/shrub_wy_laea_2016.tif', format = "GTiff")
rm(sage, pirgd_ex_sage)

######################################
###PIRGd/Spring Scale Summary Stats###
######################################

#load spring scale stack
#ss <- stack('/Volumes/SSD/climate_effects/dlc/springScale.grd')

#reproject to eMODIS grid
#ss <- ss %>% projectRaster(., pirgd)
#writeRaster(ss, filename = 'dlc/springScale_wy_laea_2000_2019.tif', format = "GTiff")

#load spring scale stack
ss <- stack('/Volumes/SSD/climate_effects/dlc/springScale_wy_laea_2000_2019.tif')

#we only want 2001-2018
ss <- ss[[2:19]]

#calc mean ss for each pixel
ss_m <- calc(ss, function(x) {mean(x,na.rm = T)})

#calc number of na values in annual layers
ss_na <- calc(ss, function(x) {sum(is.na(x))})

#only save pirgd mean value if less than half of years are missing
#we have 19years total
ss_m[ss_na > 9] <- NA
writeRaster(ss_m, filename = 'dlc/mean_springScale_wy_laea_2001_2018.tif', format = "GTiff")

#calc variance ss for each pixel
ss_var <- calc(ss, function(x) {var(x,na.rm = T)})
ss_var[ss_na > 9] <- NA
writeRaster(ss_var, filename = 'dlc/variance_springScale_wy_laea_2001_2018.tif', format = "GTiff")

#calc mean pirgd for each pixel
pirgd_m <- calc(pirgd, function(x) {mean(x,na.rm = T)}) %>% round

#calc variance pirgd for each pixel
pirgd_var <- calc(pirgd, function(x) {var(x,na.rm = T)}) %>% round

#calc number of na values in annual layers
pirgd_na <- calc(pirgd, function(x) {sum(is.na(x))})

#remove na and write
pirgd_m[pirgd_na > 9] <- NA
writeRaster(pirgd_m, filename = 'dlc/mean_maxIRGdate_wy_laea_2001_2018.tif', format = "GTiff")

pirgd_var[pirgd_na > 9] <- NA
writeRaster(pirgd_var, filename = 'dlc/variance_maxIRGdate_wy_laea_2001_2018.tif', format = "GTiff")

########################
###BURN SEVERITY DATA###
########################

#group together burned/not burned classes
burn_class <- c(2, 3, 4)
nburn_class <- c(1, 5, 6)

#find all burn files
burn_f <- list.files(path = 'fire', pattern = ".zip", recursive = T, full.names = T)

#unzip all
for(i in burn_f) unzip(i, exdir = 'fire/raw')

#find all tif files
burn_f <- list.files(path = 'fire/raw', pattern = '.tif$', full.names = T)

#load burn rasters into list, issue with extent
burn <- raster(burn_f[1])

#project pirgd_extent so we can used it to crop
pirgd_ex_burn <- pirgd_extent %>% projectRaster(., crs = crs(burn[[1]]))

#initial crop, reclassify, ag, reproj, crop
start_yr <- 1984

for(i in 1:length(burn_f)){
  #load raster
  burn <- raster(burn_f[i])
  
  #initial crop
  burn <- crop(burn, pirgd_ex_burn)
  
  #reclassify
  burn[burn %in% nburn_class] <- NA
  burn[burn %in% burn_class] <- 1
  
  #agg, reproject, and final crop
  burn <- burn %>% raster::aggregate(., fact = 8, fun = modal) %>% projectRaster(., pirgd, method = "ngb") %>% crop(., pirgd)
  
  #write annual file to disk
  writeRaster(burn, filename = str_c('fire/byyear/burn_', start_yr + i - 1, '.tif'), format = "GTiff")
  
  #clean up space
  rm(burn)
  removeTmpFiles(h = 0.0000001)
}

#clean up
rm(burn_class, nburn_class, burn_f, pirgd_ex_burn)

#now we want to convert into a single raster layer 
#with the yearly values representing the burns

#load filelist
burn_f <- list.files(path = 'fire/byyear', pattern = '.tif$', full.names = T) %>%
  str_sort(numeric = T)

#list all years of data
burn_yr <- burn_f %>% str_extract(pattern = "[:digit:]{4}") %>% as.numeric

#loop through all files
for(i in 1:length(burn_f)){
  
  #load raster
  burn <- raster(burn_f[i])
  
  #if first file, create output file with all NA values
  if(i == 1) {burn_out <- burn; burn_out[] <- NA}
  
  #fill output file with values from each year
  burn_out[burn == 1] <- burn_yr[i]
  
  #clean up
  rm(burn)
}

#write burn raster
writeRaster(burn_out, filename = 'fire/burn_wy_laea_1984_2017.tif', format = 'GTiff')

#clean up
rm(burn_out, burn_f, burn_yr, i)

##########
###PDSI###
##########

#if calculating it ourselves from Daymet data...

#install.packages('SPEI')
#install.packages('scPDSI')

#need mean temp
#need evapo
#need latitude

#use thornthwaite to calc evapo
#use pdsi to calc scPDSI

#code to organize preprocessed scPDSI from Westwide Drought Tracker
pdsi_f <- list.files(path = 'drought/raw', pattern = '.nc$', full.names = T)

#project pirgd_extent so we can used it to crop
pirgd_ex_drought <- pirgd_extent %>% projectRaster(., crs = crs(raster(pdsi_f[1])))

for(i in 8:length(pdsi_f)){
  #load raster
  pdsi <- stack(pdsi_f[i])
  
  #keep bands we want: years 2000-2019
  pdsi <- pdsi[[106:125]]
  
  #disag, reproject, crop
  pdsi <- pdsi %>% crop(., pirgd_ex_drought) %>% raster::disaggregate(., fact = 16) %>% 
    projectRaster(., pirgd, method = "bilinear") %>% crop(., pirgd)
  
  #write to disk
  writeRaster(pdsi, filename = str_c('drought/bymonth/pdsi_2000_2019_', i, '.tif'), format = 'GTiff')
  
  #clean up
  rm(pdsi)
  removeTmpFiles(h = 0.000000001)
}

#clean up
rm(pdsi_f, pirgd_ex_drought)

#files are now organized by month. change them to by year.
#load files and make sure in correct order
pdsi_f <- list.files(path = 'drought/bymonth', pattern = '.tif$', full.names = T) %>%
  str_sort(numeric = T)

#load in files by band (month) and write out
for(i in 13:20){
  pdsi <- stack(pdsi_f, bands = i)
  writeRaster(pdsi, filename = str_c('drought/pdsi_wy_laea_', 1999 + i, '.tif'),
              format = 'GTiff')
  rm(pdsi)
}

#clean up
rm(pdsi_f)

#remove bymonth files
do.call(file.remove, list(list.files('drought/bymonth', full.names = TRUE)))



