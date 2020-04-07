#code to organize, crop, and output finalized covariates for analysis

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(dynatopmodel)

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

