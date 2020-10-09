#code to organize, crop, and output finalized covariates for analysis

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(dynatopmodel)
library(foreach)
library(spatial.tools)
library(gdalUtils)
library(lubridate)
library(RCurl)

#set wd
setwd('/Volumes/SSD/climate_effects')

#set raster options to use a bit more memory
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

#load standard deviation stack
#sd <- stack('/Volumes/SSD/climate_effects/dlc/SD_IRGMaxDate.grd')

#reproject to eMODIS grid
#sd <- sd %>% projectRaster(., pirgd)
#writeRaster(se, filename = 'dlc/sd_maxIRGdate_wy_laea_2000_2019.tif', format = "GTiff")
#rm(sd)

#load standard deviation stack
sd <- stack('/Volumes/SSD/climate_effects/dlc/sd_maxIRGdate_wy_laea_2000_2019.tif')

#we only want 2001-2018
sd <- sd[[2:19]]

#calc mean ss for each pixel
sd_med <- calc(sd, function(x) {median(x,na.rm = T)})

#calc number of na values in annual layers
sd_na <- calc(sd, function(x) {sum(is.na(x))})

#only save pirgd mean value if less than half of years are missing
#we have 19years total
sd_med[sd_na > 9] <- NA

#write and clean up
writeRaster(sd_med, filename = 'dlc/median_sd_maxIRGdate_wy_laea_2001_2018.tif', format = "GTiff")
rm(sd_med, sd_na)

#calc number of years sd > 5
sd_5 <- sd
sd_5[sd_5 < 5] <- 0
sd_5[sd_5 >= 5] <- 1

sd_5_sum <- calc(sd_5, function(x) {sum(x, na.rm = T)})

#write and clean up
writeRaster(sd_5_sum, filename = 'dlc/sd_sum_5_maxIRGdate_wy_laea_2001_2018.tif', format = "GTiff")
rm(sd, sd_5, sd_5_sum)

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

#sort so in order by month
pdsi_f <- str_sort(pdsi_f, numeric = T)

#project pirgd_extent so we can used it to crop
pirgd_ex_drought <- pirgd_extent %>% projectRaster(., crs = crs(raster(pdsi_f[1])))

for(i in 1:length(pdsi_f)){
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
for(i in 1:20){
  pdsi <- stack(pdsi_f, bands = i)
  writeRaster(pdsi, filename = str_c('drought/pdsi_wy_laea_', 1999 + i, '.tif'),
              format = 'GTiff', overwrite = T)
  rm(pdsi)
}

#clean up
rm(pdsi_f)

#remove bymonth files
do.call(file.remove, list(list.files('drought/bymonth', full.names = TRUE)))

###########################
###ANNUAL HERBACEOUS B&W###
###########################

#load annherb data
herb_files <- list.files(path = 'annual_herbaceous/raw', pattern = '.img$', full.names = T)
herb <- stack(herb_files)
rm(herb_files)

#project pirgd_extent so we can use it to crop
pirgd_ex_herb <- pirgd_extent %>% 
  projectRaster(., crs = crs(herb))

#crop, reproject, crop
herb <- herb %>% crop(., pirgd_ex_herb)  %>% 
  projectRaster(., pirgd, method = "bilinear") %>% crop(., pirgd)

herb <- round(herb)

writeRaster(herb, filename = 'annual_herbaceous/ann_herb_wy_laea_2000_2016_bw.tif',
            format = "GTiff")

#clean up
rm(herb)

#############
###DAYMET####
#############

#code to mosaic together daymet tiles and crop to study area
###RUN FOR TEMP AND PRECIP AND ALL YEARS
#set initial variables
#clim <- c("tmax", "prcp", "tmin")
clim <- "vp"
years <- 2000:2019

#project pirgd_extent so we can use it to crop
pirgd_ex_dm <- pirgd_extent %>% 
  projectRaster(., crs = crs(raster("vp/mosaic/vp_wy_2000.tif")))

#run for all variables
for(vari in clim){
  
    #run for all years  
    for(yr in years){
      
        #load mosaicked tile
        mos <- brick(str_c(vari, '/mosaic/vp_wy_', yr, '.tif'))
        
        #disag, reproject, crop
        mos <- mos %>% crop(., pirgd_ex_dm) %>% raster::disaggregate(., fact = 4) %>% 
          projectRaster(., pirgd, method = "bilinear") %>% crop(., pirgd)
        
        #write to disk
        writeRaster(mos, paste0(vari ,'/', vari, '_wy_laea_', yr, '.tif'), format = "GTiff")
        
        #clean up
        rm(mos)
        removeTmpFiles(h = 0.000000001)
      }
}


#########
###RAP###
#########

#need to find the extent in WGS84 to subset and download in terminal
#reproj pirgd_extent to wgs84 and try those coords
#pirgd_extent_rap <- pirgd_extent %>% 
#  projectRaster(., crs = CRS('+init=epsg:4326'))

#extent(pirgd_extent_rap)

#set years to process
years = 2001:2019

#loop through years to process
for(yr in years){
  
  #import layer
  rap <- stack(str_c('/Volumes/SSD/climate_effects/rap/raw/rap_', yr, '.tif'))
  
  #only keep layers we want
  rap <- rap[[1:6]]
  
  #aggregate, reproject, crop
  rap <- rap %>% raster::aggregate(., fact = 8) %>% 
    projectRaster(., pirgd, method = "bilinear") %>% crop(., pirgd) %>% round
  
  #write to disk
  writeRaster(rap, str_c('rap/rap_wy_laea_', yr, '.tif'), format = "GTiff")
  
  #clean up
  rm(rap, pirgd_extent_rap)
  removeTmpFiles(h = 0.000000001)
}

##################################
###FOREST LOSS AND COVER HANSEN###
##################################

#load forest loss rasters
loss1 <- raster('/Volumes/SSD/climate_effects/forest/Hansen_GFC-2018-v1.6_lossyear_50N_120W.tif')
loss2 <- raster('/Volumes/SSD/climate_effects/forest/Hansen_GFC-2018-v1.6_lossyear_50N_110W.tif')

#merge together
loss <- merge(loss1, loss2)

#project pirgd_extent so we can use it to crop
pirgd_ex_loss <- pirgd_extent %>% 
  projectRaster(., crs = crs(loss))

#crop, reproject, crop
loss <- loss %>% crop(., pirgd_ex_loss)  %>% 
  raster::aggregate(., fact = 8, fun = modal) %>% 
  projectRaster(., pirgd, method = "ngb") %>% 
  crop(., pirgd)

writeRaster(loss, filename = 'forest/tree_loss_year_wy_laea_2000_2018.tif',
            format = "GTiff")

#clean up
rm(loss)

#load tree cover rasters
tc1 <- raster('forest/Hansen_GFC-2018-v1.6_treecover2000_50N_120W.tif')
tc2 <- raster('forest/Hansen_GFC-2018-v1.6_treecover2000_50N_110W.tif')

#merge together
tc <- merge(tc1, tc2)

#crop, reproject, crop
tc <- tc %>% crop(., pirgd_ex_loss)  %>% 
  raster::aggregate(., fact = 8) %>% 
  projectRaster(., pirgd, method = "bilinear") %>% 
  crop(., pirgd)

writeRaster(tc, filename = 'forest/tree_cover_year_wy_laea_2000.tif',
            format = "GTiff")

#clean up
rm(tc, pirgd_ex_loss)

############
###SNODAS###
############

#######################
###DOWNLOAD RAW DATA###
#######################

#set download directory
dir <- '/Volumes/SSD/goat_surveys/snodas/raw'

for(yr in 2001:2005){
  url<- str_c("ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/", yr, "/")
  mon <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE) #reading filenames from ftp-server
  mon <- strsplit(mon, "\n")
  mon = unlist(mon)
  mon <- mon %>% str_subset(pattern = "_")
  
  for(m in mon){
    url2 <- str_c("ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/", yr, "/", m, "/")
    filenames <- getURL(url2, ftp.use.epsv = FALSE, dirlistonly = TRUE) #reading filenames from ftp-server
    filenames <- strsplit(filenames, "\n")
    filenames = unlist(filenames)
    filenames <- filenames %>% str_subset(pattern = ".tar")
    
    for (filename in filenames) {
      download.file(paste(url2, filename, sep = ""), paste(dir, "/", filename,
                                                           sep = ""))
    }
    rm(url2, filenames)
  }
  rm(url, mon)
}

#######################
###ORGANIZE RAW DATA###
#######################

#denote where raw data is
raw <- '/Volumes/SSD/goat_surveys/snodas/raw'

#create temp dir and output files dir
dir.create('snodas/temp')
dir.create('snodas/extracted')

#get files in dir
files <- list.files(raw, '.tar', full.names = T)

#find leap year files
leap <- str_which(files, '0229')

#remove leap year files
file.remove(files[leap])

#get files in dir
files <- list.files(raw, '.tar', full.names = T)

#check to see if any files are missing
#extract dates from file names and parse as Date object
check <- str_extract(files, '[:digit:]{8}')
check <- ymd(check)

#create sequence along timeseries from 2001 to end of 2019
dates <- seq(ymd("20010101"), ymd(20191231), by = 1)

#print which are missing
#seems these files just don't exist so when I stack the data i'll
#just leave whole bands as NA
dates[!(dates %in% check)]

#we only want oct to april so only extract those months
months <- c(1, 2, 3, 4, 10, 11, 12)
files <- files[month(check) %in% months]

#set to run for three variables
vari <- matrix(c('swe', 'snowdepth', 'snowprcp',
                 'us_ssmv11034', 'us_ssmv11036', 'us_ssmv01025SlL01'), ncol = 2)

#run loop to extract files
for(file in files){
  #untar main file
  untar(file, exdir = str_c(getwd(), raw, '/temp'))
  
  #loop through, unzip, crop, and write .tif for 3 variables
  for(i in 1:nrow(vari)){
    
    #find file for variable
    f <- list.files(str_c(getwd(), raw, '/temp'), pattern = glob2rx(str_c(vari[i,2], '*.dat.gz')), full.names = T)
    f_out <- str_c(vari[i,1], '_', str_extract(file, '[:digit:]{8}'), '.dat')
    gunzip(f[1], destname = str_c(getwd(), raw, '/temp/', f_out))
    
    #create hdr
    #different hdr for before and after oct 1 2013
    if(ymd(str_extract(file, '[:digit:]{8}')) < "2013-10-01"){
      hdr <- c('nrows 3351',
               'ncols 6935',
               'nbands 1',
               'nbits 16',
               'pixeltype signedint',
               'byteorder M',
               'layout dat',
               'ulxmap -124.729583333331703',
               'ulymap 52.871249516804028',
               'xdim 0.00833333333',
               'ydim 0.00833333333')
    }else{
      hdr <- c('nrows 3351',
               'ncols 6935',
               'nbands 1',
               'nbits 16',
               'pixeltype signedint',
               'byteorder M',
               'layout dat',
               'ulxmap -124.733333333333',
               'ulymap 52.8749999999999',
               'xdim 0.00833333333',
               'ydim 0.00833333333')
    }
    
    #write hdr to disk
    writeLines(hdr, str_c(getwd(), raw, '/temp/', str_replace(f_out, '.dat', '.hdr')))
    
    #load in dat file
    sno <- raster(str_c(getwd(), raw, '/temp/', f_out))
    
    #set crs
    crs(sno) <- CRS("+init=epsg:4326")
    
    #reproj gnp outline
    gnp_sno <- spTransform(gnp, crs(sno))
    
    #crop sno and write to extracted folder
    sno <- crop(sno, gnp_sno, filename = str_c(getwd(), raw, '/extracted/', str_replace(f_out, '.dat', '.tif')),
                format = 'GTiff', datatype = 'INT2S', NAflag = -9999)
    
    #clean up
    rm(f, f_out, hdr, gnp_sno, sno)
  }
  #clear temp folder
  invisible(do.call(file.remove, list(list.files(str_c(getwd(), raw, '/temp'), full.names = TRUE))))
}

#######################################
###HOMER ANNUAL SHRUBLAND COMPONENTS###
#######################################

#load test layers
test <- raster('homer_annual/raw/nlcd_herb_2001_mos_v1_gAhIbnWn4xjidvl8NSjM.tiff')
test2 <- raster('homer_annual/raw/nlcd_herb_2001_mos_v1_s9BgbK9qx4Th9eAFLxUc.tiff')

#project pirgd_extent so we can used it to crop
pirgd_ex_homer <- pirgd_extent %>% projectRaster(., crs = crs(test))

#loop through sage, herb and shrub datasets
for(var in c('sage', 'herb', 'shrub')){
  
  #load all sage files scene 1
  s1 <- Sys.glob(str_c('homer_annual/raw/nlcd_', var, '*gAhIbnWn4xjidvl8NSjM.tiff'))
  s1 <- stack(s1)
  
  #load all sage files scene 2
  s2 <- Sys.glob(str_c('homer_annual/raw/nlcd_', var, '*s9BgbK9qx4Th9eAFLxUc.tiff'))
  s2 <- stack(s2)
  
  #mosaic datasets, initial crop, aggregate
  scene <- merge(s1, s2) %>% crop(., pirgd_ex_homer) %>% aggregate(., fact = 8, fun = modal)
  
  #rm
  rm(s1, s2)
  
  #convert na data to NA
  scene[scene == 101] <- NA
  scene[scene == 102] <- NA
  
  #reproject, crop
  scene <- scene %>% projectRaster(., pirgd, method = "bilinear") %>% round %>% crop(., pirgd)
  
  #write to disk
  writeRaster(scene, filename = str_c('homer_annual/', var, '_wy_laea_2001_2018_no_2012.tif'), format = "GTiff")
  rm(scene)
  
  removeTmpFiles(h = 0.00000001)
}
