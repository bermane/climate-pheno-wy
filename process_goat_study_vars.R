#################
###DESCRIPTION###
#################

#This code downloads and processes temperature and moisture variables from
#daymet and snodas for montana goat study surveys and area quadrants. There
#are 5 variables calculated:
#1. Precip.Neo - total precipitation from may 15  - june 15 in study area quadrants (from Daymet data)
#2. Temp.Neo- average temp from may 15 - june 15 in study area quadrants (from Daymet data)
#3. Quad.TempWint - average temp from oct-apr in study area quadrants (from Daymet data)
#4. Temp.Hour - hourly temperatures corresponding to goat survey dates spreadsheet (from Daymet data)
#5. Quad.SWE - average swe from oct-apr in study area quadrants (from SNODAS data)

#Authored by Ethan Berman, November 2020

###############################
###SET KEY VARIABLES HERE!!!###
###############################

#CHANGE VALUES HERE!!!
folder <- '/Users/Ediz/Desktop/goat_test' #main directory to process data

start_yr <- 2008 #first year you want to process data

end_yr <- 2019 #last year you want to process data. Must be at least start_yr + 1 because variables
#calculate over winter from Oct to Apr

daymet_tiles <- c(12453, 12454) #daymet tiles to download. default is 2 tiles needed to cover goat study

daymet_vars <- c('tmax', 'tmin', 'prcp') #default daymet variables needed.

gnp_file <- '/Volumes/SSD/goat_surveys/reference/boundary2003.shp' #location of shapefile with gnp boundary

quad_file <- '/Volumes/SSD/goat_surveys/reference/gnp_quad.shp' #location of shapefile with four quadrants

occu_file <- '/Volumes/SSD/goat_surveys/reference/Predicted Occupancy_for_Tab.tif' #location of occupancy file
#to mask values where occupancy is low

surv_file <- '/Volumes/SSD/goat_surveys/reference/GoatSurveyTimes_MJYedit.csv' #if using different/updated file
#ensure it is the same format as the file listed here!!! and that you are processing all the necessary years
#to fill the table!!!

viewshed_file <- '/Volumes/SSD/goat_surveys/reference/Viewsheds_clipped_by_lakes.shp' #shapefile of the viewsheds that match the survey locations
#the names don't match the survey data perfectly and the code below deals with the discrepancies. These may change
#over time/with no data. Something to be aware of.

out_file_daymet_vars <- '/Users/Ediz/Desktop/goat_test/daymet_vars.csv' #output file name for dayment variables
#in quadrants

out_file_temp_hour <- '/Users/Ediz/Desktop/goat_test/out_temp_hour.csv' #output file name for the hourly temps
#which will be a column on the survey data csv

out_file_swe <- '/Users/Ediz/Desktop/goat_test/out_swe.csv' #output file name for swe in quadrants

#####################
###SET ENVIRONMENT###
#####################

#set working directory
setwd(folder)

#set years
years <- start_yr:end_yr

#load packages
library(daymetr)
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(rgeos)
library(hms)
library(lubridate)
library(chillR)
library(RCurl)
library(R.utils)

##########################
###DOWNLOAD DAYMET DATA###
##########################
  
#create daymet dir if not already there
if(dir.exists('daymet') == F) dir.create('daymet')
  
#download daymet data to daymet directory
for(tt in daymet_tiles){
  for(yr in years){
    for(param in daymet_vars) {
      daymetr::download_daymet_tiles(tiles = tt, start = yr, end = yr, param = param, path = "./daymet")
    }
  }
}

#clean up
rm(tt, yr, param)

########################################################
###CALCULATE DAYMET VARIABLES IN STUDY AREA QUADRANTS###
########################################################

#load gnp boundary, quadrants, and occupancy
gnp <- readOGR(gnp_file)
quad <- readOGR(quad_file)
occu <- raster(occu_file)

#create output datafram
dat <- data.frame(area = c(rep(c('sw', 'se', 'nw', 'ne'), 3)),
                       variable = c(rep(c('Precip.Neo', 'Temp.Neo', 
                                          'Quad.TempWint'), each = 4)))

#loop through years to calculate variables
for(i in 1:(length(years)-1)){
  
  #load temp and prcp file years
  tmax <- list.files(path = 'daymet', pattern = str_c('tmax_', years[i]), full.names = T)
  tmin <- list.files(path = 'daymet', pattern = str_c('tmin_', years[i]), full.names = T)
  tmax2 <- list.files(path = 'daymet', pattern = str_c('tmax_', years[i+1]), full.names = T)
  tmin2 <- list.files(path = 'daymet', pattern = str_c('tmin_', years[i+1]), full.names = T)
  prcp <- list.files(path = 'daymet', pattern = str_c('prcp_', years[i]), full.names = T)
  
  if(length(tmax) != 2 & length(tmin) != 2 & length(prcp) != 2 & length(tmax2) != 2 & length(tmin2) != 2){
    print(str_c('missing files for ', years[i], '-', years[i+1], '!!! CHECK DAYMET DOWNLOADS!'))
  }
  
  #load temp and prcp years, merge multiple files
  tmax <- lapply(tmax, function(x) brick(x)) %>% do.call(merge, .)
  tmin <- lapply(tmin, function(x) brick(x)) %>% do.call(merge, .)
  tmax2 <- lapply(tmax2, function(x) brick(x)) %>% do.call(merge, .)
  tmin2 <- lapply(tmin2, function(x) brick(x)) %>% do.call(merge, .)
  prcp <- lapply(prcp, function(x) brick(x)) %>% do.call(merge, .)
  
  #stack oct to apr temp separately
  tmax_wint <- stack(tmax[[274:365]], tmax2[[1:120]])
  tmin_wint <- stack(tmin[[274:365]], tmin2[[1:120]])
  
  #keep may 15 to june 15 for current year
  prcp_neo <- prcp[[135:165]]
  tmax_neo <- tmax[[135:165]]
  tmin_neo <- tmin[[135:165]]
  
  #clear other rasters
  rm(tmax, tmax2, tmin, tmin2, prcp)
  
  #calc average daily temp
  tavg_neo <- (tmax_neo + tmin_neo) / 2
  tavg_wint <- (tmax_wint + tmin_wint) / 2
  
  #change names of tavg_neo to match prcp_neo for output df calculated below
  names(tavg_neo) <- names(prcp_neo)
  
  #change quad and gnp to crs of raster
  quad_temp <- spTransform(quad, CRSobj = crs(tavg_neo))
  
  #reproject occupany layer to match daymet
  occu_temp <- occu %>% raster::aggregate(fact = 10) %>% projectRaster(to = tavg_neo)
  
  #set mask values in occupancy layer
  occu_temp[occu_temp < 0.004] <- NA
  
  #mask prcp and temp values
  prcp_mask <- mask(prcp_neo, occu_temp)
  tavg_mask <- mask(tavg_neo, occu_temp)
  tavg_mwint <- mask(tavg_wint, occu_temp)
  
  #extract daily prcp and avg temp for neo time period and avg temp for winter
  prcp_out <- raster::extract(prcp_mask, quad_temp, fun = mean, na.rm = T, df = T, weights = T)
  temp_out <- raster::extract(tavg_mask, quad_temp, fun = mean, na.rm = T, df = T, weights = T)
  wint_out <- raster::extract(tavg_mwint, quad_temp, fun = mean, na.rm = T, df = T, weights = T)
  
  #remove id columns
  prcp_out <- prcp_out[, -(colnames(prcp_out) %in% 'ID')]
  temp_out <- temp_out[, -(colnames(temp_out) %in% 'ID')]
  wint_out <- wint_out[, -(colnames(wint_out) %in% 'ID')]
  
  #calc sums of prcp and averages of tavg
  out <- rbind(data.frame(o = rowSums(prcp_out, na.rm = T)),
               data.frame(o = rowMeans(temp_out, na.rm = T)),
               data.frame(o = rowMeans(wint_out, na.rm = T))) %>% round(2)
  colnames(out) <- years[i]
  
  #bind to dat
  dat <- cbind(dat, out[, 1, drop = F])
  
  #clean up 
  rm(quad_temp, prcp_mask, tavg_mask, out, occu_temp,tavg_mwint, prcp_out, temp_out, wint_out,
     tavg_neo, tavg_wint, tmax_neo, tmin_neo, tmax_wint, tmin_wint, prcp_neo)
  removeTmpFiles(h = 0.000001)
}

rm(i)

#write daymet variables to file
write.csv(dat, out_file_daymet_vars, row.names = F)

######################################################################
###CALCULATE HOURLY TEMPERATURE CORRESPONDING TO SURVEY DATES/TIMES###
######################################################################

#load survey data
surv <- read.csv(surv_file)
surv[surv == ""] <- NA
surv$SiteName <- as.character(surv$SiteName)

#organize dates in surv
surv$Date <- surv$Date %>% mdy
surv$year <- year(surv$Date)
surv$doy <- yday(surv$Date)

#load viewshed data and gbuffer... some issues with extent
viewshed <- readOGR(viewshed_file) %>% gBuffer(byid=TRUE, width=0)

#combine shapes of separate parts of viewshed into single shapes
viewshed <- raster::aggregate(viewshed, by = "NAME_1")

#check for matching names in goat survey data and viewshed polygons
surv$SiteName %>% unique %>% sort
viewshed$NAME_1 %>% unique %>% sort

#need to re-order before changing the names on the viewshed data
#organize names here using the results of surv_names and view_names but hard code order
SiteName <- c("Apikuni Falls", "Autumn Creek", "Avalanche Lake", "Beaver Woman Lake", "Boulder Pass", 
                "Brown Pass", "Coal Creek", "Cobalt Lake", "Cosley Lake", "Cut Bank Creek", 
                "Dry Fork Creek", "Elizabeth Lake", "Fifty Mountain", "Firebrand Pass", "Grace Lake", 
                "Gunsight Pass", "Harrison Lake", "Haystack Butte", "Hidden Lake", "Iceberg Lake", 
                "Janet Lake", "Numa Lookout", "Ole Creek", "Otokomi Lake", "Park Creek", 
                "Pitamakin Pass", "Poia Lake", "Preston Park", "Red Eagle Lake", "Scenic Point", 
                "Siyeh Pass Loop", "Sperry Chalet", "Swiftcurrent Lookout", "Triple Divide Pass", 
                "Trout Lake", "Upper Kintla", "Upper Nyack")

names(SiteName) <-  c("ApikuniFalls", "AutumnCrk", "Avalanche", "BeaverWmn", "BouldrPs", "BrownPs", 
                        "CoalCrk", "CobaltLk", "CosleyLk", "CutBnk", "Dryfork", "ElzbthLk", "FiftyMtn", 
                        "FrbrndPs", "GraceLk", "GnsghtPs", "Harrison", "HystckBt", "HiddenLk", "IcebrgLk", 
                        "JanetLk", "NumaLO", "OleCreek", "Otokomi", "ParkCrk", "Pitimakn", 
                        "PoiaLake", "PrestnPk", "RedEagle", "Scenic", "SiyehPs", "Sperry", "SwftcrntLO", 
                        "TrplDivide", "TroutLk", "UpperKntl", "NyackCrk") 

#update viewshed shapefile with site name column
viewshed$SiteName <- viewshed$NAME_1 %>% recode(., !!!SiteName)

#extract daymet max and min temp within each area and for each relavent date
#add columns for max and min temp to survey data
surv$tmax <- NA
surv$tmin <- NA

#show unique survey years. might want to ensure you download enough daymet data to cover all years
#and set start_yr and end_yr accordingly
unique(surv$year)

for(i in 1:length(unique(surv$year))){

  #load temp file years
  tmax <- list.files(path = 'daymet', pattern = str_c('tmax_', unique(surv$year)[i]), full.names = T)
  tmin <- list.files(path = 'daymet', pattern = str_c('tmin_', unique(surv$year)[i]), full.names = T)
  
  if(length(tmax) != 2 & length(tmin) != 2){
    print(str_c('missing files for ', yr))
  }
  
  #load temp years
  tmax <- lapply(tmax, function(x) brick(x)) %>% do.call(merge, .)
  tmin <- lapply(tmin, function(x) brick(x)) %>% do.call(merge, .)
  
  #transform viewshed crs
  viewshed_temp <- spTransform(viewshed, crs(tmax))
  
  for(site in unique(surv$SiteName)){
    surv_hold <- surv[surv$SiteName == site & surv$year == unique(surv$year)[i],]
    if(NROW(surv_hold) != 0){
      #loop through rows to extract data
      for(j in 1:NROW(surv_hold)){
        surv_hold$tmax[j] <- raster::extract(tmax[[surv_hold$doy[j]]], viewshed[viewshed$SiteName == site,], 
                                             fun = mean, na.rm = T, df = F, weights = T) %>% as.numeric
        surv_hold$tmin[j] <- raster::extract(tmin[[surv_hold$doy[j]]], viewshed[viewshed$SiteName == site,], 
                                             fun = mean, na.rm = T, df = F, weights = T) %>% as.numeric
        
        surv[surv$SiteName == site & surv$year == unique(surv$year)[i],] <- surv_hold
      }}
    rm(surv_hold)
  }
  rm(tmax, tmin, viewshed_temp, j)
}

rm(i)  

#calculate mean temp
surv$tmean <- (surv$tmax + surv$tmin)/2  

#derive optimal hour of survey
#first coerce all time columns into time objects then numeric, so we can find max and min values
surv$N.Start <- surv$N.Start %>% hm %>% hms %>% as.numeric
surv$N.End <- surv$N.End %>% hm %>% hms %>% as.numeric
surv$S.Start <- surv$S.Start %>% hm %>% hms %>% as.numeric
surv$S.End <- surv$S.End %>% hm %>% hms %>% as.numeric
surv$DetectTime <- surv$DetectTime %>% hm %>% hms %>% as.numeric
surv$ManualTime <- surv$ManualTime %>% hm %>% hms %>% as.numeric  

#create hour column
surv$hour <- NA

#next denote how to extract time
for(i in 1:NROW(surv)){
  if(sum(is.na(surv$N.Start[i]), is.na(surv$N.End[i]), is.na(surv$S.Start[i]), is.na(surv$S.End[i])) <= 2){
    surv$hour[i] <- (max(surv$N.Start[i], surv$N.End[i], surv$S.Start[i], surv$S.End[i], na.rm = T) + 
                       min(surv$N.Start[i], surv$N.End[i], surv$S.Start[i], surv$S.End[i], na.rm = T))/2
  } else if(sum(is.na(surv$N.Start[i]), is.na(surv$N.End[i]), is.na(surv$S.Start[i]), is.na(surv$S.End[i])) == 3){
    surv$hour[i] <- max(surv$N.Start[i], surv$N.End[i], surv$S.Start[i], surv$S.End[i], na.rm = T)
  } else if(sum(is.na(surv$N.Start[i]), is.na(surv$N.End[i]), is.na(surv$S.Start[i]), is.na(surv$S.End[i])) == 4){
    surv$hour[i] <- max(surv$DetectTime[i], surv$ManualTime[i], na.rm = T)
  }
}

#check for no missing hour values
sum(is.na(surv$hour))

#find centroid of each polygon and convert to lat lon
cent <- gCentroid(viewshed, byid = T) %>% spTransform(., "+init=epsg:4326")
cent$site_name <- viewshed$SiteName

#add latitude to surv data
surv <- surv %>% mutate(lat = apply(surv, 1, function(x) cent@coords[cent$site_name == x['SiteName'],2]))

#calculate hourly temperature and extract value for each survey
#we need days before and after date of interest so makes sense to just 
#calculate on entire time series for each site location
#and use surv data to extract what we need

#first create time series using list so only have to load each year once
site_temp <- list()

for(i in 1:length(unique(surv$SiteName))){
  site_temp[[i]] <- data.frame(SiteName = unique(surv$SiteName)[i], 
                               Year = rep(min(unique(surv$year)):max(unique(surv$year)), each = 365),
                               JDay = rep(1:365, times = 12), Tmax = NA, Tmin = NA)
  names(site_temp)[i] <- unique(surv$SiteName)[i] %>% as.character
}

#load tmax and tmin data into lists
for(yr in min(unique(surv$year)):max(unique(surv$year))){
  
  #load temp file years
  tmax <- list.files(path = 'daymet', pattern = str_c('tmax_', yr), full.names = T)
  tmin <- list.files(path = 'daymet', pattern = str_c('tmin_', yr), full.names = T)
  
  if(length(tmax) != 2 & length(tmin) != 2){
    print(str_c('missing files for ', yr))
  }
  
  #load temp years
  tmax <- lapply(tmax, function(x) brick(x)) %>% do.call(merge, .)
  tmin <- lapply(tmin, function(x) brick(x)) %>% do.call(merge, .)
  
  #transform viewshed crs
  viewshed <- spTransform(viewshed, crs(tmax))
  
  for(i in 1:length(unique(surv$SiteName))){
    site_temp[[as.character(unique(surv$SiteName)[i])]]$Tmax[site_temp[[i]]$Year == yr] <- 
      raster::extract(tmax, viewshed[viewshed$SiteName == as.character(unique(surv$SiteName)[i]),], 
                      fun = mean, na.rm = T, df = F, weights = T) %>% as.numeric
    site_temp[[as.character(unique(surv$SiteName)[i])]]$Tmin[site_temp[[i]]$Year == yr] <- 
      raster::extract(tmin, viewshed[viewshed$SiteName == as.character(unique(surv$SiteName)[i]),], 
                      fun = mean, na.rm = T, df = F, weights = T) %>% as.numeric
  }
}

#run using lapply
site_temp <- lapply(site_temp, function(x) make_hourly_temps(surv$lat[surv$SiteName == as.character(x$SiteName[1])][1],
                                                             x))

#extract hourly values we want
#convert hour column back to hms
surv$hour <- surv$hour %>% as_hms

#change hour column name to time
surv <- surv %>% rename(time = hour)

#generate rounded hour column
surv$hour <- surv$time %>% hour

#create hourly temp column
surv$Temp.Hour <- NA

#extract hourly values -- easy and fast for loop
for(i in 1:NROW(surv)){
  surv$Temp.Hour[i] <- site_temp[[surv[i,'SiteName']]] %>% filter(Year == surv[i,'year'] & JDay == surv[i,'doy']) %>%
    select(str_c('Hour_', surv[i,'hour']) %>% str_squish) %>% as.numeric
}
rm(i)

#coerce all time columns back into time objects
surv$N.Start <- surv$N.Start %>% as_hms
surv$N.End <- surv$N.End %>% as_hms
surv$S.Start <- surv$S.Start %>% as_hms
surv$S.End <- surv$S.End %>% as_hms
surv$DetectTime <- surv$DetectTime %>% as_hms
surv$ManualTime <- surv$ManualTime %>% as_hms

#clean up df
surv$tmax <- surv$tmax %>% round(digits = 2)
surv$tmin <- surv$tmin %>% round(digits = 2)
surv$tmean <- surv$tmean %>% round(digits = 2)
surv$Temp.Hour <- surv$Temp.Hour %>% round(digits = 2)

#write Temp.Hour data to file
write.csv(surv, out_file_temp_hour, row.names = F)

##########################
###DOWNLOAD SNODAS DATA###
##########################

#create snodas dir if not already there
if(dir.exists('snodas') == F) dir.create('snodas')

#create raw dir if not already there
if(dir.exists('snodas/raw') == F) dir.create('snodas/raw')

#set download directory
raw_dir <- str_c(getwd(), '/snodas/raw')

for(yr in years){
  url <- str_c("ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/", yr, "/")
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
      download.file(paste(url2, filename, sep = ""), paste(raw_dir, "/", filename,
                                                           sep = ""))
    }
    rm(url2, filenames)
  }
  rm(url, mon)
}

#load gnp shape
gnp <- readOGR(gnp_file)

#create temp dir and output files dir
if(dir.exists('snodas/temp') == F) dir.create('snodas/temp')
if(dir.exists('snodas/extracted') == F) dir.create('snodas/extracted')

#create dir paths
temp_dir <- str_c(getwd(), '/snodas/temp')
extr_dir <- str_c(getwd(), '/snodas/extracted')

#get files in raw_dir
files <- list.files(raw_dir, '.tar', full.names = T)

#find leap year files
leap <- str_which(files, '0229')

#remove leap year files
file.remove(files[leap])

#get files in dir
files <- list.files(raw_dir, '.tar', full.names = T)

#check to see if any files are missing
#extract dates from file names and parse as Date object
check <- str_extract(files, '[:digit:]{8}')
check <- ymd(check)

#create sequence along timeseries
dates <- seq(ymd(str_c(start_yr, "0101")), ymd(str_c(end_yr, "1231")), by = 1)

#print which are missing
#seems these files just don't exist so when I stack the data i'll
#just leave whole bands as NA
dates[!(dates %in% check)]

#we only want oct to april so only extract those months
months <- c(1, 2, 3, 4, 10, 11, 12)
files <- files[month(check) %in% months]

#set to run for swe only
vari <- matrix(c('swe', 'us_ssmv11034'), ncol = 2)

#run loop to extract files
for(file in files){
  #untar main file
  untar(file, exdir = temp_dir)
  
  #loop through, unzip, crop, and write .tif for 3 variables
  for(i in 1:nrow(vari)){
    
    #find file for variable
    f <- list.files(temp_dir, pattern = glob2rx(str_c(vari[i,2], '*.dat.gz')), full.names = T)
    f_out <- str_c(vari[i,1], '_', str_extract(file, '[:digit:]{8}'), '.dat')
    gunzip(f[1], destname = str_c(temp_dir, '/', f_out))
    
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
    writeLines(hdr, str_c(temp_dir, '/', str_replace(f_out, '.dat', '.hdr')))
    
    #load in dat file
    sno <- raster(str_c(temp_dir, '/', f_out))
    
    #set crs
    crs(sno) <- CRS("+init=epsg:4326")
    
    #reproj gnp outline
    gnp_sno <- spTransform(gnp, crs(sno))
    
    #crop sno and write to extracted folder
    sno <- crop(sno, gnp_sno, filename = str_c(extr_dir, '/', str_replace(f_out, '.dat', '.tif')),
                format = 'GTiff', datatype = 'INT2S', NAflag = -9999)
    
    #clean up
    rm(f, f_out, hdr, gnp_sno, sno)
  }
  #clear temp folder
  invisible(do.call(file.remove, list(list.files(temp_dir, full.names = TRUE))))
}


###########################################
###CALCULATE SWE IN STUDY AREA QUADRANTS###
###########################################

#create output dataframe
sno_out <- data.frame(area = c('sw', 'se', 'nw', 'ne'),
                      variable = rep('Quad.SWE', each = 4))

for(i in 1:(length(years)-1)){
  
  #load list of all snodas files
  files <- list.files('snodas/extracted', pattern = '.tif', full.names = T)
  
  #convert to date object
  dates <- str_extract(files, '[:digit:]{8}') %>% ymd
  
  #create date object with dates we want
  target <- seq(ymd(str_c(years[i], '1001')), ymd(str_c(years[i+1], '0430')), by = 1)
  
  #mask files based on dates we want
  files <- files[dates %in% target]
  rm(dates, target)
  
  #load raster stack
  swe <- stack(str_subset(files, 'swe'))
  
  #transform shapes
  quad_sno <- spTransform(quad, crs(swe))
  
  #reproject occupany layer to match snodas
  occu_sno <- occu %>% raster::aggregate(fact = 10) %>% projectRaster(to = swe)
  
  #set mask values in occupancy layer
  occu_sno[occu_sno < 0.004] <- NA
  
  #mask values
  swe_mask <- mask(swe, occu_sno)
  
  #apply scale correction
  swe_mask <- swe_mask/1000
  
  #extract average values for entire winter
  swe_yr <- raster::extract(swe_mask, quad_sno, fun = mean, na.rm = T, df = T, weights = T) %>%
    .[, -1]
  
  sno_yr <- data.frame(o = rowMeans(swe_yr, na.rm = T)) %>% round(4)
  colnames(sno_yr) <- str_c(years[i])
  
  #cbind to df
  sno_out <- cbind(sno_out, sno_yr[, 1, drop = F])
  
  #clean up
  rm(gnp_sno, quad_sno, occu_sno, sno_yr, swe_mask,
     swe, files, swe_yr)
}
rm(i)

#write to disk and clean up
write.csv(sno_out, out_file_swe)


