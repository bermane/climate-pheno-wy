#this code transforms the calculated migration samples back to lines
#since I failed to save the lines themselves!
#might need to also include start and end points

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(sf)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load mig samples files and start/end files
mig <- list.files('/Volumes/SSD/climate_effects/mig_data/processed', 'mig_samples', full.names = T)
seg <- list.files('/Volumes/SSD/climate_effects/mig_data/processed', 'seg_start_end', full.names = T)
names <- list.files('/Volumes/SSD/climate_effects/mig_data/processed', 'mig_samples')

for(i in 1:length(mig)){
  
  #load mig samples and start/end
  load(mig[i])
  load(seg[i])
  
  #set count
  count <- 1
  
  #loop through unique animal IDs
  for(j in unique(samples$GlobalID)){
    
    #subset points
    samp <- subset(samples, GlobalID == j)
    hold <- subset(mig_start_end, GlobalID == j)
    
    #change hold to meet columns of samp
    hold <- hold[,!(colnames(hold@data) %in% 'Timestamp')]
    hold$SampleID <- c(0, length(samp) + 1)
    hold$SamplePerc <- samp$SamplePerc[1]
    
    #rbind together
    samp <- rbind(samp, hold)
    
    #sort by SampleID
    samp <- samp[order(samp$SampleID),]
    
    #convert to sf points
    samp <- samp %>% st_as_sf %>% st_geometry %>% st_cast(., 'POINT')
    
    #build linestrings
    linestrings <- lapply(X = 1:(length(samp) - 1), FUN = function(x) {
      
      pair <- st_combine(c(samp[x], samp[x + 1]))
      line <- st_cast(pair, "LINESTRING")
      return(line)
    })
    
    #One MULTILINESTRING object with all the LINESTRINGS
    multilinetring <- st_multilinestring(do.call("rbind", linestrings))
    
    if(count == 1){
      lines <- as(multilinetring, 'Spatial')
      lines$GlobalID <- j
    }else{
      l_hold <- as(multilinetring, 'Spatial')
      l_hold$GlobalID <- j
      lines <- rbind(lines, l_hold)
    }

    count <- count + 1
  }
  
  #create export filename string
  exp <- str_c('/Volumes/SSD/climate_effects/mig_data/processed/shp/',
               str_replace(names[i], 'samples', 'lines')) %>%
    str_replace(., 'RData', 'shp')
  
  #create sldf
  crs(lines) <- crs(samples)
  row.names(lines) <- as.character(1:length(lines))
  lines <- SpatialLinesDataFrame(lines, data.frame(GlobalID = lines$GlobalID))
  
  #export to shp
  writeOGR(lines, 
           dsn = exp,
           layer = "lines",
           driver="ESRI Shapefile",
           overwrite_layer = T)
  
}
