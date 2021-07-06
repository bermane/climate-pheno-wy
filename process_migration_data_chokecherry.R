#this code processes raw deer migration data by herd file
#first sampling individuals and than
#using Visvalingam’s algorithm for line simplification and resampling of specific
#individuals at 4 km points
#the first part is based off Ellen Aiken's R script 
#'ManuallyIdMigration_WithAppendixFig

#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(rmapshaper)
library(foreign)
library(adehabitatLT)
library(sf)
library(lubridate)
library(swaRm)
library(berryFunctions)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load WY shapefile
wy <- readOGR('reference/wyoming.shp')

###################################################
###SELECT INDIVIDUAL SEGMENTS USING ELLEN'S CODE###
###################################################

#some functions used later on 
makeTransparent <- function(..., alpha=0.5){
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha <- floor(255*alpha)  
  newColor <- col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent <- function(col, alpha){
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor <- apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
}


nsd <-function(x, y){
  ns <- length(x)
  msd <- NULL
  for(i in 1:ns){
    msd <- c(msd, ((x[i]-x[1])^2 + (y[i]-y[1])^2))
  }
  return(msd)
}

#load chokecherry data
load('mig_data/Organized Data/Mule_WY_Chokecherry_cleaned.RData')

#change loaded data name to gps.data if necessary
gps.data <- data.sub
rm(data.sub)

#load as df
gps_df <- gps.data@data

#add coordinates
gps_df <- cbind(gps_df, data.frame(x = gps.data@coords[,1], y = gps.data@coords[,2]))

#change date to POSIXct object if necessary
str(gps_df$Timestamp)
#gps_df$Timestamp <- as.POSIXct(gps_df$Timestamp, tz = 'GMT')

#load unique animal ids
ids <- unique(gps_df$AnimalID)

#allocate column for NSD
gps_df$NSD <- NA

#calculate Net Squared Displacement 
gps_df_temp <- gps_df[0, ]
for(d in 1:length(ids)){
  di<-subset(gps_df, AnimalID==ids[d])
  di$NSD <- swaRm::nsd(di$x, di$y, geo = TRUE)
  gps_df_temp<-rbind(gps_df_temp, di)
}

#put into main df name used in ellen's code
full.df <- gps_df_temp
rm(gps_df_temp)

#create year column
full.df$Year <- year(full.df$Timestamp)

#visualize each deer, and used NSD to id start/end of spring and fall migration
par(mfrow=c(1,2))

years <- sort(unique(full.df$Year))

migration.dates <- data.frame(AnimalID=ids[1], year=years[1],
                            startSpring=full.df$Timestamp[1],
                            endSpring=full.df$Timestamp[1],
                            startFall=full.df$Timestamp[1],
                            endFall=full.df$Timestamp[1])

m.dates <- migration.dates[0, ]

for(j in 1:length(ids)){
  #subset each deer
  full.i<-subset(full.df, AnimalID==ids[j])
  
  #id start and end of migration for each year
  for(y in years){
    repeat{
      
      #subset for year
      d.y<-full.i[full.i$Year==y, ]
      
      #in the case we don't have data for that year
      if(nrow(d.y)==0){
        startSpring <- NA
        endSpring <- NA
        startFall<- NA
        endFall<- NA
        break
      }
      else{
        
        trans<-makeTransparent(rainbow(nrow(d.y)))
        plot(full.i$x, full.i$y, asp=1, col="grey50", pch=16)
        points(d.y$x, d.y$y, asp=1, col=trans, pch=16)
        plot(full.i$Timestamp, full.i$NSD, type="l", col="grey50", main=ids[j])
        points(d.y$Timestamp, d.y$NSD, col=trans, pch=16)
        lines(d.y$Timestamp, d.y$NSD)
        
        pts<-as.numeric(readline("how many events do you want to id?"))
        if(pts==4){
          startSpring <- locator(1)
          endSpring <- locator(1)
          startFall<- locator(1)
          endFall<- locator(1)
        }
        else if (pts==0){
          startSpring <- NA
          endSpring <- NA
          startFall<- NA
          endFall<- NA
          break
        }
        else{
          
          ss<-readline("do you want to id the start of spring migration (Y or N)?")
          if(ss=="Y"|ss=="y"){
            startSpring <- locator(1)
          }
          else{
            startSpring <- NA
          }
          
          es<-readline("do you want to id the end of spring migration (Y or N)?")
          if(es=="Y"|es=="y"){
            endSpring <- locator(1)
          }
          else{
            endSpring <- NA
          }
          
          sf<-readline("do you want to id the start of fall migration (Y or N)?")
          if(sf=="Y"|sf=="y"){
            startFall <- locator(1)
          }
          else{
            startFall <- NA
          }
          
          ef<-readline("do you want to id the end of fall migration (Y or N)?")
          if(ef=="Y"|ef=="y"){
            endFall <- locator(1)
          }
          else{
            endFall <- NA
          }
          
        }
        
        #estop<-readline("does this deer use extended stopovers (Y/N)?")
        
        numericT <- as.numeric(d.y$Timestamp)
        
        days<-ddays(10)
        #have end of spring
        if(is.na(endSpring)!=T){
          es.idx<-which(abs(numericT-endSpring$x)==min(abs(numericT-endSpring$x)))[1]
          endS.d<-numericT[es.idx]
          endSpring<-d.y$Timestamp[es.idx]
          
          es.plus<-endSpring+days
          es.mins<-endSpring-days
          
          spring <- d.y[numericT <= endS.d,]
          
        }
        #have start of spring
        if(is.na(startSpring)!=T){
          ss.idx<-which(abs(numericT-startSpring$x)==min(abs(numericT-startSpring$x)))[1]
          startS.d<-numericT[ss.idx]
          startSpring<-d.y$Timestamp[ss.idx]
          
          
          ss.plus<-startSpring+days
          ss.mins<-startSpring-days
          #have start of spring and end of spring
          if(is.na(endSpring)!=T){
            #spring already defined above, adding in the start cut off
            spring <-spring[numericT >= startS.d,]
          }
          #have just 
          else{
            #spring just defined by the start cut off
            spring <-d.y[numericT >= startS.d,]
          }
          
        }
        
        
        #have end of fall
        if(is.na(endFall)!=T){
          ef.idx<-which(abs(numericT-endFall$x)==min(abs(numericT-endFall$x)))[1]
          endF.d<-numericT[ef.idx]
          endFall<-d.y$Timestamp[ef.idx]
          
          fall <- d.y[numericT <= endF.d,]
          
          ef.plus<-endFall+days
          ef.mins<-endFall-days
          
        }
        #have start of fall
        if(is.na(startFall)!=T){
          sf.idx<-which(abs(numericT-startFall$x)==min(abs(numericT-startFall$x)))[1]
          startF.d<-numericT[sf.idx]
          startFall<-d.y$Timestamp[sf.idx]
          
          sf.plus<-startFall+days
          sf.mins<-startFall-days
          
          #have start and end of fall
          if(is.na(endFall)!=T){
            #fall already defined above, adding start of fall cut-off
            fall <- fall[numericT >= startF.d,]
          }
          else if(exists("startF.d")){
            #fall only defined by start of fall cut off
            fall <- d.y[numericT >= startF.d,]
          }
          
          
        }
        
        #checking your choices
        
        #zoom to each selection
        
        #only plot if spring data exists
        if(is.na(startSpring)!=T){
          plot(d.y$Timestamp, d.y$NSD, type="l", main="NDS")
          
          points(spring$Timestamp, spring$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(spring)))) 
          abline(v=c(ss.mins, ss.plus), col="red", lty=2)
          
          plot(d.y$Timestamp, d.y$NSD, type="l", xlim=c(ss.mins, ss.plus), main="Start Spring Migration")
          points(spring$Timestamp, spring$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(spring)))) 
          startSpring <- locator(1)
          ss.idx<-which(abs(numericT-startSpring$x)==min(abs(numericT-startSpring$x)))[1]
          
          startS.d<-numericT[ss.idx]
          startSpring<-d.y$Timestamp[ss.idx]
          
        }
        
        if(is.na(endSpring)!=T){
          
          plot(d.y$Timestamp, d.y$NSD, type="l", main="NDS")
          
          points(spring$Timestamp, spring$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(spring)))) 
          abline(v=c(es.mins, es.plus), col="red", lty=2)
          
          plot(d.y$Timestamp, d.y$NSD, type="l", xlim=c(es.mins, es.plus), main="End Spring Migration")
          points(spring$Timestamp, spring$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(spring)))) 
          endSpring <- locator(1)
          es.idx<-which(abs(numericT-endSpring$x)==min(abs(numericT-endSpring$x)))[1]
          
          endS.d<-numericT[es.idx]
          endSpring<-d.y$Timestamp[es.idx]
          spring <- d.y[numericT <= endS.d,]
        }
        
        #redefine spring with more exact dates 
        if(is.na(endSpring)!=T){
          spring <-spring[numericT >= startS.d,]
        }
        else if(exists("startS.d")){
          spring <-d.y[numericT >= startS.d,]
        }
        
        
        
        if(is.na(startFall)!=T){
          plot(d.y$Timestamp, d.y$NSD, type="l", main="NDS")
          
          points(fall$Timestamp, fall$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(fall)))) 
          abline(v=c(sf.mins, sf.plus), col="red", lty=2)
          
          plot(d.y$Timestamp, d.y$NSD, type="l", xlim=c(sf.mins, sf.plus), main="Start Fall Migration")
          points(fall$Timestamp, fall$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(fall)))) 
          startFall <- locator(1)
          sf.idx<-which(abs(numericT-startFall$x)==min(abs(numericT-startFall$x)))[1]
          
          startF.d<-numericT[sf.idx]
          startFall<-d.y$Timestamp[sf.idx]
          
        }
        
        if(is.na(endFall)!=T){
          plot(d.y$Timestamp, d.y$NSD, type="l", main="NDS")
          
          points(fall$Timestamp, fall$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(fall)))) 
          abline(v=c(ef.mins, ef.plus), col="red", lty=2)
          
          plot(d.y$Timestamp, d.y$NSD, type="l", xlim=c(ef.mins, ef.plus), main="End Fall Migration")
          points(fall$Timestamp, fall$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(fall)))) 
          endFall <- locator(1)
          ef.idx<-which(abs(numericT-endFall$x)==min(abs(numericT-endFall$x)))[1]
          
          endF.d<-numericT[ef.idx]
          endFall<-d.y$Timestamp[ef.idx]
          fall <- d.y[numericT <= endF.d,]
        }
        
        
        if(is.na(endFall)!=T){
          fall <-fall[numericT >= startF.d,]
        }
        else{
          #fall <-d.y[numericT >= startF.d,]
        }
        
        
        #situation when you have at least start dates for both fall and spring
        #par(mfrow=c(1,2))
        if(is.na(startSpring)!=T && is.na(startFall)!=T){
          plot(d.y$x, d.y$y, asp=1, 
               col=makeTransparent("black"), pch=16)
          points(spring$x, spring$y, 
                 asp=1, col=makeTransparent(rainbow(nrow(spring))), 
                 pch=16)
          
          points(fall$x, fall$y, 
                 asp=1, col=makeTransparent(topo.colors(nrow(fall))), 
                 pch=16)
          
          plot(d.y$Timestamp, d.y$NSD, type="l")
          points(spring$Timestamp, spring$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(spring))))
          
          points(fall$Timestamp, fall$NSD, pch=16, 
                 col=makeTransparent(topo.colors(nrow(fall))))
        }
        
        #situation with only have spring date(s)
        else if(is.na(startSpring)!=T){
          plot(d.y$x, d.y$y, asp=1, 
               col=makeTransparent("black"), pch=16)
          
          points(spring$x, spring$y, 
                 asp=1, col=makeTransparent(rainbow(nrow(spring))), 
                 pch=16)
          
          
          plot(d.y$Timestamp, d.y$NSD, type="l")
          
          points(spring$Timestamp, spring$NSD, pch=16, 
                 col=makeTransparent(rainbow(nrow(spring))))
          
        }
        
        #situation with only have fall date(s)
        else if(is.na(startFall)!=T){
          plot(d.y$x, d.y$y, asp=1, 
               col=makeTransparent("black"), pch=16)
          
          points(fall$x, fall$y, 
                 asp=1, col=makeTransparent(topo.colors(nrow(fall))), 
                 pch=16)
          
          plot(d.y$Timestamp, d.y$NSD, type="l")
          
          points(fall$Timestamp, fall$NSD, pch=16, 
                 col=makeTransparent(topo.colors(nrow(fall))))
        }
        
      }
      
      repe<-readline("Repeat selection for this animal (Y or N)?")
      if(repe=="N"|repe=="n"){
        break
      }
    }
    
    
    startSpring<-as.POSIXct(startSpring, origin="1970-01-01")
    endSpring<-as.POSIXct(endSpring, origin="1970-01-01")
    startFall<-as.POSIXct(startFall, origin="1970-01-01")
    endFall<-as.POSIXct(endFall, origin="1970-01-01")
    m.row<-data.frame(AnimalID=as.character(ids[j]), year=y, 
                      startSpring=startSpring, endSpring=endSpring,  
                      startFall=startFall, endFall=endFall)
    
    m.dates<-rbind(m.dates, m.row)
  }
  
}

m.dates  

#save image of initial m.dates output to make sure don't lose it!
#save.image(file='/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/segment selection/mig_segment_selection_chokecherry.RData')

#remove rows with missing spring start and end
m.dates <- m.dates[is.na(m.dates$startSpring) == F,]
m.dates <- m.dates[is.na(m.dates$endSpring) == F,]

#remove fall columns don't need them
m.dates <- m.dates[, !(colnames(m.dates) %in% c('startFall', 'endFall'))]

#allocate column for migration distance
m.dates$migDist <- as.numeric(NA)

#loop through segments to identify distance between start and end date
for(i in 1:NROW(m.dates)){
  
  #create points object of start and end
  if(is.error(start_end <- st_sfc(st_point(as.numeric(full.df[full.df$AnimalID == m.dates$AnimalID[i] & 
                                                              full.df$Timestamp == m.dates$startSpring[i], c('x', 'y')])),
                                  st_point(as.numeric(full.df[full.df$AnimalID == m.dates$AnimalID[i] & 
                                                              full.df$Timestamp == m.dates$endSpring[i], c('x', 'y')])),
                                  crs = as.character(crs(gps.data)))) == F){
    start_end <- st_sfc(st_point(as.numeric(full.df[full.df$AnimalID == m.dates$AnimalID[i] & 
                                                      full.df$Timestamp == m.dates$startSpring[i], c('x', 'y')])),
                        st_point(as.numeric(full.df[full.df$AnimalID == m.dates$AnimalID[i] & 
                                                      full.df$Timestamp == m.dates$endSpring[i], c('x', 'y')])),
                        crs = as.character(crs(gps.data)))
    
    #calculate distance between first and last point and output as km
    m.dates$migDist[i] <- (st_distance(start_end)/1000)[1,2] %>% as.numeric
  }
  
}

#remove NA values
m.dates <- na.omit(m.dates)

#check df for distances and use 15 km cut off if possible
m.dates <- m.dates[m.dates$migDist > 15,]

#randomly sample animal years (no individual overlaps)
#first group by animal id and sample one segment from each id
m.smpl <- m.dates %>% group_by(AnimalID) %>% sample_n(1)

#next sample 15 individuals
m.smpl <- m.smpl %>% ungroup %>% sample_n(15)

#create spdf of start and end points to use for duration metrics
#lopo through each row to get start and end
for(i in 1:NROW(m.smpl)){
  if(i == 1){
    
    #load relevant rows from full.df
    row_start <- full.df[full.df$AnimalID == m.smpl$AnimalID[i] &
                     full.df$Timestamp == m.smpl$startSpring[i],]
    row_end <- full.df[full.df$AnimalID == m.smpl$AnimalID[i] &
                         full.df$Timestamp == m.smpl$endSpring[i],]
      
    mig_start_end <- rbind(SpatialPointsDataFrame(coords = matrix(c(row_start$x, row_start$y), ncol = 2),
                                     proj4string = crs(gps.data),
                                     data = data.frame(AnimalID = row_start$AnimalID,
                                                       Timestamp = row_start$Timestamp,
                                                       Population = row_start$Population,
                                                       AY_ID = row_start$AY_ID,
                                                       GlobalID = row_start$GlobalID)),
                    SpatialPointsDataFrame(coords = matrix(c(row_end$x, row_end$y), ncol = 2),
                                           proj4string = crs(gps.data),
                                           data = data.frame(AnimalID = row_end$AnimalID,
                                                             Timestamp = row_end$Timestamp,
                                                             Population = row_end$Population,
                                                             AY_ID = row_end$AY_ID,
                                                             GlobalID = row_end$GlobalID)))
  } else{
    
    #load relevant rows from full.df
    row_start <- full.df[full.df$AnimalID == m.smpl$AnimalID[i] &
                           full.df$Timestamp == m.smpl$startSpring[i],]
    row_end <- full.df[full.df$AnimalID == m.smpl$AnimalID[i] &
                         full.df$Timestamp == m.smpl$endSpring[i],]
    
    mig_start_end <- rbind(mig_start_end, 
                    SpatialPointsDataFrame(coords = matrix(c(row_start$x, row_start$y), ncol = 2),
                                           proj4string = crs(gps.data),
                                           data = data.frame(AnimalID = row_start$AnimalID,
                                                             Timestamp = row_start$Timestamp,
                                                             Population = row_start$Population,
                                                             AY_ID = row_start$AY_ID,
                                                             GlobalID = row_start$GlobalID)),
                    SpatialPointsDataFrame(coords = matrix(c(row_end$x, row_end$y), ncol = 2),
                                           proj4string = crs(gps.data),
                                           data = data.frame(AnimalID = row_end$AnimalID,
                                                             Timestamp = row_end$Timestamp,
                                                             Population = row_end$Population,
                                                             AY_ID = row_end$AY_ID,
                                                             GlobalID = row_end$GlobalID)))
    
  }
}

#output file of start and end dates
#save(mig_start_end, file = ('mig_data/processed/seg_start_end_chokecherry.RData'))


#########################
###LINE SIMPLIFICATION###
#########################

#run the line simplification as a loop that lets you
#choose the % the keep, view the results, and iterate through 
#until happy with the output. Then save all final sample points 
#in a spdf together

#start count
count <- 1

#change back to single plot
par(mfrow=c(1,1))

#loop through sample rows
for(i in 1:NROW(m.smpl)){
  
  #load segment data
  seg <- full.df[full.df$AnimalID == m.smpl$AnimalID[i] & 
                   full.df$Timestamp >= m.smpl$startSpring[i] & full.df$Timestamp <= m.smpl$endSpring[i], ]
  
  #convert to spdf
  seg <- SpatialPointsDataFrame(coords = matrix(c(seg$x, seg$y), ncol = 2),
                                proj4string = crs(gps.data),
                                data = data.frame(AnimalID = seg$AnimalID,
                                                  Timestamp = seg$Timestamp,
                                                  Population = seg$Population,
                                                  AY_ID = seg$AY_ID,
                                                  GlobalID = seg$GlobalID))
  
  #convert to spLines
  lines <- spLines(seg)
  
  #repeat
  repeat{
    
    #plot
    plot(lines)
    
    #choose simplification percentage
    perc <- readline('What percentage of points to retain? ') %>% as.numeric
    
    #simplify migration route
    simp <- ms_simplify(lines, keep = perc)
    
    #check plot
    plot(simp, col = 'red', add = T)
    
    #check if line is sufficient
    check <- readline('Re-run line simplification (Y/N)? ')
    if(check == 'N' | check == 'n') break
    
  }
  
  #repeat
  repeat{
    check <- readline('Use normal sampling interval (Y/N)? ')
    if(check == 'Y' | check == 'y'){
      int <- round(SpatialLinesLengths(simp, longlat = T)/5) + 1
    } else int <- readline('Specify numeric sampling interval: ')
    #sample points along line at 5 km intervals. include start and finish
    sample <- spsample(simp, int, type = 'regular')
    
    #plot
    plot(sample, add = T, cex =2, pch = 1)
    
    check <- readline('Re-run sampling (Y/N)? ')
    if(check == 'N' | check == 'n') break
  }
  
  #if first iteration
  #save samples in spdf
  if(count == 1){
    samples <- SpatialPointsDataFrame(coords = sample@coords,
                                      proj4string = crs(gps.data),
                                      data = data.frame(AnimalID = seg$AnimalID[1],
                                                        SampleID = 1:length(sample),
                                                        Population = as.character(seg$Population[1]),
                                                        AY_ID = as.character(seg$AY_ID[1]),
                                                        GlobalID = as.character(seg$GlobalID[1]),
                                                        SamplePerc = perc))
    
    #up counter
    count <- count + 1
    
  } else{ #rbind if not first iteration
    samples <- rbind(samples, 
                     SpatialPointsDataFrame(coords = sample@coords,
                                            proj4string = crs(gps.data),
                                            data = data.frame(AnimalID = seg$AnimalID[1],
                                                              SampleID = 1:length(sample),
                                                              Population = as.character(seg$Population[1]),
                                                              AY_ID = as.character(seg$AY_ID[1]),
                                                              GlobalID = as.character(seg$GlobalID[1]),
                                                              SamplePerc = perc)))
    
  }
}

#output file of migration samples
#save(samples, file = ('mig_data/processed/mig_samples_chokecherry.RData'))

