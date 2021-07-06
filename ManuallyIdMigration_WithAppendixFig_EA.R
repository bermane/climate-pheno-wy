require(lubridate)

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


nsd<-function(x, y){
  ns <- length(x)
  msd <- NULL
  for(i in 1:ns){
    msd <- c(msd, ((x[i]-x[1])^2 + (y[i]-y[1])^2))
  }
  return(msd)
}

setwd("C:/Users/eaikens/Box Sync/Wyoming Range_EA/Data/GPS collars/Tidy")

full.df<-read.csv("AllDeer2015_08_26.csv")

full.df$timestamp<-ymd_hms(full.df$POSIX) 

n<-names(full.df)

#just need to rename column names, so that everything works in the code below
n[9:10]<-c("utm.x", "utm.y")
names(full.df)<-n
#full.df<-full.df[order(full.df$AnimalID, full.df$timestamp), ]
ids<-unique(full.df$AnimalID)

full.df$NSD<-NA

#calculate Net Squared Displacement 
full.df.temp<-full.df[0, ]
for(d in 1:length(ids)){
  di<-subset(full.df, AnimalID==ids[d])
  di$NSD<-nsd(di$utm.x, di$utm.y)
  full.df.temp<-rbind(full.df.temp, di)
}
full.df<-full.df.temp

#visualize each deer, and used NSD to id start/end of spring and fall migration
par(mfrow=c(1,2))

years<-sort(unique(full.df$Year))

migration.dates<-data.frame(AnimalID=ids[1], year=years[1],
                startSpring=full.df$timestamp[1],
                endSpring=full.df$timestamp[1],
                startFall=full.df$timestamp[1],
                endFall=full.df$timestamp[1])
m.dates<-migration.dates[0, ]
j<-61
for(j in 1:length(ids)){
  #subset each deer
  full.i<-subset(full.df, AnimalID==ids[j])

  #id start and end of migration for each year
  for(y in years){
    repeat{
      
    #subest for year
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
      plot(full.i$utm.x, full.i$utm.y, asp=1, col="grey50", pch=16)
      points(d.y$utm.x, d.y$utm.y, asp=1, col=trans, pch=16)
      plot(full.i$timestamp, full.i$NSD, type="l", col="grey50", main=ids[j])
      points(d.y$timestamp, d.y$NSD, col=trans, pch=16)
      lines(d.y$timestamp, d.y$NSD)

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
      
      numericT <- as.numeric(d.y$timestamp)
      
      days<-ddays(10)
      #have end of spring
      if(is.na(endSpring)!=T){
        es.idx<-which(abs(numericT-endSpring$x)==min(abs(numericT-endSpring$x)))[1]
        endS.d<-numericT[es.idx]
        endSpring<-d.y$timestamp[es.idx]
        
        es.plus<-endSpring+days
        es.mins<-endSpring-days
        
        spring <- d.y[numericT <= endS.d,]
        
      }
      #have start of spring
      if(is.na(startSpring)!=T){
        ss.idx<-which(abs(numericT-startSpring$x)==min(abs(numericT-startSpring$x)))[1]
        startS.d<-numericT[ss.idx]
        startSpring<-d.y$timestamp[ss.idx]
        
        
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
        endFall<-d.y$timestamp[ef.idx]
        
        fall <- d.y[numericT <= endF.d,]
        
        ef.plus<-endFall+days
        ef.mins<-endFall-days
        
      }
      #have start of fall
      if(is.na(startFall)!=T){
        sf.idx<-which(abs(numericT-startFall$x)==min(abs(numericT-startFall$x)))[1]
        startF.d<-numericT[sf.idx]
        startFall<-d.y$timestamp[sf.idx]
        
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
        plot(d.y$timestamp, d.y$NSD, type="l", main="NDS")
        
        points(spring$timestamp, spring$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(spring)))) 
        abline(v=c(ss.mins, ss.plus), col="red", lty=2)
        
        plot(d.y$timestamp, d.y$NSD, type="l", xlim=c(ss.mins, ss.plus), main="Start Spring Migration")
        points(spring$timestamp, spring$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(spring)))) 
        startSpring <- locator(1)
        ss.idx<-which(abs(numericT-startSpring$x)==min(abs(numericT-startSpring$x)))[1]
        
        startS.d<-numericT[ss.idx]
        startSpring<-d.y$timestamp[ss.idx]
        
      }
      
      if(is.na(endSpring)!=T){
        
        plot(d.y$timestamp, d.y$NSD, type="l", main="NDS")
        
        points(spring$timestamp, spring$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(spring)))) 
        abline(v=c(es.mins, es.plus), col="red", lty=2)
        
        plot(d.y$timestamp, d.y$NSD, type="l", xlim=c(es.mins, es.plus), main="End Spring Migration")
        points(spring$timestamp, spring$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(spring)))) 
        endSpring <- locator(1)
        es.idx<-which(abs(numericT-endSpring$x)==min(abs(numericT-endSpring$x)))[1]
        
        endS.d<-numericT[es.idx]
        endSpring<-d.y$timestamp[es.idx]
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
        plot(d.y$timestamp, d.y$NSD, type="l", main="NDS")
        
        points(fall$timestamp, fall$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(fall)))) 
        abline(v=c(sf.mins, sf.plus), col="red", lty=2)
        
        plot(d.y$timestamp, d.y$NSD, type="l", xlim=c(sf.mins, sf.plus), main="Start Fall Migration")
        points(fall$timestamp, fall$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(fall)))) 
        startFall <- locator(1)
        sf.idx<-which(abs(numericT-startFall$x)==min(abs(numericT-startFall$x)))[1]
        
        startF.d<-numericT[sf.idx]
        startFall<-d.y$timestamp[sf.idx]
        
      }
      
      if(is.na(endFall)!=T){
        plot(d.y$timestamp, d.y$NSD, type="l", main="NDS")
        
        points(fall$timestamp, fall$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(fall)))) 
        abline(v=c(ef.mins, ef.plus), col="red", lty=2)
        
        plot(d.y$timestamp, d.y$NSD, type="l", xlim=c(ef.mins, ef.plus), main="End Fall Migration")
        points(fall$timestamp, fall$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(fall)))) 
        endFall <- locator(1)
        ef.idx<-which(abs(numericT-endFall$x)==min(abs(numericT-endFall$x)))[1]
        
        endF.d<-numericT[ef.idx]
        endFall<-d.y$timestamp[ef.idx]
        fall <- d.y[numericT <= endF.d,]
      }
      
      
       if(is.na(endFall)!=T){
         fall <-fall[numericT >= startF.d,]
       }
       else{
         fall <-d.y[numericT >= startF.d,]
       }
      
      
      #situation when you have at least start dates for both fall and spring
      #par(mfrow=c(1,2))
      if(is.na(startSpring)!=T && is.na(startFall)!=T){
        plot(d.y$utm.x, d.y$utm.y, asp=1, 
             col=makeTransparent("black"), pch=16)
        points(spring$utm.x, spring$utm.y, 
               asp=1, col=makeTransparent(rainbow(nrow(spring))), 
               pch=16)
        
        points(fall$utm.x, fall$utm.y, 
               asp=1, col=makeTransparent(topo.colors(nrow(fall))), 
               pch=16)
        
        plot(d.y$timestamp, d.y$NSD, type="l")
        points(spring$timestamp, spring$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(spring))))
        
        points(fall$timestamp, fall$NSD, pch=16, 
               col=makeTransparent(topo.colors(nrow(fall))))
      }
      
      #situation with only have spring date(s)
      else if(is.na(startSpring)!=T){
        plot(d.y$utm.x, d.y$utm.y, asp=1, 
             col=makeTransparent("black"), pch=16)
        
        points(spring$utm.x, spring$utm.y, 
               asp=1, col=makeTransparent(rainbow(nrow(spring))), 
               pch=16)
        
        
        plot(d.y$timestamp, d.y$NSD, type="l")
        
        points(spring$timestamp, spring$NSD, pch=16, 
               col=makeTransparent(rainbow(nrow(spring))))
        
      }
      
      #situation with only have fall date(s)
      else if(is.na(startFall)!=T){
        plot(d.y$utm.x, d.y$utm.y, asp=1, 
             col=makeTransparent("black"), pch=16)
        
        points(fall$utm.x, fall$utm.y, 
               asp=1, col=makeTransparent(topo.colors(nrow(fall))), 
               pch=16)
        
        plot(d.y$timestamp, d.y$NSD, type="l")
        
        points(fall$timestamp, fall$NSD, pch=16, 
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


#example for appendix 
d.y<-full.df[full.df$Year==2013 &full.df$AnimalID==33, ]

startS.d<-1368612900
endS.d<-1371341700
startF.d<-1381564980
endF.d<-1385345700

require("RColorBrewer")
col.spring<-colorRampPalette(c("#a6bddb", "#67a9cf", "#3690c0", "#02818a", 
                               "#016c59"))(nrow(spring))

col.fall<-colorRampPalette(c("#fec44f", "#fe9929", "#ec7014", "#cc4c02", 
                             "#993404"))(nrow(fall))

numericT <- as.numeric(d.y$timestamp)
fall <- d.y[numericT >= startF.d & numericT<=endF.d,]
spring <- d.y[numericT >= startS.d & numericT<=endS.d,]


setwd("C:/Users/eaikens/Box Sync/Wyoming Range_EA/Figures/For Manuscript/Phenology Tracking Variability/high resolution final")
tiff(filename = "NSDtoIDMigration.tif",  width = 970*6, 
     height = 917*6, units = "px", res=600,
     compression ="none", bg = "white")
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  par(mar=c(5, 5, 0.5, 0.5))
  
  plot(d.y$timestamp, d.y$NSD, type="l", lwd=3, 
       xlab="Time", ylab="NSD", bty="L", cex.lab=2, 
       cex.axis=1.5, xlim=c(min(d.y$timestamp), max(d.y$timestamp)+months(1)), 
       ylim=c(0, max(d.y$NSD)+500000000))
  
  points(spring$timestamp, spring$NSD, pch=16, 
         col=col.spring, cex=1.5)
  
  points(fall$timestamp, fall$NSD, pch=16, 
         col=col.fall, cex=1.5)
  arrows(1366915401, 355036586, 1368218565, 83735044, lwd=2, length = 0.15)
  text(1365515008, 590540820, "Start spring migration", cex=1.5)
  
  arrows(1386519514, 355036586, 1385443987, 83735044, lwd=2, length = 0.15)
    text(1388016129, 590540820, "End fall migration", cex=1.5)
  
  arrows(1368671839, 5509765883, 1370258299, 5509765883, lwd=2, length = 0.15)
  text(1365515008, 5509765883, "End spring migration", cex=1.5)
  
  arrows(1383686549, 5464548959, 1382100089, 5464548959, lwd=2, length = 0.15)
  text(1386600089, 5464548959, "Start fall migration", cex=1.5)
  
  text(1363015008, 6191765883, "(a)", cex=2, font=2)
  
  plot(d.y$utm.x, d.y$utm.y, asp=1, xlab="UTM Easting", ylab="UTM Northing",
       col=makeTransparent("black"), pch=16, cex.lab=2, cex=1.5,
       cex.axis=1.5, bty="L", ylim=c(min(d.y$utm.y), 4752500))
  
  points(spring$utm.x, spring$utm.y, 
         asp=1, col=col.spring, pch=16, cex=1.5)
  text(512500, 4752500, "(b)", cex=2, font=2)
  
  plot(d.y$utm.x, d.y$utm.y, asp=1, xlab="UTM Easting", ylab="UTM Northing",
       col=makeTransparent("black"), pch=16, cex.lab=2, cex=1.5,
       cex.axis=1.5, bty="L", ylim=c(min(d.y$utm.y), 4752500))
  
  points(fall$utm.x, fall$utm.y, cex=1.5,
         asp=1, col=col.fall, pch=16)
  text(512500, 4752500, "(c)", cex=2, font=2)

dev.off()




#for right now, just deal with cases when you have start and end of both spring and fall migrations
dates<-m.dates[complete.cases(m.dates[, c(4:6)]), ]
ids<-unique(as.character(dates$AnimalID))
winter <- spring <- fall <- summer <-full.df[0, ]
for(i in 1:length(ids)){
  years<-unique(dates[dates$AnimalID==ids[i],]$year) 
  for(y in 1:length(years)){
    d.i<-subset(full.df, AnimalID==ids[i])
    
    ss<-dates[dates$AnimalID==ids[i]& dates$year==years[y], ]$startSpring
    es<-dates[dates$AnimalID==ids[i]& dates$year==years[y], ]$endSpring
    sf<-dates[dates$AnimalID==ids[i]& dates$year==years[y], ]$startFall
    ef<-dates[dates$AnimalID==ids[i]& dates$year==years[y], ]$endFall
    
    springM<-d.i[d.i$timestamp>=ss & d.i$timestamp<=es, ]
    summerR<-d.i[d.i$timestamp>es & d.i$timestamp<sf, ]
    fallM<-d.i[d.i$timestamp>=sf & d.i$timestamp<=ef, ]
    #if we have data from previous year
    if(length(years)>=2 & years[y]>=years[2]){
        sw<-dates[dates$AnimalID==ids[i]& dates$year==years[y-1],]$endFall
        winterR<-d.i[d.i$timestamp>sw & d.i$timestamp<ss, ]
    }
    else{
    #start of data collection - start of spring migration of that year (either 2013 or 2014 if first year of data collection)
        winterR<-d.i[d.i$timestamp<ss, ]
    }
    #jpeg(paste("deer_", 
    #           ifelse(nchar(ids[i])==1, paste("00", ids[i], sep=""), 
    #                  ifelse(nchar(ids[i])==2, paste("0", ids[i], sep=""), ids[i])),"_", y, ".jpg", sep=""), 
    #     960, 960, quality=100)
    par(mfrow=c(2, 2))
    #check everything by plotting
    
    #winter
    plot(d.i$utm.x, d.i$utm.y, asp=1, col=rgb(0,0,0,.25), main=paste("Winter deer", ids[i], years[y]), pch=16)
    points(winterR$utm.x, winterR$utm.y, col=rgb(1, 0, 0, .25), pch=16)
    #spring migration
    plot(d.i$utm.x, d.i$utm.y, asp=1, col=rgb(0,0,0,.25), main=paste("Spring deer", ids[i], years[y]), pch=16)
    points(springM$utm.x, springM$utm.y, col=rgb(0, 1, 0, .25), pch=16)
    #Summer range
    plot(d.i$utm.x, d.i$utm.y, asp=1, col=rgb(0,0,0,.25), main=paste("Summer deer", ids[i], years[y]), pch=16)
    points(summerR$utm.x, summerR$utm.y, col=rgb(0, 0, 1, .25), pch=16)
    #Fall Migration
    plot(d.i$utm.x, d.i$utm.y, asp=1, col=rgb(0,0,0,.25), main=paste("Fall deer", ids[i], years[y]), pch=16)
    points(fallM$utm.x, fallM$utm.y, col=rgb(1, 0, 1, .25), pch=16)
    #dev.off()
    #save data 
    winter<-rbind(winter, winterR)
    spring<-rbind(spring, springM)
    summer<-rbind(summer, summerR)
    fall<-rbind(fall, fallM)
  }
}


write.csv(winter, "2015_08_26WinterRange.csv")
write.csv(spring, "2015_09_15SpringMigration.csv")
write.csv(summer, "2015_08_26SummerRange.csv")
write.csv(fall, "2015_08_26FallMigration.csv")

setwd("C:/Users/eaikens/Box Sync/Wyoming Range_EA/Data/GPS collars")
write.csv(m.dates, "MigrationDatesAll.csv")

