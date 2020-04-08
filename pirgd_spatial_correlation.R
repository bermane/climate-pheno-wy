#load packages
library(rgdal)
library(sp)
library(raster)
library(tidyverse)
library(gstat)
library(usdm)
library(foreach)
library(doParallel)
library(tictoc)
library(viridis)


#set wd
setwd('/Volumes/SSD/climate_effects')

#set parallel env
getDoParWorkers()
registerDoParallel(cores = 8)
getDoParWorkers()

#load image
load(file='/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/semivario.RData')

###LOAD DATA###
#load dlc data
pirgd <- stack('dlc/maxIRGdate_wy_laea_2000_2019.tif')

#test running different variogram sizes on an individual year
r <- pirgd[[1]]
r_var_20 <- usdm::Variogram(r, cutoff = 20000, size = 100)
r_var_50 <- usdm::Variogram(r, cutoff = 50000, size = 100)
r_var_100 <- usdm::Variogram(r, cutoff = 100000, size = 100)
r_var_200 <- usdm::Variogram(r, cutoff = 200000, size = 20)

#save.image(file='/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/semivario.RData')

#using cutoff = 100000 (100 km) and size = 100, run for every year
#loop through each year to create list and joint df
years = 2000:2019

#generate annual results using dopar
tic()
var_100_ls <- foreach(i = 1:length(years)) %dopar% {
  usdm::Variogram(pirgd[[i]], cutoff = 100000, size = 100)
}
toc()

#save.image(file='/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/semivario.RData')

#extract dataframes from vario objects and add year variable
var_100_df <- lapply(seq_along(var_100_ls), function(x) data.frame(distance = var_100_ls[[x]]@variogram$distance, 
                                                                   gamma = var_100_ls[[x]]@variogram$gamma, year = 1999 + x))

#convert to single dataframe for plotting dataframe
var_df <- do.call("rbind", var_100_df)

#plot results
ggplot(var_df, aes(x = distance, y = gamma, color = year)) +
  geom_point(size = 0.5) + theme_bw() + 
  scale_x_continuous(name = "Distance (km)", labels = c("0" = "0", "25000" = "25", 
                                                   "50000" = "50", "75000" = "75", 
                                                   "100000" = "100")) +
  scale_y_continuous(name = "Gamma") + ggtitle("Semivariogram of PIRGd up to 100 km with 100 Samples/Year") +
  scale_color_viridis()


#playing with other spatial autocorrelation indices in the raster package

#Moran(r) #this is the global index of autocorrelation
#x1  <-  MoranLocal(r) #local measure of autocorr as a raster object that can be plotted
#plot(x1) #this will plot the autocorrelation raster results
#Geary(r)
