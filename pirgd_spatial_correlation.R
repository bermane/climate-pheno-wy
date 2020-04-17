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
library(Rmisc)


#set wd
setwd('/Volumes/SSD/climate_effects')

#set parallel env
#getDoParWorkers()
#registerDoParallel(cores = 8)
#getDoParWorkers()

#load image
load(file='/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/semivario.RData')

#set seed for random sampling
set.seed(3)

###LOAD DATA###
#load dlc data
#pirgd <- stack('dlc/maxIRGdate_wy_laea_2000_2019.tif')

#subset from 2001-2018
#pirgd <- pirgd[[2:19]]

#test running different variogram sizes on an individual year
#r <- pirgd[[1]]
#r_var_20 <- usdm::Variogram(r, cutoff = 20000, size = 100)
#r_var_50 <- usdm::Variogram(r, cutoff = 50000, size = 100)
#r_var_100 <- usdm::Variogram(r, cutoff = 100000, size = 100)
#r_var_200 <- usdm::Variogram(r, cutoff = 200000, size = 20)

#save.image(file='/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/semivario.RData')

#using cutoff = 100000 (100 km) and size = 100, run for every year
#loop through each year to create list and joint df
#years = 2000:2019

#generate annual results using dopar
#tic()
#var_100_ls <- foreach(i = 1:length(years)) %dopar% {
#  usdm::Variogram(pirgd[[i]], cutoff = 100000, size = 100)
#}
#toc()

#save.image(file='/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/semivario.RData')

#extract dataframes from vario objects and add year variable
#var_100_df <- lapply(seq_along(var_100_ls), function(x) data.frame(distance = var_100_ls[[x]]@variogram$distance, 
#                                                                   gamma = var_100_ls[[x]]@variogram$gamma, year = 1999 + x))

#remove 2000 and 2019
#var_100_df <- var_100_df[2:19]

#convert to single dataframe for plotting dataframe
#var_df <- do.call("rbind", var_100_df)

#plot results
ggplot(var_df, aes(x = distance, y = gamma, color = year)) +
  geom_point(size = 0.5) + theme_classic() + 
  scale_x_continuous(name = "Distance (km)", labels = c("0" = "0", "25000" = "25", 
                                                   "50000" = "50", "75000" = "75", 
                                                   "100000" = "100")) +
  scale_y_continuous(name = "Gamma") + ggtitle("Semivariogram of PIRGd (2001-2018) up to 100 km with 100 Unique Samples/Year") +
  scale_color_viridis() + labs(color = "Year") + geom_vline(xintercept = 16000, linetype = "dotted", size = 1) +
  annotate(geom = "text", x=21000, y=650, label = "16 km", fontface = 2)
ggsave("output/semivariogram.png")

#what would it look like to sample landscape at 16 km?
#16 km is every 250 m * 4 * 16 = every 64th cell

#find the number of points per row and column
#NROW(r)/(4*16)
#NCOL(r)/(4*16)

#take a sample at regular intervals
#rsp <- sampleRegular(r, size = 31*38, ext = extent(r), cells=FALSE, xy=FALSE, asRaster=FALSE, 
#              sp=TRUE, useGDAL=FALSE)

#load WY shapefile and transform to proj of dlc
#wy <- readOGR('./reference/wyoming.shp')
#wy <- spTransform(wy, crs(pirgd))

#plot over PIRGd layer and WY shape
plot(r, main = "PIRGd Year 2000 (DOY) Overlaid with WY Boundary \n and 1147 Sample Points (~16 km grid)")
plot(wy, add = T)
plot(rsp, add = T)

#export sampling points
#writeOGR(obj = rsp, dsn = './reference', layer = 'sampling_points', driver = "ESRI Shapefile")

#extract pirgd sample from all layers into one vector
pirgd_ex <- raster::extract(pirgd, rsp)
pirgd_ex <- as.vector(pirgd_ex)

#make sample of PIRGd from whole study area
pirgd_ma <- as.array(pirgd)
pirgd_ma <- as.vector(pirgd_ma)

#plot histogram of sampled points
p1 <- ggplot(as.data.frame(pirgd_ex), aes(x = pirgd_ex)) +
  geom_histogram(color = "darkblue", fill = "lightblue", breaks = seq(0,250, by = 5)) + theme_classic() +
  labs(title = "PIRGd of Sampling Grid 2001-2018 (20646 samples)", x = "PIRGd (DOY)", y = "Count")
#ggsave("output/pirgd_data_sample_hist.png")

#plot histogram of total points
#take sample of 100000 points
pirgd_ma <- sample(pirgd_ma, 100000)

p2 <- ggplot(as.data.frame(pirgd_ma), aes(x = pirgd_ma)) +
  geom_histogram(color = "darkblue", fill = "lightblue", breaks = seq(0,250, by = 5)) + theme_classic() +
  labs(title = "PIRGd of Landscape 2001-2018 (100000 samples)", x = "PIRGd (DOY)", y = "Count")
#ggsave("output/pirgd_landscape_hist.png")

#load elevation and generate histograms
dem <- raster('/Volumes/SSD/climate_effects/dem/srtm/dem_wy_laea.tif')

#extract dem samples
dem_ex <- raster::extract(dem, rsp)

#make sample of DEM from whole study area
dem_ma <- as.vector(dem)

#plot histogram of sampled points
p3 <- ggplot(as.data.frame(dem_ex), aes(x = dem_ex)) +
  geom_histogram(color = "black", fill = "brown", breaks = seq(0,4000, by = 100)) + theme_classic() +
  labs(title = "Elevation of Sampling Grid (1147 samples)", x = "Elevation (m)", y = "Count")
#ggsave("output/elevation_data_sample_hist.png")

#plot histogram of total points
#take sample of 100000 points
dem_ma <- sample(dem_ma, 100000)

p4 <- ggplot(as.data.frame(dem_ma), aes(x = dem_ma)) +
  geom_histogram(color = "black", fill = "brown", breaks = seq(0,4000, by = 100)) + theme_classic() +
  labs(title = "Elevation of Landscape (100000 samples)", x = "Elevation (m)", y = "Count")
#ggsave("output/elevation_landscape_hist.png")

#multiplot
png("output/sampling_grid_landscape_hist.png", width = 2700, height = 1800, res = 250)
multiplot(p1, p2, p3, p4, cols = 2)
dev.off()

#load landcover and generate barplots
lc <- raster('landcover/five_class_landcover_wy_laea_2016.tif')

#extract lc samples into one vector
lc_ex <- raster::extract(lc, rsp) 

#organize to plot
lc_ex <- as.data.frame(lc_ex[is.na(lc_ex) == F])
colnames(lc_ex) <- "lc"
lc_ex$lc <- as.factor(lc_ex$lc)

#make sample of lc from whole study area
lc_ma <- as.vector(lc)

#plot histogram of sampled points
p5 <- ggplot(lc_ex, aes(x = lc, fill = lc)) +
  geom_bar() + theme_classic() + 
  scale_x_discrete(labels = c("1" = "Shrub", "2" = "Grassland", "3" = "Deciduous Forest", 
                              "4" = "Evergreen Forest", "5" = "Wetlands")) +
  scale_fill_manual(values = c('orangered3', 'goldenrod4', 'maroon3', 'green4', 'dodgerblue4')) +
  labs(title = "Landcover of Sampling Grid (1147 samples)", x = "", y = "Count") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), title = element_text(size = 14))
#ggsave("output/landcover_data_sample_bar.png")

#plot histogram of total points
#take sample of 100000 points
lc_ma <- sample(lc_ma[is.na(lc_ma) == F], 100000)

#organize to plot
lc_ma <- as.data.frame(lc_ma[is.na(lc_ma) == F])
colnames(lc_ma) <- "lc"
lc_ma$lc <- as.factor(lc_ma$lc)

p6 <- ggplot(lc_ma, aes(x = lc, fill = lc)) +
  geom_bar() + theme_classic() + 
  scale_x_discrete(labels = c("1" = "Shrub", "2" = "Grassland", "3" = "Deciduous Forest", 
                              "4" = "Evergreen Forest", "5" = "Wetlands")) +
  scale_fill_manual(values = c('orangered3', 'goldenrod4', 'maroon3', 'green4', 'dodgerblue4')) +
  labs(title = "Landcover of Landscape (100000 samples)", x = "", y = "Count") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), title = element_text(size = 14))
#ggsave("output/landcover_landscape_bar.png")

#multiplot
png("output/lc_sampling_grid_landscape_bar.png", width = 2000, height = 1800, res = 250)
multiplot(p5, p6)
dev.off()

#playing with other spatial autocorrelation indices in the raster package

#Moran(r) #this is the global index of autocorrelation
#x1  <-  MoranLocal(r) #local measure of autocorr as a raster object that can be plotted
#plot(x1) #this will plot the autocorrelation raster results
#Geary(r)
