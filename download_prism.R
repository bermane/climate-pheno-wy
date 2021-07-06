#code to download and extract climate variables from PRISM
#across the WY study area

#load packages
library(tidyverse)

#set download dir
#if different variables needed run in different windows
#to increase downloading speed
setwd('/Volumes/SSD/climate_effects/prism/prcp')

#set daily time series
time <- seq(from = as.Date("2000-01-01"), to = as.Date("2019-12-31"), by = 1)

#remove leap year dates
time <- time[!grepl('02-29', time)]

#convert to format we want
time <- time %>% as.character %>% str_replace_all(pattern = '-', replacement = '') 

for(t in time){
  #get URL
  url <- str_c('http://services.nacse.org/prism/data/public/4km/ppt/', t)
  
  #download
  download.file(url, destfile = str_c(t, '.zip'), method = 'wget')
}


