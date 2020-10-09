#this analysis builds a model for the drivers of pirgd in WY
#it is run on the data extracted with the extract_data_table function 
#plus additional variables added in individual codes

#load packages
library(rgdal)
library(sp)
library(raster)
library(dynatopmodel)
library(tidyverse)
library(data.table)
library(lme4)
library(tictoc)
library(scales)
library(nlme)
library(berryFunctions)
library(effects)
library(ggplot2)
library(car)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load source functions
source("reference/HighstatLibV10.R")

#load workspace
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_09_17.RData')

##################################
###REMOVE MISSING/UNWANTED DATA###
##################################
#put the variables with different elevation levels into single columns
dat$rain_elev <- rowSums(data.frame(dat$rain_low_mar_apr, dat$rain_med_mmar_mmay, dat$rain_hi_apr_may), na.rm = T)
dat$gdd_elev <- rowSums(data.frame(dat$gdd_low_apr, dat$gdd_med_mapr_mmay, dat$gdd_hi_may), na.rm = T)
dat$vp_min_elev <- rowSums(data.frame(dat$vp_min_low_mar_apr, dat$vp_min_med_mmar_mmay, dat$vp_min_hi_apr_may), na.rm = T)
dat$vp_avg_elev <- rowSums(data.frame(dat$vp_avg_low_mar_apr, dat$vp_avg_med_mmar_mmay, dat$vp_avg_hi_apr_may), na.rm = T)

#change evergreen to 3 instead of 4
dat$lc[dat$lc == 4] <- 3

#load into new dataframe
dat2 <- dat

#remove columns with missing data due to elevation levels
dat2 <- select(dat2, -c(rain_low_mar_apr, rain_med_mmar_mmay, rain_hi_apr_may,
                        gdd_low_apr, gdd_med_mapr_mmay, gdd_hi_may,
                        vp_min_low_mar_apr, vp_min_med_mmar_mmay, vp_min_hi_apr_may,
                        vp_avg_low_mar_apr, vp_avg_med_mmar_mmay, vp_avg_hi_apr_may))

#remove homer annual layers since we are moving forward without them
#remove columns with missing data due to elevation levels
dat2 <- select(dat2, -c(herb_homer_ann, sage_homer_ann, shrub_homer_ann))

#remove any remaining missing data points for PIRGd or snowmelt
#2 points for pirgd
dat2 <- dat2[is.na(dat2$pirgd) == 0,]
#dat2 <- dat2[is.na(dat2$ss) == 0,]
#200 points for snowmelt
dat2 <- dat2[is.na(dat2$snowmelt) == 0,]

#any other missing points?
#may need to remove the variables that include snowmelt and pirgd...
#16% of samples have PIRGd before snowmelt...
sum(is.na(dat2))

#################################
###TEST VARIABLES FOR OUTLIERS###
#################################

plot(x = dat2$dem, y = dat2$pirgd)

#####################################################
###TEST VARIABLES FOR NORMALITY AND EQUAL VARIANCE###
#####################################################

#id col and lc as factors
dat2$id <- as.factor(dat2$id)
dat2$lc <- as.factor(dat2$lc)

#loop through all covariates in table
#check to make sure correct starting column -- start with dem!!
colnames(dat2)[7]
ncol(dat2)
for(col_id in 7:53){
  
  #set varible name
  vari <- colnames(dat2)[col_id]
  
  #fit a linear mixed model with landcover, random effect, and variable
  m <- eval(substitute(lme(pirgd ~ lc + variable, random = ~ 1 | id, data = dat2),
                         list(variable = as.name(vari))))

  #plot residuals
  e2 <- resid(m, type = 'normalized')
  f2 <- fitted(m)
  
  #set up utput plot
  jpeg(str_c('output/norm_equal_var_test_plots/' ,vari, '.jpeg'), width = 900, height = 600)
  
  #plot graphs to test for normality and equal variance
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  plot(x = f2, y = e2, xlab = 'Fitted values', ylab = 'Residuals', main = str_c(vari, ' equal variance?'))
  hist(e2, xlab = 'Residuals', main = str_c(vari, ' normality?'))
  plot(dat2[,col_id], e2, xlab = vari, ylab = 'Residuals', main = str_c(vari, ' independence?'))
  
  #dev off
  dev.off()
  
  #clean up
  rm(m, e2, f2, op, vari, col_id)
  
}

###############################################
###TEST VARIABLES FOR QUADRATIC RELATIONSHIP###
###############################################

#allocate output dataframe
quad_df <- data.frame(var = character(), quad_coef = numeric(), quad_95 = numeric(),
                      quad_aic = numeric(), lin_aic = numeric())

#loop through all covariates in table
#check to make sure correct starting column -- start with dem!!
colnames(dat2)[7]
ncol(dat2)
for(col_id in 7:53){
  
  #set varible name
  vari <- colnames(dat2)[col_id]
  
  #fit a linear mixed model with landcover, random effect, and variable
  m <- eval(substitute(lme(pirgd ~ lc + variable, random = ~ 1 | id, data = dat2),
                       list(variable = as.name(vari))))
  
  #get summary
  s <- summary(m)
    
  #fit with quadratic term  
  m2 <- eval(substitute(lme(pirgd ~ lc + poly(variable, 2, raw = TRUE), random = ~ 1 | id, data = dat2),
                       list(variable = as.name(vari))))
  
  #get summary
  s2 <- summary(m2)
  
  quad_df <- rbind(quad_df, data.frame(var = vari, quad_coef = s2$tTable[str_c('poly(', vari, ', 2, raw = TRUE)2'),"Value"],
                                       quad_95 = (s2$tTable[str_c('poly(', vari, ', 2, raw = TRUE)2'),"Std.Error"])*1.96,
                                       quad_aic = s2$AIC, lin_aic = s$AIC))
  
}

#set low and high confidence intervals to test overlap with zero
quad_df$quad_low <- quad_df$quad_coef - quad_df$quad_95
quad_df$quad_hi <- quad_df$quad_coef + quad_df$quad_95

#subset dataframe to variables where quad does not overlap zero
quad_df2 <- quad_df[!(quad_df$quad_low < 0 & quad_df$quad_hi > 0),]

#subset dataframe to variables where quad AIC is lower than lin AIC
quad_df2 <- quad_df2[quad_df2$quad_aic < quad_df2$lin_aic,]

#move forward with quadratic terms from these variables!
quad_var <- quad_df2$var

#clean up
rm(quad_df, quad_df2, m, m2, s, s2, vari, col_id)

#####################
###RUN TEST MODELS###
#####################

#create new df to rescale covariates
dat3 <- dat2

#which columns to rescale? Want everything after dem
colnames(dat3)[7]
ncol(dat3)
#7 to 50

#rescale variables between 0 and 1
for(i in 7:53){
  dat3[,i] <- scale(dat3[,i])
}

#test fit a linear model with only elevation
m <- lm(pirgd ~ lc + dem, dat3, method = 'ML')

#test fit a linear mixed model with only elevation
m2 <- lme(pirgd ~ lc + dem, random = ~ 1 | id, data = dat3, method = 'ML')

#check AIC of the two
AIC(m)
AIC(m2)

###########################
###LOOK INTO CORRELATION###
###########################

#check into correlation between all variables
#compute correlation matrix
#make sure to check including all columns!!!
cor_m <- cor(dat3[,7:53])
cor_m <- round(cor_m, 2)

#create output table of overall correlations that are significant
#allocate df
cor_df <- data.frame(var1 = character(), var2 = character(), cor = numeric())

#loop through
for(i in 1:nrow(cor_m)){
  for(j in 1:ncol(cor_m)){
    if(abs(cor_m[i,j]) >= .7 & rownames(cor_m)[i] != colnames(cor_m)[j]){
      cor_df <- rbind(cor_df, data.frame(var1 = rownames(cor_m)[i],
                                         var2 = colnames(cor_m)[j],
                                         cor = cor_m[i,j]))
    }
  }
}

#write to csv file
write.csv(cor_df, file = 'output/correlation results - 2020-09-17.csv')

#check correlation between bare ground and the sum of all other rap variables.
#use raw values instead of rescaled
bg <- dat2$bare_ground_rap
sum_rap <- dat2$shrub_rap + dat2$ann_forb_rap + dat2$perenn_forb_rap +
  dat2$tree

#calc correlation
cor(bg, sum_rap)

#very high correlation, remove bare ground from analysis below
rm(bg, sum_rap)

#check into correlation between pirgd and ss
pairs(dat3[, c('pirgd', 'ss')], upper.panel = panel.smooth,
      lower.panel = panel.cor)

#correlation across different temperature variables and their relationship to pirgd
pairs(dat3[, c("pirgd", "gdd_jan_pirgd", "gdd_jan_earlyavg_pirgd","gdd_elev",
               "chill_oct_apr", "chill_jan_apr", "chill_jan_may", "solar_jan_pirgd",
               "solar_jan_apr", "dem")], 
      upper.panel = panel.smooth, lower.panel = panel.cor, cex.labels = 1.5)

#correlation across different precip variables and their relationship to pirgd
pairs(dat3[, c('pirgd', 'snow_oct_apr', 'rain_oct_pirgd', 'rain_jan_pirgd', 'rain_elev',
               'tot_prcp_oct_apr', 'dem', 'twi', 'vp_avg_jan_apr', 'vp_avg_elev', 'snowmelt')], 
      upper.panel = panel.smooth, lower.panel = panel.cor, cex.labels = 1.5)

#correlation across different drought variables and their relationship to pirgd
pairs(dat3[, c('pirgd', 'pdsi_apr_apr_min', 'pdsi_mar_apr_min', 'pdsi_jan_pirgd_min',
               'pdsi_apr_apr_mean', 'pdsi_mar_apr_mean', 'pdsi_jan_pirgd_mean')], 
      upper.panel = panel.smooth, lower.panel = panel.cor, cex.labels = 1.5)

#correlation across different veg component variables and their relationship to pirgd
pairs(dat3[, c('pirgd', 'herb_homer', 'sage_homer', 'shrub_homer',
               'ann_forb_rap', 'bare_ground_rap', 'perenn_forb_rap',
               'shrub_rap', 'tree_rap', 'ann_perenn_forb_rap')], 
      upper.panel = panel.smooth, lower.panel = panel.cor, cex.labels = 1.5)

########################################################
###TEST VARIABLES FOR NON-ZERO and NON-SIGN SWITCHING###
########################################################

#which columns do we want to test? everything after landcover
colnames(dat3)[6]
ncol(dat3)
#7 to 47

#allocate dataframe to save coefficients
uni_df <- data.frame(var = character(), beta = numeric(), beta_sq = numeric(), se_95 = numeric(), se_95_sq = numeric(),
                     sign = character(), sign_sq = character(), aic = numeric())

#loop through variables
for(i in 7:53){
  
  #get column to test
  var <- colnames(dat3)[i]
  
  #check if we want to model quadratic, otherwise linear term only
  if(var %in% quad_var){
    
    #fit a linear mixed model with base + quad variables -- only if converges!
    
    if(is.error(eval(substitute(lme(pirgd ~ lc + poly(variable, 2, raw = TRUE), random = ~ 1 | id, 
                                    data = dat3, method = 'ML'),
                                list(variable = as.name(var))))) == F){

      m <- eval(substitute(lme(pirgd ~ lc + poly(variable, 2, raw = TRUE), random = ~ 1 | id, 
                               data = dat3, method = 'ML'),
                           list(variable = as.name(var))))
      
      #get summary of model
      s <- summary(m)
      
      #extract values we want from model summary
      beta <- s$tTable[str_c('poly(', var, ', 2, raw = TRUE)1'),"Value"]
      beta_sq <- s$tTable[str_c('poly(', var, ', 2, raw = TRUE)2'),"Value"]
      se_95 <- s$tTable[str_c('poly(', var, ', 2, raw = TRUE)1'),"Std.Error"]*1.96
      se_95_sq <- s$tTable[str_c('poly(', var, ', 2, raw = TRUE)2'),"Std.Error"]*1.96
      aic <- s$AIC
      
      #find whether signs are zero/positive/negative
      if((beta - se_95 < 0) & (beta + se_95 < 0)) sign <- 'negative'
      if((beta - se_95 > 0) & (beta + se_95 > 0)) sign <- 'positive'
      if((beta - se_95 < 0) & (beta + se_95 > 0)) sign <- 'zero'
      
      if((beta_sq - se_95_sq < 0) & (beta_sq + se_95_sq < 0)) sign_sq <- 'negative'
      if((beta_sq - se_95_sq > 0) & (beta_sq + se_95_sq > 0)) sign_sq <- 'positive'
      if((beta_sq - se_95_sq < 0) & (beta_sq + se_95_sq > 0)) sign_sq <- 'zero'
      
      uni_df <- rbind(uni_df, data.frame(var = var, beta = beta, beta_sq = beta_sq, se_95 = se_95, 
                                         se_95_sq = se_95_sq, sign = sign, sign_sq = sign_sq, aic = aic)) 
      
    } else print(str_c(var, ' did not converge')) #save names of models that did not converge
    
  } else{ #run as linear term only!
    
    #fit a linear mixed model with base + linear variables if no convergence error!
    if(is.error(eval(substitute(lme(pirgd ~ lc + variable, random = ~ 1 | id, 
                                    data = dat3, method = 'ML'),
                                list(variable = as.name(var))))) == F){
      m <- eval(substitute(lme(pirgd ~ lc + variable, random = ~ 1 | id, 
                               data = dat3, method = 'ML'),
                           list(variable = as.name(var))))
      
      #get summary of model
      s <- summary(m)
      
      #extract values we want from model summary
      beta <- s$tTable[var,"Value"]
      se_95 <- s$tTable[var,"Std.Error"]*1.96
      aic <- s$AIC
      
      #find whether signs are zero/positive/negative
      if((beta - se_95 < 0) & (beta + se_95 < 0)) sign <- 'negative'
      if((beta - se_95 > 0) & (beta + se_95 > 0)) sign <- 'positive'
      if((beta - se_95 < 0) & (beta + se_95 > 0)) sign <- 'zero'
      
      uni_df <- rbind(uni_df, data.frame(var = var, beta = beta, beta_sq = NA, se_95 = se_95, 
                                         se_95_sq = NA, sign = sign, sign_sq = 'NA', aic = aic)) 

      
    } else print(str_c(var, ' did not converge')) #save names of models that did not converge
  }
}

#clean up
rm(m, s, beta, beta_sq, se_95, aic, sign, sign_sq, var, i, se_95_sq)

#write to csv file
write.csv(uni_df, file = 'output/univariate results - 2020-09-25.csv')

###############################
###RUN MULTIVARIATE ANALYSIS###
###############################

#the code below runs model and variable combinations based on the 
#results of the univariate analysis. Results were quantitatively, and qualitatively
#assessed to pick the best combinations and variables to move forward with and test.
#first step is to try all combinations within categories: temp, precip, veg.

############################
###RUN TEMPERATURE COMBOS###
############################

#create table of all temp combos to evaluate
combos <- as.data.frame(expand.grid(gdd_elev = 0:1, chill_jan_apr = 0:1, dem = 0:1))

#remove first row with no variables and combos with dem except dem + gdd_elev
combos <- combos[-c(1, 7, 8),]

#allocate dataframe for AIC results
temp_df <- data.frame(var = character(), aic = numeric())

#loop through all combos and run model
for(i in 1:nrow(combos)){
  
  #load row with combination
  df <- combos[i,]
  
  #allocate character vector for variable names
  vars = rep("NA", ncol(df))
  
  #loop through columns to extract variable names
  for(j in 1:ncol(df)) {
    if(df[1,j] == 1) {vars[j] <- colnames(df)[j]}
  }
  
  #remove missing names
  vars <- vars[!(vars %in% "NA")]
  
  #save variable names in readable format
  vars_name <- vars
  
  #change names to quadratic if applicable
  for(j in 1:length(vars)){
    if(vars[j] %in% quad_var) vars[j] <- str_c('poly(', vars[j], ', 2, raw = TRUE)')
  }
  
  #change names to quadratic readable if applicable
  for(j in 1:length(vars_name)){
    if(vars_name[j] %in% quad_var) vars_name[j] <- str_c(vars_name[j], '^2')
  }
  
  #add variable names together to load into model
  vars <- paste(vars, collapse = " + ")
  
  #add variable names together for readability
  vars_name <- paste(vars_name, collapse = " + ")
  
  #run model
  m <- eval(parse(text = str_c('lme(pirgd ~ lc + ', vars, ', random = ~ 1 | id, 
                               data = dat3, method = "ML")')))
  
  #save AIC
  temp_df <- rbind(temp_df, data.frame(var = vars_name, aic = AIC(m)))

}

#clean up
rm(m, i, j, vars, combos, df, vars_name)

#######################
###RUN PRECIP COMBOS###
#######################

#create table of all precip combos to evaluate
combos <- as.data.frame(expand.grid(rain_elev = 0:1, pdsi_mar_apr_min = 0:1, twi = 0:1,
                                    vpd_tmax_mean_jan_apr = 0:1, snow_oct_apr = 0:1, tot_prcp_oct_apr = 0:1,
                                    snowmelt = 0:1))

#remove first row with no variables
combos <- combos[-1,]

#allocate dataframe for AIC results
precip_df <- data.frame(var = character(), aic = numeric())

#loop through all combos and run model
for(i in 1:nrow(combos)){
  
  #load row with combination
  df <- combos[i,]
  
  #allocate character vector for variable names
  vars = rep("NA", ncol(df))
  
  #loop through columns to extract variable names
  for(j in 1:ncol(df)) {
    if(df[1,j] == 1) {vars[j] <- colnames(df)[j]}
  }
  
  #remove missing names
  vars <- vars[!(vars %in% "NA")]
  
  #save variable names in readable format
  vars_name <- vars
  
  #only run if no correlated variables
  if(sum(vars %in% c("snow_oct_apr", 'tot_prcp_oct_apr', "snowmelt")) <= 1) {
  
  #change names to quadratic if applicable
  for(j in 1:length(vars)){
    if(vars[j] %in% quad_var) vars[j] <- str_c('poly(', vars[j], ', 2, raw = TRUE)')
  }
  
  #change names to quadratic readable if applicable
  for(j in 1:length(vars_name)){
    if(vars_name[j] %in% quad_var) vars_name[j] <- str_c(vars_name[j], '^2')
  }
    
    #add variable names together to load into model
    vars <- paste(vars, collapse = " + ")
    
    #add variable names together for readability
    vars_name <- paste(vars_name, collapse = " + ")
    
    #run model if no convergence error!
    if(is.error(eval(parse(text = str_c('lme(pirgd ~ lc + ', vars, ', random = ~ 1 | id, 
                                        data = dat3, method = "ML")')))) == F){
      m <- eval(parse(text = str_c('lme(pirgd ~ lc + ', vars, ', random = ~ 1 | id, 
                                   data = dat3, method = "ML")')))
      
      #save AIC
      precip_df <- rbind(precip_df, data.frame(var = vars_name, aic = AIC(m)))
      
  } else print(str_c(vars_name, ' did not converge')) #save names of models that did not converge
  }
}

#clean up
rm(m, i, j, vars, combos, df, vars_name)

####################
###RUN VEG COMBOS###
####################

#create table of all temp combos to evaluate
#remove bare ground due to high correlation with other RAP variables
combos <- as.data.frame(expand.grid(ann_forb_rap = 0:1, perenn_forb_rap = 0:1, tree_rap = 0:1))

#remove first row with no variables
combos <- combos[-1,]

#allocate dataframe for AIC results
veg_df <- data.frame(var = character(), aic = numeric())

#loop through all combos and run model
for(i in 1:nrow(combos)){
  
  #load row with combination
  df <- combos[i,]
  
  #allocate character vector for variable names
  vars = rep("NA", ncol(df))
  
  #loop through columns to extract variable names
  for(j in 1:ncol(df)) {
    if(df[1,j] == 1) {vars[j] <- colnames(df)[j]}
  }
  
  #remove missing names
  vars <- vars[!(vars %in% "NA")]
  
  #save variable names in readable format
  vars_name <- vars
  
  #change names to quadratic if applicable
  for(j in 1:length(vars)){
    if(vars[j] %in% quad_var) vars[j] <- str_c('poly(', vars[j], ', 2, raw = TRUE)')
  }
  
  #change names to quadratic readable if applicable
  for(j in 1:length(vars_name)){
    if(vars_name[j] %in% quad_var) vars_name[j] <- str_c(vars_name[j], '^2')
  }
  
  #add variable names together to load into model
  vars <- paste(vars, collapse = " + ")
  
  #add variable names together for readability
  vars_name <- paste(vars_name, collapse = " + ")
  
  #run model
  m <- eval(parse(text = str_c('lme(pirgd ~ lc + ', vars, ', random = ~ 1 | id, 
                               data = dat3, method = "ML")')))
    
  #save AIC
  veg_df <- rbind(veg_df, data.frame(var = vars_name, aic = AIC(m)))
    
}

#clean up
rm(m, i, j, vars, combos, df, vars_name)


#########################
###RUN ALL BEST COMBOS###
#########################

#for this one I just removed sage_homer b/c the best veg model drops that variable
#but now we've also dropped year 2012...

#test null model and the combinations of the best category models

m_null <- lme(pirgd ~ lc, random = ~ 1 | id, data = dat3, method = 'ML')

m_temp <- lme(pirgd ~ lc + poly(gdd_elev, 2, raw = TRUE), 
              random = ~ 1 | id, data = dat3, method = 'ML')

m_precip <- lme(pirgd ~ lc + poly(rain_elev, 2, raw = TRUE) + poly(pdsi_mar_apr_min, 2, raw = TRUE) + 
                  poly(snowmelt, 2, raw = TRUE), random = ~ 1 | id, data = dat3, method = 'ML')

m_veg <- lme(pirgd ~ lc + poly(ann_forb_rap, 2, raw = TRUE) + 
               perenn_forb_rap + tree_rap, random = ~ 1 | id, data = dat3, method = 'ML')

m_temp_precip <- lme(pirgd ~ lc + poly(gdd_elev, 2, raw = TRUE) +
                       poly(rain_elev, 2, raw = TRUE) + poly(pdsi_mar_apr_min, 2, raw = TRUE) + 
                      poly(snowmelt, 2, raw = TRUE), 
                     random = ~ 1 | id, data = dat3, method = 'ML')

m_temp_veg <- lme(pirgd ~ lc + poly(gdd_elev, 2, raw = TRUE) + poly(chill_jan_apr, 2, raw = TRUE) +
                    poly(ann_forb_rap, 2, raw = TRUE) + 
                    perenn_forb_rap + tree_rap, random = ~ 1 | id, data = dat3, method = 'ML')

m_precip_veg <- lme(pirgd ~ lc + poly(rain_elev, 2, raw = TRUE) + poly(pdsi_mar_apr_min, 2, raw = TRUE) + 
                      poly(snowmelt, 2, raw = TRUE) + 
                      poly(ann_forb_rap, 2, raw = TRUE) + 
                      perenn_forb_rap + tree_rap, random = ~ 1 | id, data = dat3, method = 'ML')

m_temp_precip_veg <- lme(pirgd ~ lc + poly(gdd_elev, 2, raw = TRUE) +
                           poly(rain_elev, 2, raw = TRUE) + poly(pdsi_mar_apr_min, 2, raw = TRUE) + 
                           poly(snowmelt, 2, raw = TRUE) +
                           poly(ann_forb_rap, 2, raw = TRUE) + 
                           perenn_forb_rap + tree_rap, random = ~ 1 | id, data = dat3, method = 'ML')

#create output df
combo_df <- data.frame(model = character(), aic = numeric())

#rbind all results
combo_df <- rbind(combo_df, data.frame(model = 'm_null', aic = AIC(m_null)),
                  data.frame(model = 'm_temp', aic = AIC(m_temp)),
                  data.frame(model = 'm_precip', aic = AIC(m_precip)),
                  data.frame(model = 'm_veg', aic = AIC(m_veg)),
                  data.frame(model = 'm_temp_precip', aic = AIC(m_temp_precip)),
                  data.frame(model = 'm_temp_veg', aic = AIC(m_temp_veg)),
                  data.frame(model = 'm_precip_veg', aic = AIC(m_precip_veg)),
                  data.frame(model = 'm_temp_precip_veg', aic = AIC(m_temp_precip_veg)))

#re-check correlations as well!!!!

###########################
###BEST MODEL VALIDATION###
###########################

#the best model includes variables from all three categories
#now need to check the model for homogeneity of variance and normality.

#re run model using REML
m_temp_precip_veg <- lme(pirgd ~ lc + poly(gdd_elev, 2, raw = TRUE) +
                           poly(rain_elev, 2, raw = TRUE) + poly(pdsi_mar_apr_min, 2, raw = TRUE) + 
                           poly(snowmelt, 2, raw = TRUE) +
                           poly(ann_forb_rap, 2, raw = TRUE) + 
                           perenn_forb_rap + tree_rap, random = ~ 1 | id, data = dat3)

#plot residuals of best model
e2 <- resid(m_temp_precip_veg, type = 'normalized')
f2 <- fitted(m_temp_precip_veg)

#set up output plot
jpeg('output/model_validation/m_temp_precip_veg_resid_dist.jpeg', width = 900, height = 600)

#test for normality
hist(e2, xlab = 'Residuals', main = 'Normal Distribution?', breaks = 100)

#dev off
dev.off()

#set up output plot
jpeg('output/model_validation/m_temp_precip_veg_resid_vs_fitted.jpeg', width = 900, height = 600)

#test for equal variance
plot(x = f2, y = e2, xlab = 'Fitted values', ylab = 'Residuals', main = 'Homogenity of Variance?')

#dev off
dev.off()

#plot residuals against each explanatory variable to check for patterns

#set up output plot
jpeg('output/model_validation/m_temp_precip_veg_resid_lc.jpeg', width = 900, height = 600)

#box plot of landcover
boxplot(e2 ~ lc, data = dat3, main = 'Landcover', ylab = 'Residuals', xlab = 'Landcover')

#dev off
dev.off()

#plot all other explanatory variables against residuals
#create list of variables
vars <- c('gdd_elev', 'rain_elev', 'pdsi_mar_apr_min',
          'snowmelt', 'ann_forb_rap',
          'perenn_forb_rap', 'tree_rap')

#loop through
for(var in vars){
  
  #set up output plot
  jpeg(str_c('output/model_validation/m_temp_precip_veg_resid_', var, '.jpeg'), width = 900, height = 600)
  
  #plot variable against residuals
  plot(x = dat3[, var], y = e2, ylab = 'Residuals', xlab = var, main = var)
  
  
  #dev off
  dev.off()
}

#run anova on the variance per data point to ensure the variance
#is equal for each data point

#extracts the residuals and places them in a new column in our original data table
dat3$resid <- e2

#creates a new column with the absolute value of the residuals
dat3$abs_resid <-abs(dat3$resid)

#squares the absolute values of the residuals to provide the more robust estimate
dat3$resid_sq <- dat3$abs_resid^2 

#ANOVA of the squared residuals
m_resid <- lm(resid_sq ~ id, data = dat3) 

#displays the results
anova(m_resid) 

#very significant means we don't have equal variance...

#compute r^2 as a measure of goodness of fit. use lme4 package fit to run correct functions
m_temp_precip_veg_lme4 <- lmer(pirgd ~ lc + poly(gdd_elev, 2, raw = TRUE) +
                           poly(rain_elev, 2, raw = TRUE) + poly(pdsi_mar_apr_min, 2, raw = TRUE) + 
                           poly(snowmelt, 2, raw = TRUE) +
                           poly(ann_forb_rap, 2, raw = TRUE) + 
                           perenn_forb_rap + tree_rap + (1 | id), dat3)

#compute model checks
performance::check_model(m_temp_precip_veg_lme4)
performance::model_performance(m_temp_precip_veg_lme4)

#calculate VIF
vif <- vif(m_temp_precip_veg)

#square the GVIF DF value to ensure < 5
vif <- vif[,3]^2
vif <- data.frame(vif)

############################
###STORE MODEL PARAMETERS###
############################

s <- summary(m_temp_precip_veg)
fitted <- as.data.frame(s$tTable)

#################
###MODEL PLOTS###
#################

#set up output plot
jpeg('output/model_validation/predicted_vs_observed_pirgd.jpeg', width = 600, height = 600)

#actual vs. predicted pirgd values
plot(x = dat3$pirgd, y = f2, xlab = 'Observed', ylab = 'Predicted', main = 'Predicted vs. Observed PIRGd',
     xlim = c(25, 250), ylim = c(25, 250))
abline(1,1, col = 'blue')

#dev.off
dev.off()

#re-run model with original values to see correct effect size
m_original <- lme(pirgd ~ lc + poly(gdd_elev, 2, raw = TRUE) +
                           poly(rain_elev, 2, raw = TRUE) + poly(pdsi_mar_apr_min, 2, raw = TRUE) + 
                           poly(snowmelt, 2, raw = TRUE) +
                           poly(ann_forb_rap, 2, raw = TRUE) + 
                           perenn_forb_rap + tree_rap, random = ~ 1 | id, data = dat2)

#set up output plot
jpeg('output/model_validation/variable_effect_size.jpeg', width = 1600, height = 902)

#plot all effects using correct effect size
plot(effects::predictorEffects(m_original))

#dev.off
dev.off()

#set up output plot
jpeg('output/model_validation/pirgd_vs_landcover.jpeg', width = 900, height = 600)

#plot distribution of LC classes vs. pirgd to make sure the LC effect makes sense
boxplot(pirgd ~lc, data = dat3, main = 'PIRGd vs Landcover', ylab = 'PIRGd', xlab = 'Landcover (Shrub, Herb, Evergreen)')

#dev.off
dev.off()

#create plot of residual size vs point on grid
#first, average residual size by point
#create new df with ID and resid
resid_df <- dat3[,c('id', 'resid')]

#get mean residuals per point
resid_df <- aggregate(resid_df$resid, list(resid_df$id), mean)

#add mean residuals to grid of points
grid@data$resid <- resid_df$x

#load wy shapefile
wy <- readOGR('/Volumes/SSD/climate_effects/reference/wyoming.shp')

#load lc
lc <- raster('/Volumes/SSD/climate_effects/landcover/five_class_landcover_wy_laea_2016.tif')
lc[lc == 3] <- NA
lc[lc == 5] <- NA
lc[lc == 4] <- 3

#set up output plot
jpeg('output/model_validation/resid_size_sampling_grid.jpeg', width = 715, height = 495)

#plot grid points using residual size
plot(lc, legend=FALSE, col=c("coral3", "papayawhip", "forestgreen"), xaxt='n', yaxt='n',
     main = "Average Residual Size on Sampling Grid") 
par(xpd=TRUE)
legend("bottom", legend=c("Shrub", "Herb", "Evergreen"),
       fill=c("coral3","papayawhip", "forestgreen"),horiz = TRUE, inset=-0.175)
plot(wy, add = T)
plot(grid, pch = 1, add = T, cex = .2)
plot(grid, pch = 19, cex = abs(grid@data$resid)*20, add = T)

#dev.off
dev.off()

#create a new df to look at elevation categories
dat_pdsi <- data.frame(pdsi = dat2$pdsi_mar_apr_min, pirgd = dat2$pirgd, dem = dat2$dem, lc = dat2$lc)

#separate elevations
elev <- quantile(dat_pdsi$dem, probs = c(1/3, 2/3))

#create factor variable for elevations
dat_pdsi$elev <- dat_pdsi$dem
dat_pdsi$elev[dat_pdsi$elev < elev[1]] <- 1
dat_pdsi$elev[dat_pdsi$elev >= elev[1] & dat_pdsi$elev < elev[2]] <- 2
dat_pdsi$elev[dat_pdsi$elev >= elev[2]] <- 3
dat_pdsi$elev <- as.factor(dat_pdsi$elev)

#set elev levels to elevation ranges
levels(dat_pdsi$elev) <- c('Less than 1645 m', '1645 - 2298 m', 'Above 2298 m')

#set lc levels to categories
levels(dat_pdsi$lc) <- c('Shrub', 'Herb', 'Evergreen')

#set up output plot
jpeg('output/model_validation/pdsi_pirgd_lc.jpeg', width = 715, height = 495)

#why does pdsi have a strong quadratic curve?
#check to see if pdsi varies at different elevation, or by landcovers...
ggplot(dat_pdsi, aes(x = pdsi, y = pirgd, color = lc)) + geom_point(size = .2) +
  geom_smooth() + theme_classic()

#dev.off
dev.off()

#set up output plot
jpeg('output/model_validation/pdsi_pirgd_elev.jpeg', width = 715, height = 495)

#plot!
ggplot(dat_pdsi, aes(x = pdsi, y = pirgd, color = elev)) + geom_point(size = .2) +
  geom_smooth() + theme_classic()

#dev.off
dev.off()


#code to rescale scaled values back to normal
#don't need anymore
#r*attr(xs,'scaled:scale') + attr(xs, 'scaled:center')

