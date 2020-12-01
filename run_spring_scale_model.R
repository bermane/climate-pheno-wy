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
library(performance)

#set wd
setwd('/Volumes/SSD/climate_effects')

#load source functions
source("reference/HighstatLibV10.R")

#load workspace
load('/Users/Ediz/Team Braintree Dropbox/Ethan Berman/R Projects/climate-pheno-wy/pirgd_dat_2020_10_18.RData')

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

#remove all veg layers since we are moving forward without them
dat2 <- select(dat2, -c(herb_homer_ann, sage_homer_ann, shrub_homer_ann, herb_homer,
                        sage_homer, shrub_homer, ann_forb_rap, bare_ground_rap,
                        perenn_forb_rap, shrub_rap, tree_rap, ann_perenn_forb_rap))

#remove vapor pressure variables
dat2 <- select(dat2, -c(vp_min_elev, vp_avg_elev, vp_min_jan_apr, vp_avg_jan_apr))

#remove time1-PIRGd variables
dat2 <- select(dat2, -c(solar_jan_pirgd, rain_oct_pirgd, rain_jan_pirgd,
                        gdd_jan_pirgd, pdsi_jan_pirgd_mean, pdsi_jan_pirgd_min,
                        pdsi_jan_pirgd_med))

#remove vpd by elev. keep in values from jan-apr
dat2 <- select(dat2, -c(vpd_tmax_mean_elev, vpd_tmin_mean_elev, vpd_tavg_mean_elev))


#remove any remaining missing data points for spring scale or snowmelt
#2 points for pirgd
#dat2 <- dat2[is.na(dat2$pirgd) == 0,]
dat2 <- dat2[is.na(dat2$ss) == 0,]
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
levels(dat2$lc) <- c('shrub', 'herb', 'evergreen')

#loop through all covariates in table
#check to make sure correct starting column -- start with dem!!
colnames(dat2)[7]
ncol(dat2)
for(col_id in 7:35){
  
  #set varible name
  vari <- colnames(dat2)[col_id]
  
  #fit a linear mixed model with landcover, random effect, and variable
  m <- eval(substitute(lme(ss ~ lc + variable, random = ~ 1 | id, data = dat2),
                         list(variable = as.name(vari))))

  #plot residuals
  e2 <- resid(m, type = 'normalized')
  f2 <- fitted(m)
  
  #set up output plot
  jpeg(str_c('output/spring_scale/norm_equal_var_test_plots/' ,vari, '.jpeg'), width = 900, height = 600)
  
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
for(col_id in 7:35){
  
  #set variable name
  vari <- colnames(dat2)[col_id]
  
  #fit a linear mixed model with landcover, random effect, and variable
  m <- eval(substitute(lme(ss ~ lc + variable, random = ~ 1 | id, 
                           data = dat2, method = 'ML'),
                       list(variable = as.name(vari))))
  
  #get summary
  s <- summary(m)
    
  #fit with quadratic term  
  m2 <- eval(substitute(lme(ss ~ lc + poly(variable, 2, raw = TRUE), random = ~ 1 | id, 
                            data = dat2, method = 'ML'),
                       list(variable = as.name(vari))))
  
  #get summary
  s2 <- summary(m2)
  
  quad_df <- rbind(quad_df, data.frame(var = vari, quad_coef = s2$tTable[str_c('poly(', vari, ', 2, raw = TRUE)2'),"Value"],
                                       quad_95 = (s2$tTable[str_c('poly(', vari, ', 2, raw = TRUE)2'),"Std.Error"])*1.96,
                                       quad_aic = s2$AIC, lin_aic = s$AIC))
  
}

#FOR NOW ONLY USE AIC TO DECIDE WHICH VARIABLES TO CARRY FORWARD AS QUADRATICS
#set low and high confidence intervals to test overlap with zero
#quad_df$quad_low <- quad_df$quad_coef - quad_df$quad_95
#quad_df$quad_hi <- quad_df$quad_coef + quad_df$quad_95

#subset dataframe to variables where quad does not overlap zero
#quad_df2 <- quad_df[!(quad_df$quad_low < 0 & quad_df$quad_hi > 0),]

#get null model AIC to test against quadratic and linear AIC
AIC(lme(ss ~ lc, random = ~ 1 | id, 
        data = dat2, method = 'ML'))

#subset dataframe to variables where quad AIC is lower than lin AIC but diff > 2
quad_df2 <- quad_df[quad_df$quad_aic < quad_df$lin_aic & abs(quad_df$quad_aic - quad_df$lin_aic) > 2,]

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
#7 to 35

#rescale variables between 0 and 1
for(i in 7:35){
  dat3[,i] <- scale(dat3[,i])
}

#test fit a linear model with only elevation
m <- lm(ss ~ lc + dem, dat3)

#test fit a linear mixed model with only elevation
m2 <- lme(ss ~ lc + dem, random = ~ 1 | id, data = dat3, method = 'ML')

#check AIC of the two
AIC(m)
AIC(m2)

###########################
###LOOK INTO CORRELATION###
###########################

#check into correlation between all variables
#compute correlation matrix
#make sure to check including all columns!!!
cor_m <- cor(dat3[,7:35])
cor_m <- round(cor_m, 2)

#create output table of overall correlations that are significant
#allocate df
cor_df <- data.frame(var1 = character(), var2 = character(), cor = numeric())

#write to csv file full correlation results
write.csv(cor_df, file = 'output/spring_scale/full correlation results - 2020-10-18.csv')

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

#write to csv file significant correlation results
write.csv(cor_df, file = 'output/spring_scale/significant correlation results - 2020-10-18.csv')

########################################################
###TEST VARIABLES FOR NON-ZERO and NON-SIGN SWITCHING###
########################################################

#which columns do we want to test? everything after landcover
colnames(dat3)[7]
ncol(dat3)
#7 to 35

#allocate dataframe to save coefficients
uni_df <- data.frame(var = character(), beta = numeric(), beta_sq = numeric(), se_95 = numeric(), se_95_sq = numeric(),
                     sign = character(), sign_sq = character(), aic = numeric(), mar_r2 = numeric())

#loop through variables
for(i in 7:35){
  
  #get column to test
  var <- colnames(dat3)[i]
  
  #check if we want to model quadratic, otherwise linear term only
  if(var %in% quad_var){
    
    #fit a linear mixed model with base + quad variables -- only if converges!
    
    if(is.error(eval(substitute(lme(ss ~ lc + poly(variable, 2, raw = TRUE), random = ~ 1 | id, 
                                    data = dat3, method = 'ML'),
                                list(variable = as.name(var))))) == F){

      m <- eval(substitute(lme(ss ~ lc + poly(variable, 2, raw = TRUE), random = ~ 1 | id, 
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
      mar_r2 <- suppressWarnings(performance::r2(m)$R2_marginal)
      
      #find whether signs are zero/positive/negative
      if((beta - se_95 < 0) & (beta + se_95 < 0)) sign <- 'negative'
      if((beta - se_95 > 0) & (beta + se_95 > 0)) sign <- 'positive'
      if((beta - se_95 < 0) & (beta + se_95 > 0)) sign <- 'zero'
      
      if((beta_sq - se_95_sq < 0) & (beta_sq + se_95_sq < 0)) sign_sq <- 'negative'
      if((beta_sq - se_95_sq > 0) & (beta_sq + se_95_sq > 0)) sign_sq <- 'positive'
      if((beta_sq - se_95_sq < 0) & (beta_sq + se_95_sq > 0)) sign_sq <- 'zero'
      
      uni_df <- rbind(uni_df, data.frame(var = var, beta = beta, beta_sq = beta_sq, se_95 = se_95, 
                                         se_95_sq = se_95_sq, sign = sign, sign_sq = sign_sq, 
                                         aic = aic, mar_r2 = mar_r2)) 
      
    } else print(str_c(var, ' did not converge')) #save names of models that did not converge
    
  } else{ #run as linear term only!
    
    #fit a linear mixed model with base + linear variables if no convergence error!
    if(is.error(eval(substitute(lme(ss ~ lc + variable, random = ~ 1 | id, 
                                    data = dat3, method = 'ML'),
                                list(variable = as.name(var))))) == F){
      m <- eval(substitute(lme(ss ~ lc + variable, random = ~ 1 | id, 
                               data = dat3, method = 'ML'),
                           list(variable = as.name(var))))
      
      #get summary of model
      s <- summary(m)
      
      #extract values we want from model summary
      beta <- s$tTable[var,"Value"]
      se_95 <- s$tTable[var,"Std.Error"]*1.96
      aic <- s$AIC
      mar_r2 <- suppressWarnings(performance::r2(m)$R2_marginal)
      
      #find whether signs are zero/positive/negative
      if((beta - se_95 < 0) & (beta + se_95 < 0)) sign <- 'negative'
      if((beta - se_95 > 0) & (beta + se_95 > 0)) sign <- 'positive'
      if((beta - se_95 < 0) & (beta + se_95 > 0)) sign <- 'zero'
      
      uni_df <- rbind(uni_df, data.frame(var = var, beta = beta, beta_sq = NA, se_95 = se_95, 
                                         se_95_sq = NA, sign = sign, sign_sq = 'NA', 
                                         aic = aic, mar_r2 = mar_r2)) 

      
    } else print(str_c(var, ' did not converge')) #save names of models that did not converge
  }
}

#clean up
rm(m, s, beta, beta_sq, se_95, aic, sign, sign_sq, var, i, se_95_sq, mar_r2)

#write to csv file
write.csv(uni_df, file = 'output/spring_scale/univariate results - 2020-10-18.csv')

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
#at this point only moving forward with gdd_jan_apr so actually not needed for this category
combos <- as.data.frame(expand.grid(gdd_mar_may = 0:1, chill_jan_may = 0:1))

#remove first row with no variables and combos with dem except dem + gdd_elev
combos <- combos[-1,]

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
  m <- eval(parse(text = str_c('lme(ss ~ lc + ', vars, ', random = ~ 1 | id, 
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
combos <- as.data.frame(expand.grid(pdsi_mar_apr_min = 0:1, rain_mar_may = 0:1,
                                    snow_oct_apr = 0:1, vpd_tavg_mean_jan_apr = 0:1))

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
    if(is.error(eval(parse(text = str_c('lme(ss ~ lc + ', vars, ', random = ~ 1 | id, 
                                        data = dat3, method = "ML")')))) == F){
      m <- eval(parse(text = str_c('lme(ss ~ lc + ', vars, ', random = ~ 1 | id, 
                                   data = dat3, method = "ML")')))
      
      #save AIC
      precip_df <- rbind(precip_df, data.frame(var = vars_name, aic = AIC(m)))
      
  } else print(str_c(vars_name, ' did not converge')) #save names of models that did not converge
  }
}

#clean up
rm(m, i, j, vars, combos, df, vars_name)

#########################
###RUN ALL BEST COMBOS###
#########################

#test null model and the combinations of the best category models

m_null <- lme(ss ~ lc, random = ~ 1 | id, data = dat3, method = 'ML')

m_precip <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                  vpd_tavg_mean_jan_apr, random = ~ 1 | id, data = dat3, method = 'ML')

#create output df
combo_df <- data.frame(model = character(), aic = numeric())

#rbind all results
combo_df <- rbind(combo_df, data.frame(model = 'm_null', aic = AIC(m_null)),
                  data.frame(model = 'm_precip', aic = AIC(m_precip)))

#re-check correlations as well!!!!

###########################
###TEST FOR INTERACTIONS###
###########################

m_base <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                vpd_tavg_mean_jan_apr, random = ~ 1 | id, data = dat3, method = 'ML')

m_intaxn1 <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                   vpd_tavg_mean_jan_apr + poly(pdsi_mar_apr_min, 2, raw = TRUE):poly(rain_mar_may, 2, raw = TRUE), 
                 random = ~ 1 | id, data = dat3, method = 'ML')

m_intaxn2 <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                   vpd_tavg_mean_jan_apr + vpd_tavg_mean_jan_apr:poly(rain_mar_may, 2, raw = TRUE), 
                 random = ~ 1 | id, data = dat3, method = 'ML')

m_intaxn3 <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                   vpd_tavg_mean_jan_apr + vpd_tavg_mean_jan_apr:poly(pdsi_mar_apr_min, 2, raw = TRUE), 
                 random = ~ 1 | id, data = dat3, method = 'ML')

m_intaxn_1_2 <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                      vpd_tavg_mean_jan_apr + poly(pdsi_mar_apr_min, 2, raw = TRUE):poly(rain_mar_may, 2, raw = TRUE) +
                      vpd_tavg_mean_jan_apr:poly(rain_mar_may, 2, raw = TRUE), 
                    random = ~ 1 | id, data = dat3, method = 'ML')

AIC(m_base)
AIC(m_intaxn1)
AIC(m_intaxn2)
AIC(m_intaxn3)
AIC(m_intaxn_1_2)

###########################
###BEST MODEL VALIDATION###
###########################

#the best model includes variables from the two categories
#now need to check the model for homogeneity of variance and normality.

#re run model using REML
m_final <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                 vpd_tavg_mean_jan_apr + poly(pdsi_mar_apr_min, 2, raw = TRUE):poly(rain_mar_may, 2, raw = TRUE) +
                 vpd_tavg_mean_jan_apr:poly(rain_mar_may, 2, raw = TRUE), 
               random = ~ 1 | id, data = dat3)

#plot residuals of best model
e2 <- resid(m_final, type = 'normalized')
f2 <- fitted(m_final)

#set up output plot
jpeg('output/spring_scale/model_validation/resid_dist.jpeg', width = 900, height = 600)

#test for normality
hist(e2, xlab = 'Residuals', main = 'Normal Distribution?', breaks = 100)

#dev off
dev.off()

#set up output plot
jpeg('output/spring_scale/model_validation/m_final_resid_vs_fitted.jpeg', width = 900, height = 600)

#test for equal variance
plot(x = f2, y = e2, xlab = 'Fitted values', ylab = 'Residuals', main = 'Homogenity of Variance?')

#dev off
dev.off()

#plot residuals against each explanatory variable to check for patterns

#set up output plot
jpeg('output/spring_scale/model_validation/m_final_resid_lc.jpeg', width = 900, height = 600)

#box plot of landcover
boxplot(e2 ~ lc, data = dat3, main = 'Landcover', ylab = 'Residuals', xlab = 'Landcover')

#dev off
dev.off()

#plot all other explanatory variables against residuals
#create list of variables
vars <- c('pdsi_mar_apr_min', 'rain_mar_may',
          'vpd_tavg_mean_jan_apr')

#loop through
for(var in vars){
  
  #set up output plot
  jpeg(str_c('output/spring_scale/model_validation/m_final_resid_', var, '.jpeg'), width = 900, height = 600)
  
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
m_final_lme4 <- lmer(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                       vpd_tavg_mean_jan_apr + 
                       poly(pdsi_mar_apr_min, 2, raw = TRUE):poly(rain_mar_may, 2, raw = TRUE) + 
                       vpd_tavg_mean_jan_apr:poly(rain_mar_may, 2, raw = TRUE) +
                       (1 | id), dat3)

#compute model checks
performance::check_model(m_final_lme4)
performance::model_performance(m_final_lme4)

#calculate VIF
vif <- vif(m_final)

#square the GVIF DF value to ensure < 5
vif <- vif[,3]^2
vif <- data.frame(vif)

############################
###STORE MODEL PARAMETERS###
############################

s <- summary(m_final)
fitted <- as.data.frame(s$tTable)

#################
###MODEL PLOTS###
#################

#set up output plot
jpeg('output/spring_scale/model_validation/predicted_vs_observed_ss.jpeg', width = 600, height = 600)

#actual vs. predicted pirgd values
plot(x = dat3$ss, y = f2, xlab = 'Observed', ylab = 'Predicted', main = 'Predicted vs. Observed Spring Scale',
     xlim = c(0, 50), ylim = c(0, 50))
abline(1,1, col = 'blue')

#dev.off
dev.off()

#re-run model with original values to see correct effect size
m_original <- lme(ss ~ lc + poly(pdsi_mar_apr_min, 2, raw = TRUE) + poly(rain_mar_may, 2, raw = TRUE) + 
                    vpd_tavg_mean_jan_apr + 
                    poly(pdsi_mar_apr_min, 2, raw = TRUE):poly(rain_mar_may, 2, raw = TRUE) + 
                    vpd_tavg_mean_jan_apr:poly(rain_mar_may, 2, raw = TRUE), 
                  random = ~ 1 | id, data = dat2)

#set up output plot
jpeg('output/spring_scale/model_validation/variable_effect_size.jpeg', width = 1600, height = 902)

#plot all effects using correct effect size
plot(effects::predictorEffects(m_original))

#dev.off
dev.off()

#plot each effect separately!
#create list of effects variables
vars_effect <- c('lc', 'pdsi_mar_apr_min', 'rain_mar_may',
                 'vpd_tavg_mean_jan_apr')

#plot each effect using correct effect size
#manual saving is working better
#code for var 1

#plot
plot(effects::predictorEffect(vars_effect[1], m_original),
     axes = list(y = list(lim = c(0, 60), lab = 'Spring Scale'),
                 x = list(lc = list(lab = 'Landcover'))),
     main = 'Landcover Predictor Effect',
     lines=list(multiline=TRUE), confint=list(style="auto"))

#code for var 2
#calc quantiles and mean and use nice rounded values
quantile(dat2$rain_mar_may, 0.25)
mean(dat2$rain_mar_may)
quantile(dat2$rain_mar_may, 0.75)

#plot
plot(effects::predictorEffect(vars_effect[2], m_original,
                              xlevels=list(rain_mar_may = c(100, 150, 200))),
     axes = list(y = list(lim = c(0, 60), lab = 'Spring Scale'),
                 x = list(pdsi_mar_apr_min = list(lab = 'Min PDSI Mar-Apr'))),
     main = 'Min PDSI Mar-Apr Predictor Effect',
     lines=list(multiline=TRUE), confint=list(style="auto"))

#code for var 3
#calc quantiles and mean and use nice rounded values
quantile(dat2$pdsi_mar_apr_min, 0.25)
mean(dat2$pdsi_mar_apr_min)
quantile(dat2$pdsi_mar_apr_min, 0.75)

quantile(dat2$vpd_tavg_mean_jan_apr, 0.25)
mean(dat2$vpd_tavg_mean_jan_apr)
quantile(dat2$vpd_tavg_mean_jan_apr, 0.75)

#plot with 3 vpd graphs
plot(effects::predictorEffect(vars_effect[3], m_original,
                              xlevels=list(pdsi_mar_apr_min = c(-2, 0, 2),
                                           vpd_tavg_mean_jan_apr = c(260, 380, 490))),
     axes = list(y = list(lim = c(0, 60), lab = 'Spring Scale'),
                 x = list(rain_mar_may = list(lab = 'Rain Mar-May'))),
     main = 'Rain Mar-May Predictor Effect',
     lines=list(multiline=TRUE), confint=list(style="auto"))

#plot with 3 rain graphs
plot(effects::predictorEffect(vars_effect[3], m_original,
                              xlevels=list(pdsi_mar_apr_min = c(-2, 0, 2),
                                           vpd_tavg_mean_jan_apr = c(260, 380, 490))),
     axes = list(y = list(lim = c(0, 60), lab = 'Spring Scale'),
                 x = list(rain_mar_may = list(lab = 'Rain Mar-May'))),
     main = 'Rain Mar-May Predictor Effect',
     lines=list(multiline=TRUE, z.var = 'vpd_tavg_mean_jan_apr'), confint=list(style="auto"))

#code for var 4
#calc quantiles and mean and use nice rounded values
quantile(dat2$rain_mar_may, 0.25)
mean(dat2$rain_mar_may)
quantile(dat2$rain_mar_may, 0.75)

#plot
plot(effects::predictorEffect(vars_effect[4], m_original,
                              xlevels=list(rain_mar_may = c(100, 150, 200))),
     axes = list(y = list(lim = c(0, 60), lab = 'Spring Scale'),
                 x = list(vpd_tavg_mean_jan_apr = list(lab = 'Mean VPD Jan-Apr'))),
     main = 'Mean VPD Jan-Apr Predictor Effect',
     lines=list(multiline=TRUE), confint=list(style="auto"))

#set up output plot
jpeg('output/spring_scale/model_validation/ss_vs_landcover.jpeg', width = 900, height = 600)

#plot distribution of LC classes vs. ss to make sure the LC effect makes sense
boxplot(ss ~ lc, data = dat3, main = 'SS vs Landcover', ylab = 'SS', xlab = 'Landcover (Shrub, Herb, Evergreen)')

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
jpeg('output/spring_scale/model_validation/resid_size_sampling_grid.jpeg', width = 715, height = 495)

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
dat_pdsi <- data.frame(pdsi = dat2$pdsi_mar_apr_min, ss = dat2$ss, dem = dat2$dem, lc = dat2$lc)

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
jpeg('output/spring_scale/model_validation/pdsi_ss_lc.jpeg', width = 715, height = 495)

#why does pdsi have a strong quadratic curve?
#check to see if pdsi varies at different elevation, or by landcovers...
ggplot(dat_pdsi, aes(x = pdsi, y = ss, color = lc)) + geom_point(size = .2) +
  geom_smooth() + theme_classic()

#dev.off
dev.off()

#set up output plot
jpeg('output/spring_scale/model_validation/pdsi_ss_elev.jpeg', width = 715, height = 495)

#plot!
ggplot(dat_pdsi, aes(x = pdsi, y = ss, color = elev)) + geom_point(size = .2) +
  geom_smooth() + theme_classic()

#dev.off
dev.off()

#plot dot and whisker of coefficient results
#fit data frame of fitted values
m_df <- data.frame(term = rownames(fitted), estimate = fitted$Value, std.error = fitted$Std.Error)

#change names of terms
#create index first
index <- c('(Intercept)', 'herb', 'evergreen', 'rain_mar_may', 'rain_mar_may^2',
           'pdsi_mar_apr_min', 'pdsi_mar_apr_min^2', 'vpd_tavg_mean_jan_apr', 'pdsi_mar_apr_min:rain_mar_may', 
           'pdsi_mar_apr_min:rain_mar_may^2',
           'pdsi_mar_apr_min^2:rain_mar_may', 'pdsi_mar_apr_min^2:rain_mar_may^2',
           'rain_mar_may:vpd_tavg_mean_jan_apr', 'rain_mar_may^2:vpd_tavg_mean_jan_apr')

names(index) <- c('(Intercept)', 'lcherb', 'lcevergreen',
                  'poly(rain_mar_may, 2, raw = TRUE)1', 'poly(rain_mar_may, 2, raw = TRUE)2',
                  'poly(pdsi_mar_apr_min, 2, raw = TRUE)1', 'poly(pdsi_mar_apr_min, 2, raw = TRUE)2',
                  'vpd_tavg_mean_jan_apr',
                  'poly(pdsi_mar_apr_min, 2, raw = TRUE)1:poly(rain_mar_may, 2, raw = TRUE)1',
                  'poly(pdsi_mar_apr_min, 2, raw = TRUE)1:poly(rain_mar_may, 2, raw = TRUE)2',
                  'poly(pdsi_mar_apr_min, 2, raw = TRUE)2:poly(rain_mar_may, 2, raw = TRUE)1',
                  'poly(pdsi_mar_apr_min, 2, raw = TRUE)2:poly(rain_mar_may, 2, raw = TRUE)2',
                  'poly(rain_mar_may, 2, raw = TRUE)1:vpd_tavg_mean_jan_apr',
                  'poly(rain_mar_may, 2, raw = TRUE)2:vpd_tavg_mean_jan_apr')

#replace values of term variable
m_df$term <- dplyr::recode(m_df$term, !!!index)

#also want to create categories
m_df$category <- NA
m_df$category[m_df$term %in% c('(Intercept)', 'herb', 'evergreen')] <- 'Landcover'
m_df$category[m_df$term %in% c('rain_mar_may', 'rain_mar_may^2', 'pdsi_mar_apr_min',
                               'pdsi_mar_apr_min^2', 'vpd_tavg_mean_jan_apr')] <- 'Moisture'
m_df$category[is.na(m_df$category)] <- 'Interaction'

#order category variable
m_df$category <- ordered(m_df$category, levels = c('Landcover', 'Moisture', 'Interaction'))

#sort big to small effect
m_df <- m_df[order(-abs(m_df$estimate)),]

dotwhisker::dwplot(m_df,
       vline = geom_vline(xintercept = 0, colour = 'grey60', linetype = 2),
       dot_args = list(aes(color = category), size = 3),
       whisker_args = list(aes(color = category), size = 1)) + 
  theme_bw() + xlab('Coefficient Estimate') + ylab('') +
  ggtitle('Effect Size of Spring Scale Drivers') + 
  scale_colour_manual(name = 'Category', values = c('#EE6677', '#4477AA', '#AA3377')) +
  theme(legend.position = (c(0.72, 0.03)),
        legend.justification = c(0, 0),
        legend.background = element_rect(colour="grey80"),
        plot.title = element_text(size = 15), axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))


#code to rescale scaled values back to normal
#don't need anymore
#r*attr(xs,'scaled:scale') + attr(xs, 'scaled:center')

