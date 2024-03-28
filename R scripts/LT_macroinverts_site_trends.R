rm(list=ls())

library(lubridate)
library(pacman)

source("Additional functions/HighstatLibV10.R")

#load data
d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

# choose which country for this task
TaskID <- read.csv("Data/LT_ResponseTrends_TaskIDs.csv", as.is = T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myCountry <- TaskID$country[which(TaskID$TaskID==task.id)]

#choose which response for this task
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]
allYrs$Response <- allYrs[, myResponse]
hist(allYrs$Response)

### Check distributions of response variables
# Main responses
hist(allYrs$abundance)
hist(log10(allYrs$abundance+1))

hist(allYrs$E10)

hist(allYrs$F_turnover)
hist(log10(allYrs$F_turnover+0.01))

hist(allYrs$FDis)
hist(allYrs$FDis^3)

hist(allYrs$FEve)
hist(allYrs$FEve^2)

hist(allYrs$FRed)

hist(allYrs$FRic)
summary(allYrs$FRic)
hist(log10(allYrs$FRic + 0.28))

hist(allYrs$shannonsH)

hist(allYrs$spp_rich_rare)
hist(sqrt(allYrs$spp_rich_rare))

hist(allYrs$spp_richness)
hist(sqrt(allYrs$spp_richness))

hist(allYrs$turnover)
hist(allYrs$turnover^3)

hist(allYrs$FRic.SES)
hist(allYrs$FEve.SES)
hist(allYrs$FDis.SES)

# taxonomic level responses
hist(allYrs$ept_spp_richness)
hist(sqrt(allYrs$ept_spp_richness))

hist(allYrs$ept_abundance)
hist(log10(allYrs$ept_abundance+1))

hist(allYrs$diptera_spp_richness)
hist(sqrt(allYrs$diptera_spp_richness))

hist(allYrs$diptera_abundance)
hist(log10(allYrs$diptera_abundance + 1))

hist(allYrs$insect_spp_richness)
hist(sqrt(allYrs$insect_spp_richness))

hist(allYrs$insect_abundance)
hist(log10(allYrs$insect_abundance + 1))

hist(allYrs$mollusc_spp_richness)
hist(sqrt(allYrs$mollusc_spp_richness))

hist(allYrs$mollusc_abundance)
hist(log10(allYrs$mollusc_abundance + 1))

hist(allYrs$annelid_spp_richness)
hist(sqrt(allYrs$annelid_spp_richness))

hist(allYrs$annelid_abundance)
hist(log10(allYrs$annelid_abundance+1))

hist(allYrs$crustacea_spp_richness)
hist(sqrt(allYrs$crustacea_spp_richness))

hist(allYrs$crustacea_abundance)
hist(log10(allYrs$crustacea_abundance + 1))

# transform variables that are right-skewed
# first run transformation
# if(myResponse %in% c("abundance", "spp_richness", "spp_rich_rare",
#                      "ept_spp_richness", "ept_abundance",
#                      "diptera_spp_richness", "diptera_abundance",
#                      "insect_spp_richness", "insect_abundance",
#                      "mollusc_spp_richness", "mollusc_abundance",
#                      "annelid_spp_richness", "annelid_abundance")){
#   
#   allYrs$Response <- log10(allYrs$Response + 1)
#   
# }else if(myResponse %in% c("E10","F_turnover", "FRed", "FRic")){
#   
#   allYrs$Response <- log10(allYrs$Response + 0.01)
#   
# }else if(myResponse %in% c("FDis", "turnover")){
#   
#   allYrs$Response <- allYrs$Response^2
#   
# }

# # second run transformation
# if(myResponse %in% c("abundance", "ept_abundance", "diptera_abundance", "insect_abundance", "mollusc_abundance", "annelid_abundance",
#                      "annelid_spp_richness")){
#   
#   allYrs$Response <- log10(allYrs$Response + 1)
#   
# }else if(myResponse %in% c("F_turnover", "FRic")){
#   
#   allYrs$Response <- log10(allYrs$Response + 0.01)
#   
# }else if(myResponse %in% c("FDis", "turnover")){
#   
#   allYrs$Response <- allYrs$Response^2
#   
# }

# # third run transformation
# if(myResponse %in% c("abundance", "ept_abundance", "diptera_abundance", "insect_abundance", "mollusc_abundance", "annelid_abundance")){
#   
#   allYrs$Response <- log10(allYrs$Response + 1)
#   
# }else if(myResponse %in% c("F_turnover", "FRic")){
#   
#   allYrs$Response <- log10(allYrs$Response + 0.01)
#   
# }else if(myResponse %in% c("FDis", "turnover")){
#   
#   allYrs$Response <- allYrs$Response^2
#   
# }else if(myResponse %in% c("ept_spp_richness", "diptera_spp_richness", "insect_spp_richness", "mollusc_spp_richness", "annelid_spp_richness")){
#   
#   allYrs$Response <- sqrt(allYrs$Response)
#   
# }

# fourth run transformation
if(myResponse %in% c("abundance", "ept_abundance", "diptera_abundance", 
                     "insect_abundance", "mollusc_abundance", 
                     "annelid_abundance", "crustacea_abundance")){
  
  allYrs$Response <- log10(allYrs$Response + 1)
  
}else if(myResponse %in% "F_turnover"){
  
  allYrs$Response <- log10(allYrs$Response + 0.01)
  
}else if(myResponse %in% "FRic") {
  
  allYrs$Response <- log10(allYrs$Response + 0.28)
  
}else if(myResponse %in% "FEve"){
  
  allYrs$Response <- allYrs$Response^2
  
}else if(myResponse %in% c("FDis", "turnover")){
  
  allYrs$Response <- allYrs$Response^3
  
}else if(myResponse %in% c("spp_richness", "spp_rich_rare", 
                           "insect_spp_richness","ept_spp_richness", 
                           "diptera_spp_richness", "mollusc_spp_richness", 
                           "annelid_spp_richness", "crustacea_spp_richness")){
  
  allYrs$Response <- sqrt(allYrs$Response)
  
}

#not transformed: FEve, shannonH, FRic.SES, FEve.SES, FDis.SES
hist(allYrs$Response)

#### two-stage models ####
### fitting gls #####
library(nlme)
library(mgcv)

# write a function to consider day of year if sampling is more than 30 days apart
fitGLSModel <- function(my_data){
  
  #centre Year - helps model convergence to center variables for the model
  my_data$cYear <- my_data$year_wMissing - median(my_data$year_wMissing)
  
  #or just have as an index starting from 1
  my_data$iYear <- my_data$year_wMissing - min(my_data$year_wMissing)+1
  
  #scale day of year
  my_data$cday_of_year <- (my_data$day_of_year - mean(my_data$day_of_year))/sd(my_data$day_of_year)
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(my_data$day_of_year)-min(my_data$day_of_year)
  
#   #Response sum
#   positiveData <- subset(my_data, !is.na(my_data$Response)|my_data$Response!=0)
#   sumResponse = sum(my_data$Response)
#   nuYears = length(unique(positiveData$year_wMissing))
#   
#   if(all(is.na(my_data$Response))|sumResponse==0|nuYears==1) { # from Ellen's paper
  
  if(length(my_data$Response[is.na(my_data$Response)])>0) { # code made with Francesca
    my_data <- subset(my_data, !is.na(Response))
  } else { 
    my_data <- my_data
  }
  
  if(maxDiffDays < 30) {
    myformula <- gls(Response ~ cYear, correlation = corAR1(form =~ iYear), data = my_data)
  } else {
    myformula <- gls(Response ~ cday_of_year + cYear, correlation = corAR1(form =~ iYear), data = my_data)
  }
  
  #fit model with gls
  gls_model <- myformula
  
  #extract model fits
  modelSummary <- summary(gls_model)$tTable["cYear", c("Value", "Std.Error", "p-value")]
  modelFits <- data.frame(estimate = modelSummary["Value"],
                          se = modelSummary["Std.Error"],
                          pval = modelSummary["p-value"])
  
  return(modelFits)
}


#the model is called in the function below
#apply function to an example dataset
est <- fitGLSModel(allYrs[which(allYrs$site_id=="LTR1"),])
est

#loop for all sites
allsites <- sort(unique(allYrs$site_id))

trends <- lapply(allsites, function(x){
  fitGLSModel(subset(allYrs, site_id == x))
})

trends <- data.frame(do.call(rbind, trends))
trends$siteID <- allsites
rownames(trends) <- 1:41
trends

saveRDS(trends, file=paste0("Outputs/Site_trends/trends__",myResponse,"__",myCountry,".RDS"))

# # ONLY USE FOR CRUSTACEA
# ## Crustacea were not found at all sites and there was an issue with model convergence because of this. 
# ## Thus, use this code for crustacea only
# 
# # make new vector for crustacea because of convergence issues
# library(dplyr)
# allsites_cru <- allsites[!allsites %in% c("LTR1319", "LTR133", "LTR327")] # these sites had no crustacea
# 
# # calculate site trends
# trends <- lapply(allsites_cru, function(x){
#   fitGLSModel(subset(allYrs, site_id == x))
# })
# 
# trends <- data.frame(do.call(rbind, trends))
# trends$siteID <- allsites_cru
# rownames(trends) <- 1:38
# trends
# 
# # add missing sites to trends dataframe
# missing_sites <- data.frame(
#   estimate = rep(NA, 3),
#   se = rep(NA, 3),
#   pval = rep(NA, 3),
#   siteID = c("LTR1319", "LTR133", "LTR327")
# )
# trends <- rbind(trends, missing_sites) # combine dataframes
# trends <- trends[order(trends$siteID), ] # reorder site_id
# rownames(trends) <- NULL # reset rownames
# trends
# 
# saveRDS(trends, file=paste0("Outputs/Site_trends/trends__",myResponse,"__",myCountry,".RDS"))

##### CLEAN UP --------------------
library(pacman)
# Clear data
rm(list = ls())  # Removes all objects from environment
# Clear packages
p_unload(all)  # Remove all contributed packages
# Clear plots
graphics.off()  # Clears plots, closes all graphics devices
# Clear console
cat("\014")  # Mimics ctrl+L
# Clear mind :)
