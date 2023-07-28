rm(list=ls())

library(lubridate)
library(pacman)

source("HighstatLibV10.R")

#load data
allYrs <- read.csv("Sensitivity/Data/LT_siteYr_AllData_sensitivity.csv", header=T)

# choose which country for this task
TaskID <- read.csv("Sensitivity/Data/LT_ResponseTrends_TaskIDs_sensitivity.csv", as.is = T)
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

hist(allYrs$shannonsH)

hist(allYrs$spp_rich_rare)
hist(sqrt(allYrs$spp_rich_rare))

hist(allYrs$spp_richness)
hist(sqrt(allYrs$spp_richness))

hist(allYrs$turnover)
hist(allYrs$turnover^3)

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

# run transformation
if(myResponse %in% c("abundance", "ept_abundance", "diptera_abundance", 
                     "insect_abundance", "mollusc_abundance", 
                     "annelid_abundance", "crustacea_abundance")){
  
  allYrs$Response <- log10(allYrs$Response + 1)

}else if(myResponse %in% "turnover"){
  
  allYrs$Response <- allYrs$Response^3
  
}else if(myResponse %in% c("spp_richness", "spp_rich_rare", 
                           "insect_spp_richness","ept_spp_richness", 
                           "diptera_spp_richness", "mollusc_spp_richness", 
                           "annelid_spp_richness", "crustacea_spp_richness")){
  
  allYrs$Response <- sqrt(allYrs$Response)
  
}

#not transformed: shannonH
hist(allYrs$Response)

#### two-stage models ####
### fitting gls #####
library(nlme)
library(mgcv)

# write a function to consider day of year if sampling is more than 30 days apart
fitGLSModel <- function(my_data){
  
  #centre Year - helps model convergence to center variables for the model
  my_data$cYear <- my_data$year - median(my_data$year)
  
  #or just have as an index starting from 1
  my_data$iYear <- my_data$year - min(my_data$year)+1
  
  #scale day of year
  my_data$cday_of_year <- (my_data$day_of_year - mean(my_data$day_of_year))/sd(my_data$day_of_year)
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(my_data$day_of_year)-min(my_data$day_of_year)
  
  if(length(my_data$Response[is.na(my_data$Response)])>0) {
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

saveRDS(trends, file=paste0("Sensitivity/Outputs/Site_trends/trends__",myResponse,"__",myCountry,".RDS"))
#write.csv(trends, file=paste0("outputs/trends__",myResponse,"__",myCountry,".csv"))

# # ONLY USE FOR CRUSTACEA
# ## Crustacea were not found at all sites and there was an issue with model convergence because of this.
# ## Thus, use this code for crustacea only
# 
# # make new vector for crustacea because of convergence issues
# library(dplyr)
# allsites <- sort(unique(allYrs$site_id))
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
# saveRDS(trends, file=paste0("Sensitivity/Outputs/Site_trends/trends__",myResponse,"__",myCountry,".RDS"))

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
