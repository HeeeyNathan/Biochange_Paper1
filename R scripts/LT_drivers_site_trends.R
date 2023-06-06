##### load functions, packages, and external sources #####
rm(list=ls())

library(lubridate)
library(pacman)
library(vegan)
library(INLA)
library(car) # for logit()
library(xtable) # for tables
library(RColorBrewer)

source("HighstatLibV10.R")

##### read in data & format #####
d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T)
allYrs <- d1[!is.na(d1$site_id_wMissing),]

# Scale the covariates
# function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- lapply(df, function(x) if(is.numeric(x)){
    scale(x, center=TRUE, scale=TRUE)
  } else x)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df[,c(1, 2)], newd)
}

#apply function
sPredictors <- scaleVars(allYrs[, c(3, 4, 58:72)])
sPredictors <- subset(sPredictors, select = -c(site_id, year, ssite_id, syear)) # remove ID variable
sPredictors$ID <- rownames(sPredictors)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, sPredictors, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# define new variable for nutrients
# Construct PCA to check the environmental variables and their relationships
pairs(allYrs[, c(73:87)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity
pairs(allYrs[, c(75:87)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity
pairs(allYrs[, c(76, 79, 82:87)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity
allYrs_pca = princomp(na.omit(allYrs[, c(75:87)]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", xlab = "Axis.1", ylab = "Axis.2", las = 1) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$site_code, col = RivTypCols) #adds text for samples
text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$site_id, col = "red", cex = 0.6) #adds text for samples
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$ftype, col = "darkblue") #adds text for samples
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$fmodified, col = "forestgreen") #adds text for samples
vec = envfit(allYrs_pca, na.omit(allYrs[, c(75:87)]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "blue", cex = 1) #plots vectors (environmetal variables)
pc1_scores <- as.data.frame(allYrs_pca$scores[, 1]*-1) # reverse signs to help with interpretation later
pc2_scores <- as.data.frame(allYrs_pca$scores[, 2]) # reverse signs to help with interpretation later
colnames(pc1_scores) <- c("PC_axis1")
colnames(pc2_scores) <- c("PC_axis2")
pc1_scores$ID <- rownames(pc1_scores)
pc2_scores$ID <- rownames(pc2_scores)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, pc1_scores, by = "ID")
allYrs <- dplyr::left_join(allYrs, pc2_scores, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

#choose which country for this task
TaskID <- read.csv("Data/LT_DriverTrends_TaskIDs.csv", as.is = T)

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "6"))
myCountry <- TaskID$country[which(TaskID$TaskID==task.id)]
allYrs <- subset(allYrs,country==myCountry)

#choose which response for this task
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]
allYrs$Response <- allYrs[,myResponse]
hist(allYrs$Response)

### Check distributions of response variables
# Main responses
hist(allYrs$flow)
summary(allYrs$flow)
hist(log10(allYrs$flow))

hist(allYrs$sflow)
summary(allYrs$sflow)
hist(log10(allYrs$sflow + 0.57))

hist(allYrs$temp)
hist(allYrs$stemp)

hist(allYrs$PC_axis1)
summary(allYrs$PC_axis1)
hist(log10(allYrs$PC_axis1 + 3.2))

hist(allYrs$PC_axis2)

# # transform variables that are right-skewed
# #third run transformation
# if(myResponse %in% c("flow")){
# 
#   allYrs$Response <- log10(allYrs$Response)
# 
# }else if(myResponse %in% "sflow"){
# 
#   allYrs$Response <- log10(allYrs$Response + 0.57)
#   
# }else if(myResponse %in% "PC_axis1"){
#   
#   allYrs$Response <- sqrt(allYrs$Response + 3.2)
# 
# }

hist(allYrs$Response)

# order by site site year
# allYrs <- allYrs[order(allYrs$year_wMissing),]

#### two-stage models ####
### fitting gls #####
library(nlme)
library(mgcv)

# write a function to consider day of year if sampling is more than 30 days apart
fitGLSModel_explan <- function(my_data){
  
  #centre Year - helps model convergence to center variables for the model
  my_data$cYear <- my_data$year_wMissing - median(my_data$year_wMissing)
  
  #or just have as an index starting from 1
  my_data$iYear <- my_data$year_wMissing - min(my_data$year_wMissing)+1
  
  #scale day of year
  my_data$cday_of_year <- (my_data$day_of_year - mean(my_data$day_of_year))/sd(my_data$day_of_year)
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(my_data$day_of_year)-min(my_data$day_of_year)
  
  if(length(my_data$Response[is.na(my_data$Response)])>0) {
    
    modelFits <- data.frame(estimate = NA,
                            se = NA,
                            pval = NA)
    
  }   else {
  
  if(maxDiffDays < 30) {
    myformula <- gls(Response ~ cYear, correlation = corAR1(form =~ iYear), data = my_data)
  } else{
    myformula <- gls(Response ~ cday_of_year + cYear, correlation = corAR1(form =~ iYear), data = my_data)
  }
  
  #fit model with gls
  gls_model <- myformula
  
  #extract model fits
  modelSummary <- summary(gls_model)$tTable["cYear", c("Value", "Std.Error", "p-value")]
  modelFits <- data.frame(estimate = modelSummary["Value"],
                          se = modelSummary["Std.Error"],
                          pval = modelSummary["p-value"])
}
  return(modelFits)
}


#the model is called in the function below
#apply function to an example dataset
est <- fitGLSModel_explan(allYrs[which(allYrs$site_id=="LTR1"),])
est
summary(est)

#loop for all sites
allsites <- sort(unique(allYrs$site_id))

trends <- lapply(allsites, function(x){
  fitGLSModel_explan(subset(allYrs, site_id == x))
})

trends <- data.frame(do.call(rbind, trends))
trends$siteID <- allsites
rownames(trends) <- 1:41

saveRDS(trends, file=paste0("Outputs/Driver_trends/trends__",myResponse,"__",myCountry,".RDS"))
#write.csv(trends, file=paste0("outputs/trends__",myResponse,"__",myCountry,".csv"))

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
