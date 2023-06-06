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

# Remove the year 2019 from all non-within-site analyses
allYrs <- allYrs[allYrs$year != "2019",]

# tranform predictor variables (if necessary)
hist(allYrs$flow)
hist(allYrs$temp)

allYrs$flow <- log10(allYrs$flow)

# Change river type variables
allYrs$river_type[allYrs$river_type == 1] <- "type1"
allYrs$river_type[allYrs$river_type == 2] <- "type2"
allYrs$river_type[allYrs$river_type == 3] <- "type3"
allYrs$river_type[allYrs$river_type == 4] <- "type4"
allYrs$river_type[allYrs$river_type == 5] <- "type5"

# make factors
allYrs$fYear <- factor(allYrs$year_wMissing)
allYrs$fmodified <- as.factor(allYrs$Heavily_modified)
allYrs$ftype <- as.factor(allYrs$river_type)
allYrs$fEQC <- as.factor(allYrs$EQC)

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
sPredictors <- scaleVars(allYrs[, c(3:4, 58:72)])
sPredictors <- subset(sPredictors, select = -c(site_id, year, ssite_id, syear)) # remove ID variable
sPredictors$ID <- rownames(sPredictors)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, sPredictors, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# Construct PCA to check the environmental variables and their relationships
pairs(allYrs[, c(77:91)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity
pairs(allYrs[, c(79:91)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity
pairs(allYrs[, c(80, 83, 86:91)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity
allYrs_pca = princomp(na.omit(allYrs[, c(79:91)]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", xlab = "Axis.1", ylab = "Axis.2", las = 1) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$site_code, col = RivTypCols) #adds text for samples
text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$site_id, col = "red", cex = 0.6) #adds text for samples
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$ftype, col = "darkblue") #adds text for samples
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$fmodified, col = "forestgreen") #adds text for samples
vec = envfit(allYrs_pca, na.omit(allYrs[, c(79:91)]), choices = c(1, 2)) #computes the environmental variable vectors
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

response_lmer <- unique(allYrs[,c("site_code", "site_id", "year",
                                 "abundance", 
                                 "spp_richness", "shannonsH", "E10", "turnover", "spp_rich_rare",
                                 "ept_spp_richness", "diptera_spp_richness", "insect_spp_richness", "mollusc_spp_richness", "annelid_spp_richness",
                                 "ept_abundance", "diptera_abundance", "insect_abundance", "mollusc_abundance", "annelid_abundance",
                                 "FRic", "FEve", "FDis", "FRed", "F_turnover", 
                                 "FRic.SES", "FEve.SES", "FDis.SES", 
                                 "sflow", "stemp", "PC_axis1", "PC_axis2", 
                                 "fmodified", "ftype", "fEQC")])

### remove missing covariate data ####
response_lmer <- response_lmer[complete.cases(response_lmer[, c(28:31)]),]

### model time ####
### All sites, all years ####
library(lme4)

### get response for this task ##
TaskID <- read.csv("Data/LT_ResponseTrends_TaskIDs.csv", as.is = T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "24"))
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]

# choose which response for this task
response_lmer$Response <- response_lmer[,myResponse]
hist(response_lmer$Response)

# # transformation 
# if(myResponse %in% c("abundance", "ept_abundance", "diptera_abundance", "insect_abundance", "mollusc_abundance", "annelid_abundance")){
#   response_lmer$Response <- log10(response_lmer$Response + 1)
# }

# # second run transformation
# if(myResponse %in% c("abundance", "ept_abundance", "diptera_abundance", "insect_abundance", "mollusc_abundance", "annelid_abundance")){
#   
#   response_lmer$Response <- log10(response_lmer$Response + 1)
#   
# }else if(myResponse %in% c("F_turnover", "FRic")){
#   
#   response_lmer$Response <- log10(response_lmer$Response + 0.01)
#   
# }else if(myResponse %in% c("FDis", "turnover")){
#   
#   response_lmer$Response <- response_lmer$Response^2
#   
# }else if(myResponse %in% c("ept_spp_richness", "diptera_spp_richness", "insect_spp_richness", "mollusc_spp_richness", "annelid_spp_richness")){
#   
#   response_lmer$Response <-sqrt(response_lmer$Response)
# }

# fourth run transformation
if(myResponse %in% c("abundance", "ept_abundance", "diptera_abundance", 
                     "insect_abundance", "mollusc_abundance", 
                     "annelid_abundance", "crustacea_abundance")){
  
  response_lmer$Response <- log10(response_lmer$Response + 1)
  
}else if(myResponse %in% "F_turnover"){
  
  response_lmer$Response <- log10(response_lmer$Response + 0.01)
  
}else if(myResponse %in% "FRic") {
  
  response_lmer$Response <- log10(response_lmer$Response + 0.28)
  
}else if(myResponse %in% "FEve"){
  
  response_lmer$Response <- response_lmer$Response^2
  
}else if(myResponse %in% c("FDis", "turnover")){
  
  response_lmer$Response <- response_lmer$Response^3
  
}else if(myResponse %in% c("spp_richness", "spp_rich_rare", 
                           "insect_spp_richness","ept_spp_richness", 
                           "diptera_spp_richness", "mollusc_spp_richness", 
                           "annelid_spp_richness", "crustacea_spp_richness")){
  
  response_lmer$Response <- sqrt(response_lmer$Response)
  
}

# check response again
hist(response_lmer$Response)

### run model ###
fit1 <- lmer(Response ~ sflow + stemp + PC_axis1 + PC_axis2 + (1|site_id) + (1|year), data = response_lmer)

summ <- summary(fit1)$coefficients
summ <- data.frame(summ)
summ$covariate <- rownames(summ)
rownames(summ) <- (1:5)

CI <- confint(fit1)
summ <- cbind(summ, CI[4:8,])

# check residuals
residuals <- resid(fit1)
plot(predict(fit1), residuals)
abline(0, 0)
hist(residuals)

# #### save output ###
# 
saveRDS(summ, file=paste0("Outputs/Drivers_all_sites/drivers__",myResponse,"__Lithuania_new.RDS"))
# write.csv(summ,file=paste0("Outputs/Drivers_all_sites/drivers__",myResponse,"__Lithuania.csv"))

##### CLEAN UP --------------------
library(pacman)
# Clear data
# rm(list = ls())  # Removes all objects from environment
# Clear packages
# p_unload(all)  # Remove all contributed packages
# Clear plots
# graphics.off()  # Clears plots, closes all graphics devices
# Clear console
# cat("\014")  # Mimics ctrl+L
# Clear mind :)
