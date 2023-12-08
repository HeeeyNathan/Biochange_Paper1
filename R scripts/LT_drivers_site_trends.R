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
allYrs_pca = princomp(na.omit(allYrs[, c("salkalinity", "sEC", 
                                         "sNO3.N", "sNO2.N", "smineral.N", "sTot.N",
                                         "sPO4.P", "sTot.P")]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
# Define a pastel color palette
color_palette <- c("#84bcda", "#ecc30b", "#f37748", "#283845", "#b6244f")
# Define the two point shapes you want to use (adjust as needed)
shape_palette <- c(19, 17)
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", cex.axis = 1.5,
     main = "Principal coordinate analysis of collinear environmental data", xlab = "PC1", ylab = "PC 2", cex.main = 2, cex.lab = 1.5) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
points(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], col = color_palette[allYrs$ftype], pch = shape_palette[allYrs$fmodified], cex = 1.5)  # You can customize col, pch, and cex
vec = envfit(allYrs_pca, na.omit(allYrs[, c("salkalinity", "sEC", 
                                            "sNO3.N", "sNO2.N", "smineral.N", "sTot.N",
                                            "sPO4.P", "sTot.P")]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "#ECB246", cex = 1.5, font = 2, pos = 4) #plots vectors (environmetal variables)

# Check broken stick model
biplot(allYrs_pca)
allYrs_pca_scores <- (scores(allYrs_pca)[,1])
corr1 <- cor(na.omit(allYrs[, c(80, 83, 86:91)]), scores(allYrs_pca)[,1]) # correlation between original variables and principal components
round(corr1, 3)
corr2 <- cor(na.omit(allYrs[, c(80, 83, 86:91)]), scores(allYrs_pca)[,2]) # correlation between original variables and principal components
round(corr2, 3)
screeplot(allYrs_pca, bstick = TRUE, npcs = length(allYrs_pca$sdev))
(ev <- allYrs_pca$sdev^2)
n <- length (ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))}
bsm$p <- 100*bsm$p/n
bsm
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, main="Broken stick model", col=c("#95ccba",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch=15, col=c("#95ccba",2), bty="n")

# isolate pc axis for use as covariate in model
pc1_scores <- as.data.frame(allYrs_pca$scores[, 1]) # reverse signs to help with interpretation later
colnames(pc1_scores) <- c("PC_axis1")
pc1_scores$ID <- rownames(pc1_scores)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, pc1_scores, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# # ### remove missing covariate data ####
# allYrs_complete <- allYrs[complete.cases(allYrs[, c(58:72, 92)]),]

# choose which country for this task
TaskID <- read.csv("Data/LT_DriverTrends_TaskIDs.csv", as.is = T)

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "14"))
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
hist(log10(allYrs$flow+0.1))

hist(allYrs$temp)
hist(allYrs$pH)
hist(allYrs$o2_dis)

hist(allYrs$NH4.N)
summary(allYrs$NH4.N)
hist(log10(allYrs$NH4.N + 0.02))

hist(allYrs$PC_axis1)
summary(allYrs$PC_axis1)
hist(log10(allYrs$PC_axis1 + 3))

hist(allYrs$alkalinity)

hist(allYrs$EC)
hist(log10(allYrs$EC + 0.1))

hist(allYrs$NO3.N)
hist(log10(allYrs$NO3.N + 0.1))

hist(allYrs$NO2.N)
hist(log10(allYrs$NO2.N + 0.1))

hist(allYrs$mineral.N)
hist(log10(allYrs$mineral.N + 0.1))

hist(allYrs$Tot.N)
hist(log10(allYrs$Tot.N + 0.1))

hist(allYrs$PO4.P)
hist(log10(allYrs$PO4.P + 0.1))

hist(allYrs$Tot.P)
hist(log10(allYrs$Tot.P + 0.1))

# transform variables that are right-skewed
#third run transformation
if(myResponse %in% c("flow", "EC", "NO3.N", "NO2.N", "mineral.N", "Tot.N", "PO4.P", "Tot.P")){

  allYrs$Response <- log10(allYrs$Response + 0.1)

}else if(myResponse %in% "NH4.N"){

  allYrs$Response <- log10(allYrs$Response + 0.02)

}else if(myResponse %in% "PC_axis1"){
  
  allYrs$Response <- log10(allYrs$Response + 3)
  
}

hist(allYrs$Response)

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
trends

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
