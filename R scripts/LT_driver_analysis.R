##### load functions, packages, and external sources #####
rm(list=ls())

library(lubridate)
library(pacman)
library(vegan)
library(INLA)
library(car) # for logit()
library(xtable) # for tables
library(RColorBrewer)

source("Additional functions/HighstatLibV10.R")

##### read in data & format #####
d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T)
allYrs <- d1[!is.na(d1$site_id_wMissing),]

# Remove the year 2019 from all non-within-site analyses
allYrs <- allYrs[allYrs$year != "2019",]

# check histograms of response variables
hist(allYrs$abundance)
hist(allYrs$spp_richness)
hist(allYrs$E10)
hist(allYrs$shannonsH)
hist(allYrs$turnover)
hist(allYrs$spp_rich_rare)
hist(allYrs$FRic)
hist(allYrs$FEve)
hist(allYrs$FDis)
hist(allYrs$FRed)
hist(allYrs$F_turnover)
hist(allYrs$FRic.SES)
hist(allYrs$FEve.SES)
hist(allYrs$FDis.SES)
hist(allYrs$ept_spp_richness)
hist(allYrs$ept_abundance)
hist(allYrs$crustacea_spp_richness)
hist(allYrs$crustacea_abundance)
hist(allYrs$diptera_spp_richness)
hist(allYrs$diptera_abundance)
hist(allYrs$insect_spp_richness)
hist(allYrs$insect_abundance)
hist(allYrs$mollusc_spp_richness)
hist(allYrs$mollusc_abundance)
hist(allYrs$annelid_spp_richness)
hist(allYrs$annelid_abundance)

# transformations of response variables
allYrs$abundance <- log10(allYrs$abundance + 1)
allYrs$ept_abundance <- log10(allYrs$ept_abundance + 1)
allYrs$diptera_abundance <- log10(allYrs$diptera_abundance + 1)
allYrs$crustacea_abundance <- log10(allYrs$crustacea_abundance + 1)
allYrs$insect_abundance <- log10(allYrs$insect_abundance + 1)
allYrs$mollusc_abundance <- log10(allYrs$mollusc_abundance + 1)
allYrs$annelid_abundance <- log10(allYrs$annelid_abundance + 1)

allYrs$spp_richness <- sqrt(allYrs$spp_richness)
allYrs$spp_rich_rare <- sqrt(allYrs$spp_rich_rare)
allYrs$ept_spp_richness <- sqrt(allYrs$ept_spp_richness)
allYrs$diptera_spp_richness <- sqrt(allYrs$diptera_spp_richness)
allYrs$crustacea_spp_richness <- sqrt(allYrs$crustacea_spp_richness)
allYrs$insect_spp_richness <- sqrt(allYrs$insect_spp_richness)
allYrs$mollusc_spp_richness <- sqrt(allYrs$mollusc_spp_richness)
allYrs$annelid_spp_richness <- sqrt(allYrs$annelid_spp_richness)

allYrs$F_turnover <- log10(allYrs$F_turnover + 0.01)
allYrs$FRic <- log10(allYrs$FRic + 0.28)

allYrs$FEve <- allYrs$FEve^2

allYrs$FDis <- allYrs$FDis^3
allYrs$turnover <- allYrs$turnover^3

# transform preditor variables (if necessary)
hist(allYrs$flow)
hist(allYrs$temp)

allYrs$flow <- log10(allYrs$flow)

# transform year variable
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1
allYrs$fYear <- factor(allYrs$year_wMissing)

# Check histograms of response variables agaiin after transformation
hist(allYrs$abundance)
hist(allYrs$spp_richness)
hist(allYrs$E10)
hist(allYrs$shannonsH)
hist(allYrs$turnover)
hist(allYrs$spp_rich_rare)
hist(allYrs$FRic)
hist(allYrs$FEve)
hist(allYrs$FDis)
hist(allYrs$FRed)
hist(allYrs$F_turnover)
hist(allYrs$FRic.SES)
hist(allYrs$FEve.SES)
hist(allYrs$FDis.SES)
hist(allYrs$ept_spp_richness)
hist(allYrs$ept_abundance)
hist(allYrs$diptera_spp_richness)
hist(allYrs$diptera_abundance)
hist(allYrs$insect_spp_richness)
hist(allYrs$insect_abundance)
hist(allYrs$crustacea_spp_richness)
hist(allYrs$crustacea_abundance)
hist(allYrs$mollusc_spp_richness)
hist(allYrs$mollusc_abundance)
hist(allYrs$annelid_spp_richness)
hist(allYrs$annelid_abundance)

# Change river type variables
allYrs$river_type[allYrs$river_type == 1] <- "type1"
allYrs$river_type[allYrs$river_type == 2] <- "type2"
allYrs$river_type[allYrs$river_type == 3] <- "type3"
allYrs$river_type[allYrs$river_type == 4] <- "type4"
allYrs$river_type[allYrs$river_type == 5] <- "type5"

# change EQC variable names
allYrs$EQC[allYrs$EQC == "Poor"] <- "1Poor"
allYrs$EQC[allYrs$EQC == "Moderate"] <- "2Moderate"
allYrs$EQC[allYrs$EQC == "Good"] <- "3Good"
allYrs$EQC[allYrs$EQC == "Very good"] <- "3High"

# make factors
allYrs$fsite_id <- factor(allYrs$site_id)
allYrs$fYear <- factor(allYrs$year_wMissing)
allYrs$fmodified <- factor(allYrs$Heavily_modified)
allYrs$ftype <- factor(allYrs$river_type)
allYrs$fEQC <- factor(allYrs$EQC)

# Scale the covariates
# function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- lapply(df, function(x) if(is.numeric(x)){
    scale(x, center=TRUE, scale=TRUE)
  } else x)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df[,c(1, 2)], newd)
}

# apply function
sPredictors <- scaleVars(allYrs[, c(3:4, 58:72)])
sPredictors <- subset(sPredictors, select = -c(site_id, year, ssite_id, syear)) # remove ID variable
sPredictors$ID <- rownames(sPredictors)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, sPredictors, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# add colours for factors
RivTypCols <- brewer.pal(5, "Dark2")
ModifiedCols <- c("#1B9E77", "#D95F02")
EQCCols <- brewer.pal(4, "Dark2")
names(RivTypCols) <- levels(allYrs$ftype)
names(ModifiedCols) <- levels(allYrs$fmodified)
names(EQCCols) <- levels(allYrs$fEQC)

allYrs$RivTypCols <- RivTypCols[as.numeric(allYrs$ftype)]
allYrs$ModifiedCols <- ModifiedCols[as.numeric(allYrs$fmodified)]
allYrs$EQCCols <- EQCCols[as.numeric(allYrs$fEQC)]

colnames(allYrs)

# # Check variable colinnearity
pairs(allYrs[, c(82, 85, 88:93)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity

# Construct PCA to check the environmental variables and their relationships
allYrs_pca = princomp(na.omit(allYrs[, c(83, 86, 89:94)]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
allYrs_pca
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", xlab = "Axis.1", ylab = "Axis.2", las = 1) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
points(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], col = "red", cex = 0.6, pch = 3) #adds text for samples
vec = envfit(allYrs_pca, na.omit(allYrs[, c(83, 86, 89:94)]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "blue", cex = 1) #plots vectors (environmetal variables)

# Check broken stick model
screeplot(allYrs_pca, bstick = TRUE, npcs = length(allYrs_pca$sdev))

# isolate pc axis for use as covariate in model
pc1_scores <- as.data.frame(allYrs_pca$scores[, 1]) # reverse signs to help with interpretation later
colnames(pc1_scores) <- c("PC_axis1")
pc1_scores$ID <- rownames(pc1_scores)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, pc1_scores, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# Create subsets for responses vs predictors
AllResponses <- c("abundance", "spp_richness", "shannonsH", "E10", "turnover", "spp_rich_rare", 
                  "FRic", "FEve", "FDis", "FRed", "F_turnover",
                  "FRic.SES", "FEve.SES", "FDis.SES", 
                  "ept_spp_richness", "ept_abundance", 
                  "crustacea_spp_richness", "crustacea_abundance", 
                  "insect_spp_richness", "insect_abundance", 
                  "mollusc_spp_richness", "mollusc_abundance", 
                  "annelid_spp_richness", "annelid_abundance")
Predictors_num <- c("sflow", "spH", "stemp", "so2_dis", "sNH4.N", "PC_axis1")

VarsForFigs <- rep(1:4, times=c(6, 8, 6, 4))
VarToFit <- c("abundance", "spp_richness", "shannonsH", "E10", "turnover", "spp_rich_rare", 
              "FRic", "FEve", "FDis", "FRed", "F_turnover", 
              "FRic.SES", "FEve.SES", "FDis.SES", 
              "ept_spp_richness", "ept_abundance", 
              "crustacea_spp_richness", "crustacea_abundance", 
              "insect_spp_richness", "insect_abundance", 
              "mollusc_spp_richness", "mollusc_abundance", 
              "annelid_spp_richness", "annelid_abundance")

names(VarToFit) <- c("Total\nabundance", "Taxon\nrichness", "Shannon\ndiversity", "Shannon\nevenness", "Taxon\nturnover", "Rarified\ntaxon richness",
                     "Func.\nrichness", "Func.\nevevness", "Func.\ndispersion", "Func.\nredundancy", "Func.\ntrunover", 
                     "Standardised\nfunc. richness", "Standardised\nfunc. evenness", "Standardised\nfunc. dispersion", 
                     "EPT\nrichness", "EPT\nabundance",
                     "Crustacea\nrichness", "Crustacea\nabundance", 
                     "Insect\nrichness", "Insect\nabundance",
                     "Mollusc\nrichness", "Mollusc\nabundance",
                     "Annelid\nrichness", "Annelid\nabundance")

FixedEffects <- c("sflow", "spH", "stemp", "so2_dis", "sNH4.N", "PC_axis1")
names(FixedEffects) <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA") 

# Plot data: covariates and responses
## Plot the precipitation covariates
pairs(allYrs[,Predictors_num], col = allYrs$RivTypCols)
pairs(allYrs[,Predictors_num], col = allYrs$ModifiedCols)
pairs(allYrs[,Predictors_num], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity

### remove missing covariate data ####
allYrs_complete <- allYrs[complete.cases(allYrs[, c(80:94, 98)]),]

### Fitting models with lme4 ####
### With river type factor ####
VarToFit
FixedEffects
form.fixedS_type <- paste("(", paste(FixedEffects, collapse = " + "), ")*ftype", collapse = "")
form.randomlme <- "(1|fsite_id) + (1|fYear)" # (1|year) according to Weiss et al. (2023) & Daskalova et al. (2021)

## Make functions
FitModelContrasts <- function(resp, fixed, random, data) {
  require(lme4) # isSingular(mod, tol = 1e-4) # a logical test to determine if the fitted mixed model is (almost/near) singular
  if(!grepl("ftype", fixed)) stop("No river type variable in fixed effect")
  Levels <- levels(data$ftype)
  form <- formula(paste0(resp, " ~ ", fixed, "+", random))
  Summs <- lapply(Levels, function(lvl, dat, ff) {
    dat$ftype <- relevel(dat$ftype, ref=lvl)
    mod <- lmer(ff, data=dat)
    summ <- summary(mod)$coefficients
    summ <- summ[!grepl("ftype", rownames(summ)),]
    CI <- confint(mod) # calculates confidence intervals (CI) for each model
    # Match the rownames of summ with CI indices
    matchedCI <- CI[rownames(summ), ]
    # Bind the matched CIs to the summary
    summ <- cbind(summ, matchedCI)
    summ
  }, dat=data, ff=form)
  names(Summs) <- Levels
  Summs
}

Models.lme1 <- lapply(VarToFit, FitModelContrasts, fixed=form.fixedS_type, 
                      random=form.randomlme, data=allYrs_complete)

PlotEffects <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
  # Extract the summary information for the specified variable
  summ <- data.frame(summs[[name]])
  # If specified, remove the intercept term from the summary
  if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
  # Set custom row names for the summary dataframe
  new_rownames <- c("Discharge", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
  rownames(summ) <- new_rownames
  # Create new columns for lower and upper limits
  summ$LLim <- summ[, 4]
  summ$ULim <- summ[, 5]
  # Create a sequence for the y-axis values
  At.Y <- 1:nrow(summ)
  # Check if confidence intervals include the null hypothesis
  include_null <- (summ$LLim <= 0) & (summ$ULim >= 0)
  # Plot the effect estimates with different colors based on null hypothesis inclusion
  plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
       xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5),
       col = ifelse(include_null, "gray90", ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84")), pch = 16, cex = 2)
  # Plot segments to the left of zero in "#f2cc84" with square ends
  left_of_zero <- summ$LLim < 0
  segments(summ$LLim[left_of_zero], At.Y[left_of_zero], pmin(0, summ$ULim[left_of_zero]), At.Y[left_of_zero],
           col = ifelse(include_null[left_of_zero], "gray90", "#f2cc84"), lwd = 3, lend = 2)
  # Plot segments to the right of zero in "#95ccba" with square ends
  right_of_zero <- summ$ULim > 0
  segments(pmax(0, summ$LLim[right_of_zero]), At.Y[right_of_zero], summ$ULim[right_of_zero], At.Y[right_of_zero],
           col = ifelse(include_null[right_of_zero], "gray90", "#95ccba"), lwd = 3, lend = 2)
  # Add a vertical dashed line at zero
  abline(v=0, lty=3)
  # If specified, label the y-axis with custom names
  if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
  # If specified, add a title to the plot
  if(title) title(main=name)
}

PlotRiverTypes <- function(wh, mod) {
  AbundEsts <- mod[[wh]]
  names(AbundEsts)[1] <- "1"
  names(AbundEsts)[2] <- "2"
  names(AbundEsts)[3] <- "3"
  names(AbundEsts)[4] <- "4"
  names(AbundEsts)[5] <- "5"
  PlotEffects(name=names(AbundEsts)[1], summs=AbundEsts, 
              label.y=TRUE, title=TRUE, removeInt = TRUE)
  sapply(names(AbundEsts)[-1], PlotEffects, summs=AbundEsts, 
         label.y=FALSE, title=TRUE, removeInt = TRUE)
  mtext(wh, 4, outer=FALSE, line=2)
}

# Figure 4 - selected metrics by river type
# select taxon indices
tiff(filename = "Plots/Figure_6_RivType_Drivers.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(1:2, 4, 10, 7)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# all taxon indices
VarToPlot <- (1:6)
par(mfrow=c(length(VarToPlot),5), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# all func. indices
VarToPlot <- 6 + (1:5)
par(mfrow=c(length(VarToPlot),5), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon richness
VarToPlot <- c(15, 19, 17, 21, 23)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon abundances
VarToPlot <- c(16, 20, 18, 22, 24)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# extra plots
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

### With heavily modified factor ####
VarToFit
FixedEffects
form.fixedS_mod <- paste("(", paste(FixedEffects, collapse = " + "), ")*fmodified", collapse = "")

## Make functions
FitModelContrasts2 <- function(resp, fixed, random, data) {
  require(lme4) # isSingular(mod, tol = 1e-4) # a logical test to determine if the fitted mixed model is (almost/near) singular
  if(!grepl("fmodified", fixed)) stop("No modification variable in fixed effect")
  Levels <- levels(data$fmodified)
  form <- formula(paste0(resp, " ~ ", fixed, "+", random))
  Summs <- lapply(Levels, function(lvl, dat, ff) {
    dat$fmodified <- relevel(dat$fmodified, ref=lvl)
    mod <- lmer(ff, data=dat)
    summ <- summary(mod)$coefficients
    summ <- summ[!grepl("fmodified", rownames(summ)),]
    CI <- confint(mod) # calculates confidence intervals (CI) for each model
    # Match the rownames of summ with CI indices
    matchedCI <- CI[rownames(summ), ]
    # Bind the matched CIs to the summary
    summ <- cbind(summ, matchedCI)
    summ
  }, dat=data, ff=form)
  names(Summs) <- Levels
  Summs
}

Models.lme2 <- lapply(VarToFit, FitModelContrasts2, fixed=form.fixedS_mod, 
                     random=form.randomlme, data=allYrs_complete)


PlotEffects2 <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
  # Extract the summary information for the specified variable
  summ <- data.frame(summs[[name]])
  # If specified, remove the intercept term from the summary
  if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
  # Set custom row names for the summary dataframe
  new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
  rownames(summ) <- new_rownames
  # Create new columns for lower and upper limits
  summ$LLim <- summ[, 4]
  summ$ULim <- summ[, 5]
  # Create a sequence for the y-axis values
  At.Y <- 1:nrow(summ)
  # Check if confidence intervals include the null hypothesis
  include_null <- (summ$LLim <= 0) & (summ$ULim >= 0)
  # Plot the effect estimates with different colors based on null hypothesis inclusion
  plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
       xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5),
       col = ifelse(include_null, "gray90", ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84")), pch = 16, cex = 2)
  # Plot segments to the left of zero in "#f2cc84" with square ends
  left_of_zero <- summ$LLim < 0
  segments(summ$LLim[left_of_zero], At.Y[left_of_zero], pmin(0, summ$ULim[left_of_zero]), At.Y[left_of_zero],
           col = ifelse(include_null[left_of_zero], "gray90", "#f2cc84"), lwd = 3, lend = 2)
  # Plot segments to the right of zero in "#95ccba" with square ends
  right_of_zero <- summ$ULim > 0
  segments(pmax(0, summ$LLim[right_of_zero]), At.Y[right_of_zero], summ$ULim[right_of_zero], At.Y[right_of_zero],
           col = ifelse(include_null[right_of_zero], "gray90", "#95ccba"), lwd = 3, lend = 2)
  # Add a vertical dashed line at zero
  abline(v=0, lty=3)
  # If specified, label the y-axis with custom names
  if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
  # If specified, add a title to the plot
  if(title) title(main=name)
}

PlotModified <- function(wh, mod) {
  AbundEsts <- mod[[wh]]
  PlotEffects(name=names(AbundEsts)[1], summs=AbundEsts, 
              label.y=TRUE, title=TRUE, removeInt = TRUE)
  sapply(names(AbundEsts)[-1], PlotEffects2, summs=AbundEsts, 
         label.y=FALSE, title=TRUE, removeInt = TRUE)
  mtext(wh, 4, outer=FALSE, line=2)
}

# Selected metrics by modification
tiff(filename = "Plots/Figure_7_Modified_Drivers.tiff", width = 4, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(1:2, 4, 10, 7)
par(mfrow=c(length(VarToPlot),2), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# select taxon indices
VarToPlot <- c(1:2, 4, 3, 5)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# select func. indices
VarToPlot <- 6 + c(4, 1:3, 5)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon richness
VarToPlot <- c(15, 19, 17, 21, 23)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon abundances
VarToPlot <- c(16, 20, 18, 22, 24)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# extra plots
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

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
