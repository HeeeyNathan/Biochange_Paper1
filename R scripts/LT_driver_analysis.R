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
allYrs_pca = princomp(na.omit(allYrs[, c(82, 85, 88:93)]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
allYrs_pca
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", xlab = "Axis.1", ylab = "Axis.2", las = 1) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
points(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], col = "red", cex = 0.6, pch = 3) #adds text for samples
vec = envfit(allYrs_pca, na.omit(allYrs[, c(82, 85, 88:93)]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "blue", cex = 1) #plots vectors (environmetal variables)

# Check broken stick model
biplot(allYrs_pca)
summary(allYrs_pca)
allYrs_pca_scores <- (scores(allYrs_pca)[,1])
corr1 <- cor(na.omit(allYrs[, c(82, 85, 88:93)]), scores(allYrs_pca)[,1]) # correlation between original variables and principal components
round(corr1, 3)
corr2 <- cor(na.omit(allYrs[, c(82, 85, 88:93)]), scores(allYrs_pca)[,2]) # correlation between original variables and principal components
round(corr2, 3)
screeplot(allYrs_pca, bstick = TRUE, npcs = length(allYrs_pca$sdev))
(ev <- allYrs_pca$sdev^2)
n <- length (ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))}
bsm$p <- 100*bsm$p/n
bsm
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, main="Broken stick model", col=c("blue",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch=15, col=c("blue",2), bty="n")

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
#                  "diptera_spp_richness", "diptera_abundance", 
                  "crustacea_spp_richness", "crustacea_abundance", 
                  "insect_spp_richness", "insect_abundance", 
                  "mollusc_spp_richness", "mollusc_abundance", 
                  "annelid_spp_richness", "annelid_abundance")
Predictors_num <- c("sflow", "spH", "stemp", "so2_dis", "sNH4.N", "PC_axis1")
# Predictors_num <- c("flow", "temp", "sus_solid", "o2_dis", "pH", "BOD7", "NH4.N", "Nut.PC1")
# Predictors <- c("site_id", "year", "flow", "temp", "sus_solid", "o2_dis", "pH", "BOD7", "NH4.N", "Nut.PC1", "fmodified", "ftype")
# Identifiers <- c("country", "study_id", "site_id", "year", "site_code", "season", "sample_id", "day", "month", "date", "day_of_year", "latitude", "longitude", "cYear", "iYear", "fYear", "syear", "RivTypCols", "ModifiedCols")

VarsForFigs <- rep(1:4, times=c(6, 8, 6, 4))
VarToFit <- c("abundance", "spp_richness", "shannonsH", "E10", "turnover", "spp_rich_rare", 
              "FRic", "FEve", "FDis", "FRed", "F_turnover", 
              "FRic.SES", "FEve.SES", "FDis.SES", 
              "ept_spp_richness", "ept_abundance", 
#              "diptera_spp_richness", "diptera_abundance", 
              "crustacea_spp_richness", "crustacea_abundance", 
              "insect_spp_richness", "insect_abundance", 
              "mollusc_spp_richness", "mollusc_abundance", 
              "annelid_spp_richness", "annelid_abundance")

names(VarToFit) <- c("Total\nabundance", "Taxon\nrichness", "Shannon\ndiversity", "Shannon\nevenness", "Taxon\nturnover", "Rarified\ntaxon richness",
                     "Func.\nrichness", "Func.\nevevness", "Func.\ndispersion", "Func.\nredundancy", "Func.\ntrunover", 
                     "Standardised\nfunc. richness", "Standardised\nfunc. evenness", "Standardised\nfunc. dispersion", 
                     "EPT\nrichness", "EPT\nabundance",
#                     "Diptera\nrichness", "Diptera\nabundance",
                     "Crustacea\nrichness", "Crustacea\nabundance", 
                     "Insect\nrichness", "Insect\nabundance",
                     "Mollusc\nrichness", "Mollusc\nabundance",
                     "Annelid\nrichness", "Annelid\nabundance")

# FixedEffects <- c("sflow", "stemp", "ssus_solid", "so2_dis", "spH", "sBOD7", "sNH4.N", "Nut.PC1")
FixedEffects <- c("sflow", "spH", "stemp", "so2_dis", "sNH4.N", "PC_axis1")
# names(FixedEffects) <- c("Flow [m3/s]", "Temp. [°C]", "Suspended solids", "Dissolved oxygen [mg/L]", "pH", "Biol. oxygen demand [mg/L]", "Ammonium [mg/L]", "Nutrients") 
names(FixedEffects) <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA") 

# Plot data: covariates and responses
## Plot the precipitation covariates
pairs(allYrs[,Predictors_num], col = allYrs$RivTypCols)
pairs(allYrs[,Predictors_num], col = allYrs$ModifiedCols)
pairs(allYrs[,Predictors_num], col = allYrs$EQCCols)
pairs(allYrs[,Predictors_num], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity

## Plot the Responses
# coloured by river type
PlotResp <- function(respname, data, Cols) {
  plot(data$year, data[,respname], type="n", yaxt="n", xaxt="n", main=respname, 
       xlab="", ylab="")
  by(data, list(data$ftype, data$site_id), function(dat, resp, cols) {
    lines(dat$year, dat[,resp], col=cols[dat$ftype[1]])
  }, resp=respname, cols=Cols)
}

par(mfrow=c(6,4), mar=c(0.5,2,2,0.2), oma=c(2,2,0,0))
sapply(VarToFit, PlotResp, data=allYrs, Cols=RivTypCols)
mtext("Year", 1, outer=TRUE)
mtext("Response", 2, outer=TRUE)

# coloured by modified or not
PlotResp <- function(respname, data, Cols) {
  plot(data$year, data[,respname], type="n", yaxt="n", xaxt="n", main=respname, 
       xlab="", ylab="")
  by(data, list(data$fmodified, data$site_id), function(dat, resp, cols) {
    lines(dat$year, dat[,resp], col=cols[dat$fmodified[1]])
  }, resp=respname, cols=Cols)
}

par(mfrow=c(6,4), mar=c(0.5,2,2,0.2), oma=c(2,2,0,0))
sapply(VarToFit, PlotResp, data=allYrs, Cols=ModifiedCols)
mtext("Year", 1, outer=TRUE)
mtext("Response", 2, outer=TRUE)

# coloured by modified or not
PlotResp <- function(respname, data, Cols) {
  plot(data$year, data[,respname], type="n", yaxt="n", xaxt="n", main=respname, 
       xlab="", ylab="")
  by(data, list(data$fEQC, data$site_id), function(dat, resp, cols) {
    lines(dat$year, dat[,resp], col=cols[dat$fEQC[1]])
  }, resp=respname, cols=Cols)
}

par(mfrow=c(6,4), mar=c(0.5,2,2,0.2), oma=c(2,2,0,0))
sapply(VarToFit, PlotResp, data=allYrs, Cols=EQCCols)
mtext("Year", 1, outer=TRUE)
mtext("Response", 2, outer=TRUE)

# site trends, not coloured
PlotResp <- function(respname, data, Cols) {
  plot(data$year, data[,respname], type="n", yaxt="n", xaxt="n", main=respname, 
       xlab="", ylab="")
  by(data, list(data$site_id, data$site_id), function(dat, resp, cols) {
    lines(dat$year, dat[,resp], col="black")
  }, resp=respname, cols=Cols)
}

par(mfrow=c(6,4), mar=c(0.5,2,2,0.2), oma=c(2,2,0,0))
sapply(VarToFit, PlotResp, data=allYrs, Cols="black")
mtext("Year", 1, outer=TRUE)
mtext("Response", 2, outer=TRUE)

pairs(allYrs[,VarToFit], col = allYrs$RivTypCols)
pairs(allYrs[,VarToFit], col = allYrs$ModifiedCols)
pairs(allYrs[,VarToFit], col = allYrs$EQCCols)
pairs(allYrs[,VarToFit], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity

# Corrs <- cor(allYrs[,AllResponses])
# BigCorrs <- which(Corrs>0.9 & Corrs<1, arr.ind = TRUE)
# pairs(allYrs[,unique(rownames(BigCorrs))])
# pairs(allYrs[,unique(rownames(BigCorrs))], col = allYrs$RivTypCols)
# pairs(allYrs[,unique(rownames(BigCorrs))], col = allYrs$ModifiedCols)
# pairs(allYrs[,unique(rownames(BigCorrs))], col = allYrs$EQCCols)
# pairs(allYrs[,unique(rownames(BigCorrs))], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix") # Check env data for collinearity

### remove missing covariate data ####
allYrs_complete <- allYrs[complete.cases(allYrs[, c(79:93, 97)]),]

### Fitting models with lme4 ####
### With river type factor ####
VarToFit
FixedEffects
form.fixedS_type <- paste("(", paste(FixedEffects, collapse = " + "), ")*ftype", collapse = "")
form.randomlme <- "(1|site_id) + (1|year)" # (1|year) according to Weiss et al. (2023) & Daskalova et al. (2021)

## Make functions
FitModelContrasts <- function(resp, fixed, random, data) {
  require(lme4)
  if(!grepl("ftype", fixed)) stop("No river type variable in fixed effect")
  Levels <- levels(data$ftype)
  form <- formula(paste0(resp, " ~ ", fixed, "+", random))
  Summs <- lapply(Levels, function(lvl, dat, ff) {
    dat$ftype <- relevel(dat$ftype, ref=lvl)
    mod <- lmer(ff, data=dat)
    summ <- summary(mod)$coefficients
    summ <- summ[!grepl("ftype", rownames(summ)),]
    CI <- confint(mod) # calculates confidence intervals (CI) for each model
    summ <- cbind(summ, CI[4:10,]) # binds the CI to each summary table
    summ
    # isSingular(mod, tol = 1e-4) # a logical test to determine if the fitted mixed model is (almost/near) singular
  }, dat=data, ff=form)
  names(Summs) <- Levels
  Summs
}

Models.lme1 <- lapply(VarToFit, FitModelContrasts, fixed=form.fixedS_type, 
                      random=form.randomlme, data=allYrs_complete)

# Models.lme1$`Total\nabundance`$type1[, 4]

# # extract variables & generate functions
# GetEsts <- function(mod, var) {
#   require(plyr)
#   GetSumm <- function(summ, vr) return(summ[vr,])
#   ests <- ldply(.data=mod, .fun=GetSumm, vr=var)
#   return(ests)
# }
# 
# # (Ests <- GetEsts(mod=Models.lme1[[1]], var="sflow"))
# 
# PlotEffects <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
#   summ <- data.frame(summs[[name]])
#   if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
#   summ$LLim <- summ[, 5]
#   summ$ULim <- summ[, 6]
#   # summ$LLim <- summ$Estimate - summ$Std..Error
#   # summ$ULim <- summ$Estimate + summ$Std..Error
#   At.Y <- 1:nrow(summ)
#   plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
#        xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"))
#   segments(summ$LLim, At.Y, summ$ULim, At.Y, col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"))
#   abline(v=0, lty=3)
#   if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
#   if(title) title(main=name)
# }
# 
# PlotEsts <- function(var, models, nrows=6) {
#   Estimates <- lapply(Models.lme1, GetEsts, var=var)
#   NEsts <-  length(Estimates)
#   if(NEsts%%nrows!=0) warning("Plot not nice and square")
#   par(mfcol=c(nrows,ceiling(NEsts/nrows)), mar=c(2,1,2,0), oma=c(2,3,0,0))
#   lapply(names(Estimates)[1:nrows], PlotEffects, summs=Estimates, label.y=TRUE)
#   lapply(names(Estimates)[(nrows+1):NEsts], PlotEffects, 
#          summs=Estimates, label.y=FALSE)
#   mtext("Estimated Effect", 1, outer=TRUE)
#   mtext("River type", 2, outer=TRUE)
# }
# 
# # Flow effects
# PlotEsts(var="sflow", models=Models.lme1)
# 
# # Temperature effects
# PlotEsts(var="stemp", models=Models.lme1)
# 
# # # Suspended solids effects
# # PlotEsts(var="ssus_solid", models=Models.lme1)
# # # Dissolved oxygen effects
# # PlotEsts(var="so2_dis", models=Models.lme1)
# # # pH effects
# # PlotEsts(var="spH", models=Models.lme1)
# # # Biological Oxygen Demand effects
# # PlotEsts(var="sBOD7", models=Models.lme1)
# # # Ammonium effects
# # PlotEsts(var="sNH4.N", models=Models.lme1)
# 
# # PC axis 1 effects
# PlotEsts(var="PC_axis1", models=Models.lme1)
# 
# # PC axis 3 effects
# PlotEsts(var="PC_axis2", models=Models.lme1)

## Effects Plotted By variable
# The effects are now plotted by variable, with different countries in dfferen columns. 
# This should give an idea about the realtive magnitudes. Note that the covariates have been standardised, 
# so they represent the effects of changing the covariate by one standard deviation.
# PlotEffects <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
#   summ <- data.frame(summs[[name]])
#   if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
#   # Set the rownames to the specified values
#   new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
#   rownames(summ) <- new_rownames
#   summ$LLim <- summ[, 4]
#   summ$ULim <- summ[, 5]
#   At.Y <- 1:nrow(summ)
#   plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
#        xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), cex = 2)
#   segments(summ$LLim, At.Y, summ$ULim, At.Y, col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), lwd = 2)
#   abline(v=0, lty=3)
#   if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
#   if(title) title(main=name)
# }

PlotEffects <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
  summ <- data.frame(summs[[name]])
  if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
  # Set the rownames to the specified values
  new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
  rownames(summ) <- new_rownames
  summ$LLim <- summ[, 4]
  summ$ULim <- summ[, 5]
  At.Y <- 1:nrow(summ)
  # # Set line end style to square globally
  # par(lend = 2)
  plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
       xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), cex = 2)
  # Plot segments to the left of zero in "#f2cc84" with square ends
  left_of_zero <- summ$LLim < 0
  segments(summ$LLim[left_of_zero], At.Y[left_of_zero], pmin(0, summ$ULim[left_of_zero]), At.Y[left_of_zero], col = "#f2cc84", lwd = 2, lend = 2)  # Use lend = 2 for square ends
  # Plot segments to the right of zero in "#95ccba" with square ends
  right_of_zero <- summ$ULim > 0
  segments(pmax(0, summ$LLim[right_of_zero]), At.Y[right_of_zero], summ$ULim[right_of_zero], At.Y[right_of_zero], col = "#95ccba", lwd = 2, lend = 2)  # Use lend = 2 for square ends
  abline(v=0, lty=3)
  if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
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

# all taxon indices
VarToPlot <- (1:6)
par(mfrow=c(length(VarToPlot),5), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# select taxon indices
tiff(filename = "Plots/LT_Driver_Est_RivTyp_TaxoIndices.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(1:2, 4, 3, 5)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# all func. indices
VarToPlot <- 6 + (1:8)
par(mfrow=c(length(VarToPlot),5), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# select func. indices
tiff(filename = "Plots/LT_Driver_Est_RivTyp_FuncIndices.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- 6 + c(4, 1:3, 5)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# mixed taxon groups richness and abundance pt. 1
VarToPlot <- 14 + (1:6)
par(mfrow=c(length(VarToPlot),5), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon richness
tiff(filename = "Plots/LT_Driver_Est_RivTyp_TaxoGroupsRich.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(15, 19, 17, 21, 23)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# mixed taxon groups richness and abundance pt. 2
VarToPlot <- 20 + (1:4)
par(mfrow=c(length(VarToPlot),5), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon abundances
tiff(filename = "Plots/LT_Driver_Est_RivTyp_TaxoGroupAbund.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(16, 20, 18, 22, 24)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# extra plots
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),5), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

tiff(filename = "Plots/LT_Driver_Est_RivTyp_Extra.tiff", width = 10, height = 8, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),5), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme1)[VarToPlot], PlotRiverTypes, mod=Models.lme1)
mtext("River Type", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

### With heavily modified factor ####
VarToFit
FixedEffects
form.fixedS_mod <- paste("(", paste(FixedEffects, collapse = " + "), ")*fmodified", collapse = "")

## Make functions
FitModelContrasts2 <- function(resp, fixed, random, data) {
  require(lme4)
  if(!grepl("fmodified", fixed)) stop("No modification variable in fixed effect")
  Levels <- levels(data$fmodified)
  form <- formula(paste0(resp, " ~ ", fixed, "+", random))
  Summs <- lapply(Levels, function(lvl, dat, ff) {
    dat$fmodified <- relevel(dat$fmodified, ref=lvl)
    mod <- lmer(ff, data=dat)
    summ <- summary(mod)$coefficients
    summ <- summ[!grepl("fmodified", rownames(summ)),]
    CI <- confint(mod) # calculates confidence intervals (CI) for each model
    summ <- cbind(summ, CI[4:10,]) # binds the CI to each summary table
    summ
    # isSingular(mod, tol = 1e-4) # a logical test to determine if the fitted mixed model is (almost/near) singular
  }, dat=data, ff=form)
  names(Summs) <- Levels
  Summs
}

Models.lme2 <- lapply(VarToFit, FitModelContrasts2, fixed=form.fixedS_mod, 
                     random=form.randomlme, data=allYrs_complete)

# # Models.lme2$`Total\nabundance`$Yes[, 4]
# 
# # extract variables & generate functions
# GetEsts2 <- function(mod, var) {
#   require(plyr)
#   GetSumm <- function(summ, vr) return(summ[vr,])
#   ests <- ldply(.data=mod, .fun=GetSumm, vr=var)
#   return(ests)
# }
# 
# # (Ests <- GetEsts(mod=Models.lme2[[1]], var="sflow"))
# 
# PlotEffects2 <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
#   summ <- data.frame(summs[[name]])
#   if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
#   new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
#   rownames(summ) <- new_rownames
#   summ$LLim <- summ[, 5]
#   summ$ULim <- summ[, 6]
#   # summ$LLim <- summ$Estimate - summ$Std..Error
#   # summ$ULim <- summ$Estimate + summ$Std..Error
#   At.Y <- 1:nrow(summ)
#   plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
#        xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"))
#   segments(summ$LLim, At.Y, summ$ULim, At.Y, col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"))
#   abline(v=0, lty=3)
#   if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
#   if(title) title(main=name)
# }
# 
# PlotEsts2 <- function(var, models, nrows=6) {
#   Estimates <- lapply(Models.lme2, GetEsts2, var=var)
#   NEsts <-  length(Estimates)
#   if(NEsts%%nrows!=0) warning("Plot not nice and square")
#   par(mfcol=c(nrows,ceiling(NEsts/nrows)), mar=c(2,1,2,0), oma=c(2,3,0,0))
#   lapply(names(Estimates)[1:nrows], PlotEffects2, summs=Estimates, label.y=TRUE)
#   lapply(names(Estimates)[(nrows+1):NEsts], PlotEffects2, 
#          summs=Estimates, label.y=FALSE)
#   mtext("Estimated Effect", 1, outer=TRUE)
#   mtext("Modified", 2, outer=TRUE)
# }
# 
# # Flow effects
# PlotEsts2(var="sflow", models=Models.lme2)
# 
# # Temperature effects
# PlotEsts2(var="stemp", models=Models.lme2)
# 
# # # Suspended solids effects
# # PlotEsts2(var="ssus_solid", models=Models.lme2)
# # # Dissolved oxygen effects
# # PlotEsts2(var="so2_dis", models=Models.lme2)
# # # pH effects
# # PlotEsts2(var="spH", models=Models.lme2)
# # # Biological Oxygen Demand effects
# # PlotEsts2(var="sBOD7", models=Models.lme2)
# # # Ammonium effects
# # PlotEsts2(var="sNH4.N", models=Models.lme2)
# 
# # PC_axis1 effects
# PlotEsts2(var="PC_axis1", models=Models.lme2)
# 
# # PC_axis2 effects
# PlotEsts2(var="PC_axis2", models=Models.lme2)
 
## Effects Plotted By variable
# The effects are now plotted by variable, with different countries in dfferen columns. 
# This should give an idea about the realtive magnitudes. Note that the covariates have been standardised, 
# so they represent the effects of changing the covariate by one standard deviation.
# PlotEffects2 <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
#   summ <- data.frame(summs[[name]])
#   if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
#   # Set the rownames to the specified values
#   new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
#   rownames(summ) <- new_rownames
#   summ$LLim <- summ[, 4]
#   summ$ULim <- summ[, 5]
#   # summ$LLim <- summ$Estimate - summ$Std..Error
#   # summ$ULim <- summ$Estimate + summ$Std..Error
#   At.Y <- 1:nrow(summ)
#   plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
#        xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), cex = 2)
#   segments(summ$LLim, At.Y, summ$ULim, At.Y, col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), lwd = 2)
#   abline(v=0, lty=3)
#   if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
#   if(title) title(main=name)
# }

PlotEffects2 <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
  summ <- data.frame(summs[[name]])
  if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
  # Set the rownames to the specified values
  new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
  rownames(summ) <- new_rownames
  summ$LLim <- summ[, 4]
  summ$ULim <- summ[, 5]
  At.Y <- 1:nrow(summ)
  # # Set line end style to square globally
  # par(lend = 2)
  plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
       xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), cex = 2)
  # Plot segments to the left of zero in "#f2cc84" with square ends
  left_of_zero <- summ$LLim < 0
  segments(summ$LLim[left_of_zero], At.Y[left_of_zero], pmin(0, summ$ULim[left_of_zero]), At.Y[left_of_zero], col = "#f2cc84", lwd = 2, lend = 2)  # Use lend = 2 for square ends
  # Plot segments to the right of zero in "#95ccba" with square ends
  right_of_zero <- summ$ULim > 0
  segments(pmax(0, summ$LLim[right_of_zero]), At.Y[right_of_zero], summ$ULim[right_of_zero], At.Y[right_of_zero], col = "#95ccba", lwd = 2, lend = 2)  # Use lend = 2 for square ends
  abline(v=0, lty=3)
  if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
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

# all taxon indices
VarToPlot <- (1:6)
par(mfrow=c(length(VarToPlot),2), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# select taxon indices
tiff(filename = "Plots/LT_Driver_Est_Modified_TaxoIndices.tiff", width = 4, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(1:2, 4, 3, 5)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# all func. indices
VarToPlot <- 6 + (1:8)
par(mfrow=c(length(VarToPlot),2), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# select func. indices
tiff(filename = "Plots/LT_Driver_Est_Modified_FuncIndices.tiff", width = 4, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- 6 + c(4, 1:3, 5)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# mixed taxon groups richness and abundance pt. 1
VarToPlot <- 14 + (1:6)
par(mfrow=c(length(VarToPlot),2), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon richness
tiff(filename = "Plots/LT_Driver_Est_Modified_TaxoGroupsRich.tiff", width = 4, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(15, 19, 17, 21, 23)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# mixed taxon groups richness and abundance pt. 2
VarToPlot <- 20 + (1:4)
par(mfrow=c(length(VarToPlot),2), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon abundances
tiff(filename = "Plots/LT_Driver_Est_Modified_TaxoGroupsAbund.tiff", width = 4, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(16, 20, 18, 22, 24)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# extra plots
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),2), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

tiff(filename = "Plots/LT_Driver_Est_Modified_Extra.tiff", width = 4, height = 8, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),2), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme2)[VarToPlot], PlotModified, mod=Models.lme2)
mtext("Heavily Modified", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

### With EQC factor ####
### fitting the model
VarToFit
FixedEffects
(form.fixedS_EQC <- paste("(", paste(FixedEffects, collapse = " + "), ")*fEQC", collapse = ""))

## Make functions
FitModelContrasts3 <- function(resp, fixed, random, data) {
  require(lme4)
  if(!grepl("fEQC", fixed)) stop("No EQC variable in fixed effect")
  Levels <- levels(data$fEQC)
  form <- formula(paste0(resp, " ~ ", fixed, "+", random))
  Summs <- lapply(Levels, function(lvl, dat, ff) {
    dat$fEQC <- relevel(dat$fEQC, ref=lvl)
    mod <- lmer(ff, data=dat)
    summ <- summary(mod)$coefficients
    summ <- summ[!grepl("fEQC", rownames(summ)),]
    CI <- confint(mod) # calculates confidence intervals (CI) for each model
    summ <- cbind(summ, CI[4:10,]) # binds the CI to each summary table
    summ
    # isSingular(mod, tol = 1e-4) # a logical test to determine if the fitted mixed model is (almost/near) singular
  }, dat=data, ff=form)
  names(Summs) <- Levels
  Summs
}

Models.lme3 <- lapply(VarToFit, FitModelContrasts3, fixed=form.fixedS_EQC, 
                      random=form.randomlme, data=allYrs_complete)

# # Models.lme3$`Total\nabundance`$Good[, 4]
# 
# # extract variables & generate functions
# GetEsts3 <- function(mod, var) {
#   require(plyr)
#   GetSumm <- function(summ, vr) return(summ[vr,])
#   ests <- ldply(.data=mod, .fun=GetSumm, vr=var)
#   return(ests)
# }
# 
# # (Ests <- GetEsts3(mod=Models.lme3[[1]], var="sflow"))
# 
# PlotEffects3 <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
#   summ <- data.frame(summs[[name]])
#   if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
#   summ$LLim <- summ[, 5]
#   summ$ULim <- summ[, 6]
#   # summ$LLim <- summ$Estimate - summ$Std..Error
#   # summ$ULim <- summ$Estimate + summ$Std..Error
#   At.Y <- 1:nrow(summ)
#   plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
#        xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"))
#   segments(summ$LLim, At.Y, summ$ULim, At.Y, col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"))
#   abline(v=0, lty=3)
#   if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
#   if(title) title(main=name)
# }
# 
# PlotEsts3 <- function(var, models, nrows=6) {
#   Estimates <- lapply(Models.lme3, GetEsts3, var=var)
#   NEsts <-  length(Estimates)
#   if(NEsts%%nrows!=0) warning("Plot not nice and square")
#   par(mfcol=c(nrows,ceiling(NEsts/nrows)), mar=c(2,1,2,0), oma=c(2,3,0,0))
#   lapply(names(Estimates)[1:nrows], PlotEffects3, summs=Estimates, label.y=TRUE)
#   lapply(names(Estimates)[(nrows+1):NEsts], PlotEffects3, 
#          summs=Estimates, label.y=FALSE)
#   mtext("Estimated Effect", 1, outer=TRUE)
#   mtext("Ecological Quality Class", 2, outer=TRUE)
# }
# 
# # Flow effects
# PlotEsts3(var="sflow", models=Models.lme3)
# 
# # Temperature effects
# PlotEsts3(var="stemp", models=Models.lme3)
# 
# # # Suspended solids effects
# # PlotEsts3(var="ssus_solid", models=Models.lme3)
# # # Dissolved oxygen effects
# # PlotEsts3(var="so2_dis", models=Models.lme3)
# # # pH effects
# # PlotEsts3(var="spH", models=Models.lme3)
# # # Biological Oxygen Demand effects
# # PlotEsts3(var="sBOD7", models=Models.lme3)
# # # Ammonium effects
# # PlotEsts3(var="sNH4.N", models=Models.lme3)
# 
# # PC_axis1 effects
# PlotEsts3(var="PC_axis1", models=Models.lme3)
# 
# # PC_axis2 effects
# PlotEsts3(var="PC_axis2", models=Models.lme3)

## Effects Plotted By variable
# The effects are now plotted by variable, with different countries in dfferen columns. 
# This should give an idea about the realtive magnitudes. Note that the covariates have been standardised, 
# so they represent the effects of changing the covariate by one standard deviation.
# PlotEffects3 <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
#   summ <- data.frame(summs[[name]])
#   if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
#   # Set the rownames to the specified values
#   new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
#   rownames(summ) <- new_rownames
#   summ$LLim <- summ[, 4]
#   summ$ULim <- summ[, 5]
#   # summ$LLim <- summ$Estimate - summ$Std..Error
#   # summ$ULim <- summ$Estimate + summ$Std..Error
#   At.Y <- 1:nrow(summ)
#   plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
#        xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), cex = 2)
#   segments(summ$LLim, At.Y, summ$ULim, At.Y, col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), lwd = 2)
#   abline(v=0, lty=3)
#   if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
#   if(title) title(main=name)
# }

PlotEffects3 <- function(name, summs, label.y=TRUE, title=TRUE, removeInt=FALSE) {
  summ <- data.frame(summs[[name]])
  if(removeInt) summ <- summ[rownames(summ)!="(Intercept)",]
  # Set the rownames to the specified values
  new_rownames <- c("Flow", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA")
  rownames(summ) <- new_rownames
  summ$LLim <- summ[, 4]
  summ$ULim <- summ[, 5]
  At.Y <- 1:nrow(summ)
  # # Set line end style to square globally
  # par(lend = 2)
  plot(summ$Estimate, At.Y, yaxt="n", ann=FALSE, 
       xlim=c(min(summ$LLim), max(summ$ULim)), ylim=c(0.5, nrow(summ)+0.5), col = ifelse(summ$Estimate >= 0, "#95ccba", "#f2cc84"), cex = 2)
  # Plot segments to the left of zero in "#f2cc84" with square ends
  left_of_zero <- summ$LLim < 0
  segments(summ$LLim[left_of_zero], At.Y[left_of_zero], pmin(0, summ$ULim[left_of_zero]), At.Y[left_of_zero], col = "#f2cc84", lwd = 2, lend = 2)  # Use lend = 2 for square ends
  # Plot segments to the right of zero in "#95ccba" with square ends
  right_of_zero <- summ$ULim > 0
  segments(pmax(0, summ$LLim[right_of_zero]), At.Y[right_of_zero], summ$ULim[right_of_zero], At.Y[right_of_zero], col = "#95ccba", lwd = 2, lend = 2)  # Use lend = 2 for square ends
  abline(v=0, lty=3)
  if(label.y) axis(2, at=At.Y, labels = rownames(summ), las=1)
  if(title) title(main=name)
}

PlotEQCs <- function(wh, mod) {
  AbundEsts <- mod[[wh]]
  names(AbundEsts)[1] <- "Bad/Poor"
  names(AbundEsts)[2] <- "Moderate"
  names(AbundEsts)[3] <- "Good"
  names(AbundEsts)[4] <- "High"
  # Define the desired order of plotting columns
  PlotEffects3(name=names(AbundEsts)[1], summs=AbundEsts, 
              label.y=TRUE, title=TRUE, removeInt = TRUE)
  sapply(names(AbundEsts)[-1], PlotEffects3, summs=AbundEsts, 
         label.y=FALSE, title=TRUE, removeInt = TRUE)
  mtext(wh, 4, outer=FALSE, line=2)
}

# all taxon indices
VarToPlot <- (1:6)
par(mfrow=c(length(VarToPlot), 4), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# select taxon indices
tiff(filename = "Plots/LT_Driver_Est_EQR_TaxoIndices.tiff", width = 8, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(1:2, 4, 3, 5)
par(mfrow=c(length(VarToPlot),4), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# all func. indices
VarToPlot <- 6 + (1:8)
par(mfrow=c(length(VarToPlot),4), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# select func. indices
tiff(filename = "Plots/LT_Driver_Est_EQR_FuncIndices.tiff", width = 8, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- 6 + c(4, 1:3, 5)
par(mfrow=c(length(VarToPlot),4), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# mixed taxon groups richness and abundance pt. 1
VarToPlot <- 14 + (1:6)
par(mfrow=c(length(VarToPlot),4), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon richness
tiff(filename = "Plots/LT_Driver_Est_EQR_TaxoGroupsRich.tiff", width = 8, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(15, 19, 17, 21, 23)
par(mfrow=c(length(VarToPlot),4), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# mixed taxon groups richness and abundance pt. 2
VarToPlot <- 20 + (1:4)
par(mfrow=c(length(VarToPlot),4), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

# taxon abundances
tiff(filename = "Plots/LT_Driver_Est_EQR_TaxoGroupsAbund.tiff", width = 8, height = 10, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(16, 20, 18, 22, 24)
par(mfrow=c(length(VarToPlot),4), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

# extra plots
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),4), mar=c(2,2,2,2), oma=c(2,4,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)

tiff(filename = "Plots/LT_Driver_Est_EQR_Extra.tiff", width = 8, height = 8, units = 'in', res = 600, compression = 'lzw')
VarToPlot <- c(6, 12, 13, 14)
par(mfrow=c(length(VarToPlot),4), mar=c(2,0.5,2,2), oma=c(2,6,2,2))
sapply(names(Models.lme3)[VarToPlot], PlotEQCs, mod=Models.lme3)
mtext("Ecological Quality Class", 3, outer=TRUE, font = 2)
mtext("Estimate", 1, outer = TRUE, line = 1)
dev.off()

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
