#### function to extract posterior distribution of the trends ####

getTrendProbability <- function(fit){
  mySamples <- posterior_samples(fit, pars="b_Intercept")
  data.frame(probIncrease = mean(mySamples>0),probDecrease = mean(mySamples<0))
}

### end of functions ############################################

### site-level trends ####

#read in model for each response and save fixed effects

library(tidyverse)
library(pacman)

response_gls <- readRDS("outputs/glsTrends_site_level_drivers.rds")

#pivot responses
response_gls_pivot <- response_gls %>%
  select(c(Response,estimate,site_id)) %>%
  pivot_wider(names_from = "Response",
              values_from = "estimate")

#get site metadata
d1 <- read.csv("Data/LT_all_site_level.csv", header=T)
d1$site_id <- d1$site
d1$River_type_fact <- as.factor(d1$River_type)
d1$Heavily_modified_fact <-as.factor(d1$Heavily_modified)
head(d1)
siteData <- unique(d1[,c("site_id","study_id","Country","season","TaxonomicRes", "Year_count", "Starting_year","Ending_year", "Sampling_years","River_type", "Heavily_modified", "Heavily_modified_code", "River_type_fact", "Heavily_modified_fact")])
response_gls_pivot <- merge(siteData,response_gls_pivot,by="site_id")
head(response_gls_pivot)
#write.csv(response_gls_pivot,"Outputs/GLS_trends_all.csv")

summaryData <- response_gls_pivot %>%
  group_by(River_type_fact, Heavily_modified_fact, Sampling_years, Year_count) %>%
  summarise(medTrends = median(flow),
            nuData = length(flow))

ggplot(summaryData)+
  geom_boxplot(aes(x=River_type_fact, y =medTrends),size=2)+
  theme_classic()

ggplot(summaryData)+
  geom_boxplot(aes(x=Heavily_modified_fact, y =medTrends),size=2)+
  theme_classic()

qplot(Year_count, medTrends, data=summaryData)

qplot(Sampling_years, medTrends, data=summaryData)

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
mycol <- t_col("pink", perc = 20, name = "lt.pink")

#####################################################
## effect of sampling years sampled on trend estimates
# tiff(filename = "Plots/Sampling_years_TrendEsts.tiff", width = 7, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mfrow=c(2,2),mar=c(4,4,0.4,0.4))
x1 <- response_gls_pivot$Sampling_years
y1 <- response_gls_pivot$flow
mod <- lm(y1 ~ x1)
summary(mod)
newx <- seq(min(x1), max(x1), length.out=54)
preds <- predict(mod, newdata = data.frame(x1=newx), interval = 'confidence')
# plot
plot(x=x1,y=y1,type="n",xlab="Years sampled", ylab="Taxon richness estimate")
polygon(x = c(1, -0.1, 2040, 2040), y = c(-100, 0, 0, -100), col ="grey80", border = NA)
box(lwd=2)
points(y1~x1)
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col =mycol, border = NA) #col = 'grey80'
# model
abline(mod,lty=1,lwd=2,col=2)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 2)
lines(newx, preds[ ,2], lty = 'dashed', col = 2)

x1 <- response_gls_pivot$Sampling_years
y1 <- response_gls_pivot$temp
mod <- lm(y1 ~ x1)
summary(mod)
newx <- seq(min(x1), max(x1), length.out=54)
preds <- predict(mod, newdata = data.frame(x1=newx), interval = 'confidence')
# plot
plot(x=x1,y=y1,type="n",xlab="Years sampled", ylab="Abundance estimate")
polygon(x = c(1, -0.1, 2040, 2040), y = c(-100, 0, 0, -100), col ="grey80", border = NA)
box(lwd=2)
points(y1~x1)
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col =mycol, border = NA) #col = 'grey80'
# model
abline(mod,lty=1,lwd=2,col=2)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 2)
lines(newx, preds[ ,2], lty = 'dashed', col = 2)

x1 <- response_gls_pivot$Sampling_years
y1 <- response_gls_pivot$PC_axis1
mod <- lm(y1 ~ x1)
summary(mod)
newx <- seq(min(x1), max(x1), length.out=54)
preds <- predict(mod, newdata = data.frame(x1=newx), interval = 'confidence')
# plot
plot(x=x1,y=y1,type="n",xlab="Years sampled", ylab="Functional richness estimate")
polygon(x = c(1, -0.1, 2040, 2040), y = c(-100, 0, 0, -100), col ="grey80", border = NA)
box(lwd=2)
points(y1~x1)
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col =mycol, border = NA) #col = 'grey80'
# model
abline(mod,lty=1,lwd=2,col=2)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 2)
lines(newx, preds[ ,2], lty = 'dashed', col = 2)

x1 <- response_gls_pivot$Sampling_years
y1 <- response_gls_pivot$PC_axis2
mod <- lm(y1 ~ x1)
summary(mod)
newx <- seq(min(x1), max(x1), length.out=54)
preds <- predict(mod, newdata = data.frame(x1=newx), interval = 'confidence')
# plot
plot(x=x1,y=y1,type="n",xlab="Years sampled", ylab="Functional redundancy estimate")
polygon(x = c(1, -0.1, 2040, 2040), y = c(-100, 0, 0, -100), col ="grey80", border = NA)
box(lwd=2)
points(y1~x1)
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col =mycol, border = NA) #col = 'grey80'
# model
abline(mod,lty=1,lwd=2,col=2)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 2)
lines(newx, preds[ ,2], lty = 'dashed', col = 2)

##
# dev.off()
##############################################################################

### meta-analysis ####
getwd()
library(rstan)
library(brms)
library(loo)

#### flow ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_noRandom_flow.rds")
loo_R2(fit)

#prob of trend
flow_prob <- getTrendProbability(fit)
flow_prob <- data.frame(Response="flow", flow_prob[,1:2])
flow_prob

#check model
plot(fit)
flow_loo <- loo(fit, cores = getOption("mc.cores", 1))
flow_loo
flow_parento <- as.list(pareto_k_table(flow_loo))
Count_flow <- rbind(flow_parento[[1]],flow_parento[[2]],flow_parento[[3]],flow_parento[[4]])
colnames(Count_flow) <- "flow"
pp_check(fit, ndraws = 100)

#pull out fixed effects
flow_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
flow_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
flow_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
flow_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
flow_fixed <- list(Response="flow", flow_fixed_995[,1:4], flow_fixed_975[,3:4],
                 flow_fixed_95[,3:4],flow_fixed_90[,3:4])
flow_fixed <-data.frame(lapply(flow_fixed, function(x) t(data.frame(x))))
flow_fixed 

#### temp ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_noRandom_temp.rds")
loo_R2(fit)

#prob of trend
temp_prob <- getTrendProbability(fit)
temp_prob <- data.frame(Response="temp", temp_prob[,1:2])
temp_prob

#check model
plot(fit)
temp_loo <- loo(fit, cores = getOption("mc.cores", 1))
temp_loo
temp_parento <- as.list(pareto_k_table(temp_loo))
Count_temp <- rbind(temp_parento[[1]],temp_parento[[2]],temp_parento[[3]],temp_parento[[4]])
colnames(Count_temp) <- "temp"
pp_check(fit, ndraws = 100)

#pull out fixed effects
temp_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
temp_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
temp_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
temp_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
temp_fixed <- list(Response="temp", temp_fixed_995[,1:4], temp_fixed_975[,3:4],
                   temp_fixed_95[,3:4],temp_fixed_90[,3:4])
temp_fixed <-data.frame(lapply(temp_fixed, function(x) t(data.frame(x))))
temp_fixed

#### PC_axis1 ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_noRandom_PC_axis1.rds")
loo_R2(fit)

#prob of trend
PC_axis1_prob <- getTrendProbability(fit)
PC_axis1_prob <- data.frame(Response="PC_axis1", PC_axis1_prob[,1:2])
PC_axis1_prob

#check model
plot(fit)
PC_axis1_loo <- loo(fit, cores = getOption("mc.cores", 1))
PC_axis1_loo
PC_axis1_parento <- as.list(pareto_k_table(PC_axis1_loo))
Count_PC_axis1 <- rbind(PC_axis1_parento[[1]],PC_axis1_parento[[2]],PC_axis1_parento[[3]],PC_axis1_parento[[4]])
colnames(Count_PC_axis1) <- "PC_axis1"
pp_check(fit, ndraws = 100)

#pull out fixed effects
PC_axis1_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
PC_axis1_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
PC_axis1_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
PC_axis1_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
PC_axis1_fixed <- list(Response="PC_axis1", PC_axis1_fixed_995[,1:4], PC_axis1_fixed_975[,3:4],
                   PC_axis1_fixed_95[,3:4],PC_axis1_fixed_90[,3:4])
PC_axis1_fixed <-data.frame(lapply(PC_axis1_fixed, function(x) t(data.frame(x))))
PC_axis1_fixed

#### PC_axis2 ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_noRandom_PC_axis2.rds")
loo_R2(fit)

#prob of trend
PC_axis2_prob <- getTrendProbability(fit)
PC_axis2_prob <- data.frame(Response="PC_axis2", PC_axis2_prob[,1:2])
PC_axis2_prob

#check model
plot(fit)
PC_axis2_loo <- loo(fit, cores = getOption("mc.cores", 1))
PC_axis2_loo
PC_axis2_parento <- as.list(pareto_k_table(PC_axis2_loo))
Count_PC_axis2 <- rbind(PC_axis2_parento[[1]],PC_axis2_parento[[2]],PC_axis2_parento[[3]],PC_axis2_parento[[4]])
colnames(Count_PC_axis2) <- "PC_axis2"
pp_check(fit, ndraws = 100)

#pull out fixed effects
PC_axis2_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
PC_axis2_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
PC_axis2_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
PC_axis2_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
PC_axis2_fixed <- list(Response="PC_axis2", PC_axis2_fixed_995[,1:4], PC_axis2_fixed_975[,3:4],
                   PC_axis2_fixed_95[,3:4],PC_axis2_fixed_90[,3:4])
PC_axis2_fixed <-data.frame(lapply(PC_axis2_fixed, function(x) t(data.frame(x))))
PC_axis2_fixed

#### assemble all model estimates from meta-analysis models #####

Yr_metaanaly_Ests <- rbind(flow_fixed, temp_fixed, PC_axis1_fixed, PC_axis2_fixed)
rownames(Yr_metaanaly_Ests) <- 1:4
write.csv(Yr_metaanaly_Ests, "Outputs/LT_Yr_metaanaly_weighted_noRandom_Ests_drivers.csv")

#### assemble all probabilities of increases/decreases from meta-analysis models #####

Yr_metaanaly_probs <- rbind(flow_prob, temp_prob, PC_axis1_prob, PC_axis2_prob)
write.csv(Yr_metaanaly_probs, "Outputs/LT_Yr_metaanaly_weighted_noRandom_probabilities_drivers.csv")

#### assemble model counts from Parento k diagnostic values from meta-analysis models #####

Yr_metaanaly_parento <- cbind(Count_flow, Count_temp, Count_PC_axis1, Count_PC_axis2)
rownames(Yr_metaanaly_parento) <- c("good[-Inf, 0.5]","ok[0.5, 0.7]","bad[0.7, 1]","verybad[1, Inf]")
write.csv(Yr_metaanaly_parento, "Outputs/LT_Yr_meta_parento_weighted_noRandom_ModelCounts_drivers.csv")

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
