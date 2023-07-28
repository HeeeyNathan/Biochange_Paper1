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

response_gls <- readRDS("Sensitivity/Outputs/glsTrends_site_level_sensitivity.rds")

#pivot responses
response_gls_pivot <- response_gls %>%
  select(c(Response, estimate, site_id)) %>%
  pivot_wider(names_from = "Response",
              values_from = "estimate")

#get site metadata
d1 <- read.csv("Sensitivity/Data/LT_site_metadata_sensitivity.csv", header=T)
d1$site_id <- d1$site
d1$River_type_fact <- as.factor(d1$River_type)
d1$Heavily_modified_fact <-as.factor(d1$Heavily_modified)
head(d1)
siteData <- unique(d1[,c("site_id","study_id","Country","season","TaxonomicRes", "Year_count", "Starting_year","Ending_year", "Sampling_years","River_type", "Heavily_modified", "Heavily_modified_code", "River_type_fact", "Heavily_modified_fact")])
response_gls_pivot <- merge(siteData,response_gls_pivot,by="site_id")
head(response_gls_pivot)
#write.csv(response_gls_pivot,"Outputs/LT_site_trends.csv")

summaryData <- response_gls_pivot %>%
  group_by(River_type_fact, Heavily_modified_fact, Sampling_years, Year_count) %>%
  summarise(medTrends = median(spp_richness),
            nuData = length(spp_richness))

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
y1 <- response_gls_pivot$spp_richness
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
y1 <- response_gls_pivot$abundance
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
y1 <- response_gls_pivot$FRic
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
y1 <- response_gls_pivot$FRed
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

#### spp_richness ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_spp_richness.rds")
loo_R2(fit)

#prob of trend
sr_prob <- getTrendProbability(fit)
sr_prob <- data.frame(Response="spp_richness", sr_prob[,1:2])
sr_prob

#check model
plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo
sr_parento <- as.list(pareto_k_table(sr_loo))
Count_sr <- rbind(sr_parento[[1]],sr_parento[[2]],sr_parento[[3]],sr_parento[[4]])
colnames(Count_sr) <- "spp_richness"
pp_check(fit, ndraws = 100)

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- list(Response="spp_richness", sr_fixed_995[,1:4], sr_fixed_975[,3:4],
                 sr_fixed_95[,3:4],sr_fixed_90[,3:4])
sr_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))
sr_fixed 

#### spp_rich_rare ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_spp_rich_rare.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_spp_rich_rare.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_spp_rich_rare.rds")

#prob of trend
srr_prob <- getTrendProbability(fit)
srr_prob <- data.frame(Response="spp_rich_rare", srr_prob[,1:2])
srr_prob

#check model
plot(fit)
srr_loo <- loo(fit, cores = getOption("mc.cores", 1))
srr_loo
srr_parento <- as.list(pareto_k_table(srr_loo))
Count_srr <- rbind(srr_parento[[1]],srr_parento[[2]],srr_parento[[3]],srr_parento[[4]])
colnames(Count_srr) <- "spp_rich_rare"
pp_check(fit, ndraws = 100)

#pull out fixed effects
srr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
srr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
srr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
srr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
srr_fixed <- list(Response="spp_rich_rare", srr_fixed_995[,1:4], srr_fixed_975[,3:4],
                  srr_fixed_95[,3:4],srr_fixed_90[,3:4])
srr_fixed <-data.frame(lapply(srr_fixed, function(x) t(data.frame(x))))
srr_fixed

#### shannonsH ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_shannonsH.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_shannonsH.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_shannonsH.rds")

#prob of trend
shH_prob <- getTrendProbability(fit)
shH_prob <- data.frame(Response="shannonsH", shH_prob[,1:2])
shH_prob

#check model
plot(fit)
shH_loo <- loo(fit, cores = getOption("mc.cores", 1))
shH_loo
shH_parento <- as.list(pareto_k_table(shH_loo))
Count_shH <- rbind(shH_parento[[1]],shH_parento[[2]],shH_parento[[3]],shH_parento[[4]])
colnames(Count_shH) <- "ShannonsH"
pp_check(fit, ndraws = 100)

#pull out fixed effects
shH_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
shH_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
shH_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
shH_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
shH_fixed <- list(Response="shannonsH", shH_fixed_995[,1:4], shH_fixed_975[,3:4],
                  shH_fixed_95[,3:4],shH_fixed_90[,3:4])
shH_fixed <-data.frame(lapply(shH_fixed, function(x) t(data.frame(x))))
shH_fixed

#### E10 ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_E10.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_E10.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_E10.rds")

#prob of trend
e10_prob <- getTrendProbability(fit)
e10_prob <- data.frame(Response="E10", e10_prob[,1:2])
e10_prob

#check model
plot(fit)
e10_loo <- loo(fit, cores = getOption("mc.cores", 1))
e10_loo
e10_parento <- as.list(pareto_k_table(e10_loo))
Count_e10 <- rbind(e10_parento[[1]],e10_parento[[2]],e10_parento[[3]],e10_parento[[4]])
colnames(Count_e10) <- "E10"
pp_check(fit, ndraws = 100)

#pull out fixed effects
e10_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
e10_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
e10_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
e10_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
e10_fixed <- list(Response="E10", e10_fixed_995[,1:4], e10_fixed_975[,3:4],
                  e10_fixed_95[,3:4],e10_fixed_90[,3:4])
e10_fixed <-data.frame(lapply(e10_fixed, function(x) t(data.frame(x))))
e10_fixed

#### abundance ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_abundance.rds")

#prob of trend
ab_prob <- getTrendProbability(fit)
ab_prob <- data.frame(Response="abundance", ab_prob[,1:2])
ab_prob

#check model
plot(fit)
ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
ab_loo
ab_parento <- as.list(pareto_k_table(ab_loo))
Count_ab <- rbind(ab_parento[[1]],ab_parento[[2]],ab_parento[[3]],ab_parento[[4]])
colnames(Count_ab) <- "abundance"
pp_check(fit, ndraws = 100)

#pull out fixed effects
abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
abund_fixed <- list(Response="abundance", abund_fixed_995[,1:4], abund_fixed_975[,3:4],
                    abund_fixed_95[,3:4],abund_fixed_90[,3:4])
abund_fixed <-data.frame(lapply(abund_fixed, function(x) t(data.frame(x))))
abund_fixed

#### turnover ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_turnover.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_turnover.rds") 
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_turnover.rds")

#prob of trend
turn_prob <- getTrendProbability(fit)
turn_prob <- data.frame(Response="turnover", turn_prob[,1:2])
turn_prob

#check model
plot(fit)
turn_loo <- loo(fit, cores = getOption("mc.cores", 1))
turn_loo
turn_parento <- as.list(pareto_k_table(turn_loo))
Count_turn <- rbind(turn_parento[[1]],turn_parento[[2]],turn_parento[[3]],turn_parento[[4]])
colnames(Count_turn) <- "turnover"
pp_check(fit, ndraws = 100)

#pull out fixed effects
turn_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
turn_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
turn_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
turn_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
turn_fixed <- list(Response="turnover", turn_fixed_995[,1:4], turn_fixed_975[,3:4],
                   turn_fixed_95[,3:4],turn_fixed_90[,3:4])
turn_fixed <-data.frame(lapply(turn_fixed, function(x) t(data.frame(x))))
turn_fixed

#### F_to ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_F_turnover.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_F_turnover.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_F_turnover.rds")

#prob of trend
fto_prob <- getTrendProbability(fit)
fto_prob <- data.frame(Response="F_turnover", fto_prob[,1:2])
fto_prob

#check model
plot(fit)
fto_loo <- loo(fit, cores = getOption("mc.cores", 1))
fto_loo
fto_parento <- as.list(pareto_k_table(fto_loo))
Count_fto <- rbind(fto_parento[[1]],fto_parento[[2]],fto_parento[[3]],fto_parento[[4]])
colnames(Count_fto) <- "F_turnover"
pp_check(fit, ndraws = 100)

#pull out fixed effects
fto_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fto_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fto_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fto_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fto_fixed <- list(Response="F_turnover", fto_fixed_995[,1:4], fto_fixed_975[,3:4],
                  fto_fixed_95[,3:4],fto_fixed_90[,3:4])
fto_fixed <-data.frame(lapply(fto_fixed, function(x) t(data.frame(x))))
fto_fixed

#### FRic ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_FRic.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_FRic.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_FRic.rds")

#prob of trend
fric_prob <- getTrendProbability(fit)
fric_prob <- data.frame(Response="FRic", fric_prob[,1:2])
fric_prob

#check model
plot(fit)
fric_loo <- loo(fit, cores = getOption("mc.cores", 1))
fric_loo
fric_parento <- as.list(pareto_k_table(fric_loo))
Count_fric <- rbind(fric_parento[[1]],fric_parento[[2]],fric_parento[[3]],fric_parento[[4]])
colnames(Count_fric) <- "FRic"
pp_check(fit, ndraws = 100)

#pull out fixed effects
fric_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fric_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fric_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fric_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fric_fixed <- list(Response="FRic", fric_fixed_995[,1:4], fric_fixed_975[,3:4],
                   fric_fixed_95[,3:4],fric_fixed_90[,3:4])
fric_fixed <-data.frame(lapply(fric_fixed, function(x) t(data.frame(x))))
fric_fixed

#### FEve ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_FEve.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_FEve.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_FEve.rds")

#prob of trend
feve_prob <- getTrendProbability(fit)
feve_prob <- data.frame(Response="FEve", feve_prob[,1:2])
feve_prob

#check model
plot(fit)
feve_loo <- loo(fit, cores = getOption("mc.cores", 1))
feve_loo
feve_parento <- as.list(pareto_k_table(feve_loo))
Count_feve <- rbind(feve_parento[[1]],feve_parento[[2]],feve_parento[[3]],feve_parento[[4]])
colnames(Count_feve) <- "FEve"
pp_check(fit, ndraws = 100)

#pull out fixed effects
feve_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
feve_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
feve_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
feve_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
feve_fixed <- list(Response="FEve", feve_fixed_995[,1:4], feve_fixed_975[,3:4],
                   feve_fixed_95[,3:4],feve_fixed_90[,3:4])
feve_fixed <-data.frame(lapply(feve_fixed, function(x) t(data.frame(x))))
feve_fixed

#### FDis ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_FDis.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_FDis.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_FDis.rds")

#prob of trend
fdis_prob <- getTrendProbability(fit)
fdis_prob <- data.frame(Response="FDis", fdis_prob[,1:2])
fdis_prob

#check model
plot(fit)
fdis_loo <- loo(fit, cores = getOption("mc.cores", 1))
fdis_loo
fdis_parento <- as.list(pareto_k_table(fdis_loo))
Count_fdis <- rbind(fdis_parento[[1]],fdis_parento[[2]],fdis_parento[[3]],fdis_parento[[4]])
colnames(Count_fdis) <- "FDis"
pp_check(fit, ndraws = 100)

#pull out fixed effects
fdis_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fdis_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fdis_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fdis_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fdis_fixed <- list(Response="FDis", fdis_fixed_995[,1:4], fdis_fixed_975[,3:4],
                   fdis_fixed_95[,3:4],fdis_fixed_90[,3:4])
fdis_fixed <-data.frame(lapply(fdis_fixed, function(x) t(data.frame(x))))
fdis_fixed

#### FRed ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_FRed.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_FRed.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_FRed.rds")

#prob of trend
FRed_prob <- getTrendProbability(fit)
FRed_prob <- data.frame(Response="FRed", FRed_prob[,1:2])
FRed_prob

#check model
plot(fit)
FRed_loo <- loo(fit, cores = getOption("mc.cores", 1))
FRed_loo
FRed_parento <- as.list(pareto_k_table(FRed_loo))
Count_FRed <- rbind(FRed_parento[[1]],FRed_parento[[2]],FRed_parento[[3]],FRed_parento[[4]])
colnames(Count_FRed) <- "FRed"
pp_check(fit, ndraws = 100)

#pull out fixed effects
FRed_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
FRed_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
FRed_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
FRed_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
FRed_fixed <- list(Response="FRed", FRed_fixed_995[,1:4], FRed_fixed_975[,3:4],
                   FRed_fixed_95[,3:4],FRed_fixed_90[,3:4])
FRed_fixed <-data.frame(lapply(FRed_fixed, function(x) t(data.frame(x))))
FRed_fixed

#### ept_spp_richness ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_ept_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_spp_richness.rds")
loo_R2(fit)

#prob of trend
ept_sr_prob <- getTrendProbability(fit)
ept_sr_prob <- data.frame(Response="ept_richness", ept_sr_prob[,1:2])
ept_sr_prob

#check model
plot(fit)
ept_sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
ept_sr_loo
ept_sr_parento <- as.list(pareto_k_table(ept_sr_loo))
Count_ept_sr <- rbind(ept_sr_parento[[1]],ept_sr_parento[[2]],ept_sr_parento[[3]],ept_sr_parento[[4]])
colnames(Count_ept_sr) <- "ept_richness"
pp_check(fit, ndraws = 100)

#pull out fixed effects
ept_sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
ept_sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
ept_sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
ept_sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
ept_sr_fixed <- list(Response="ept_richness", ept_sr_fixed_995[,1:4], ept_sr_fixed_975[,3:4],
                 ept_sr_fixed_95[,3:4],ept_sr_fixed_90[,3:4])
ept_sr_fixed <-data.frame(lapply(ept_sr_fixed, function(x) t(data.frame(x))))
ept_sr_fixed 

#### ept_abundance ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_ept_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_abundance.rds")

#prob of trend
ept_ab_prob <- getTrendProbability(fit)
ept_ab_prob <- data.frame(Response="ept_abundance", ept_ab_prob[,1:2])
ept_ab_prob

#check model
plot(fit)
ept_ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
ept_ab_loo
ept_ab_parento <- as.list(pareto_k_table(ept_ab_loo))
Count_ept_ab <- rbind(ept_ab_parento[[1]],ept_ab_parento[[2]],ept_ab_parento[[3]],ept_ab_parento[[4]])
colnames(Count_ept_ab) <- "ept_abundance"
pp_check(fit, ndraws = 100)

#pull out fixed effects
ept_abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
ept_abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
ept_abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
ept_abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
ept_abund_fixed <- list(Response="ept_abundance", ept_abund_fixed_995[,1:4], ept_abund_fixed_975[,3:4],
                    ept_abund_fixed_95[,3:4],ept_abund_fixed_90[,3:4])
ept_abund_fixed <-data.frame(lapply(ept_abund_fixed, function(x) t(data.frame(x))))
ept_abund_fixed

#### diptera_spp_richness ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_diptera_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_spp_richness.rds")
loo_R2(fit)

#prob of trend
diptera_sr_prob <- getTrendProbability(fit)
diptera_sr_prob <- data.frame(Response="diptera_richness", diptera_sr_prob[,1:2])
diptera_sr_prob

#check model
plot(fit)
diptera_sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
diptera_sr_loo
diptera_sr_parento <- as.list(pareto_k_table(diptera_sr_loo))
Count_diptera_sr <- rbind(diptera_sr_parento[[1]],diptera_sr_parento[[2]],diptera_sr_parento[[3]],diptera_sr_parento[[4]])
colnames(Count_diptera_sr) <- "diptera_richness"
pp_check(fit, ndraws = 100)

#pull out fixed effects
diptera_sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
diptera_sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
diptera_sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
diptera_sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
diptera_sr_fixed <- list(Response="diptera_richness", diptera_sr_fixed_995[,1:4], diptera_sr_fixed_975[,3:4],
                     diptera_sr_fixed_95[,3:4],diptera_sr_fixed_90[,3:4])
diptera_sr_fixed <-data.frame(lapply(diptera_sr_fixed, function(x) t(data.frame(x))))
diptera_sr_fixed 

#### diptera_abundance ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_diptera_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_abundance.rds")

#prob of trend
diptera_ab_prob <- getTrendProbability(fit)
diptera_ab_prob <- data.frame(Response="diptera_abundance", diptera_ab_prob[,1:2])
diptera_ab_prob

#check model
plot(fit)
diptera_ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
diptera_ab_loo
diptera_ab_parento <- as.list(pareto_k_table(diptera_ab_loo))
Count_diptera_ab <- rbind(diptera_ab_parento[[1]],diptera_ab_parento[[2]],diptera_ab_parento[[3]],diptera_ab_parento[[4]])
colnames(Count_diptera_ab) <- "diptera_abundance"
pp_check(fit, ndraws = 100)

#pull out fixed effects
diptera_abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
diptera_abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
diptera_abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
diptera_abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
diptera_abund_fixed <- list(Response="diptera_abundance", diptera_abund_fixed_995[,1:4], diptera_abund_fixed_975[,3:4],
                        diptera_abund_fixed_95[,3:4],diptera_abund_fixed_90[,3:4])
diptera_abund_fixed <-data.frame(lapply(diptera_abund_fixed, function(x) t(data.frame(x))))
diptera_abund_fixed

#### insect_richness ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_insect_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_richness.rds")

#prob of trend
insect_sr_prob <- getTrendProbability(fit)
insect_sr_prob <- data.frame(Response="insect_richness", insect_sr_prob[,1:2])
insect_sr_prob

#check model
plot(fit)
insect_sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
insect_sr_loo
insect_sr_parento <- as.list(pareto_k_table(insect_sr_loo))
Count_insect_sr <- rbind(insect_sr_parento[[1]],insect_sr_parento[[2]],insect_sr_parento[[3]],insect_sr_parento[[4]])
colnames(Count_insect_sr) <- "insect_richness"
pp_check(fit, ndraws = 100)

#pull out fixed effects
insect_sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
insect_sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insect_sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insect_sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insect_sr_fixed <- list(Response="insect_richness", insect_sr_fixed_995[,1:4], insect_sr_fixed_975[,3:4],
                           insect_sr_fixed_95[,3:4],insect_sr_fixed_90[,3:4])
insect_sr_fixed <-data.frame(lapply(insect_sr_fixed, function(x) t(data.frame(x))))
insect_sr_fixed

#### insect_abundance ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_insect_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_abundance.rds")

#prob of trend
insect_ab_prob <- getTrendProbability(fit)
insect_ab_prob <- data.frame(Response="insect_abundance", insect_ab_prob[,1:2])
insect_ab_prob

#check model
plot(fit)
insect_ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
insect_ab_loo
insect_ab_parento <- as.list(pareto_k_table(insect_ab_loo))
Count_insect_ab <- rbind(insect_ab_parento[[1]],insect_ab_parento[[2]],insect_ab_parento[[3]],insect_ab_parento[[4]])
colnames(Count_insect_ab) <- "insect_abundance"
pp_check(fit, ndraws = 100)

#pull out fixed effects
insect_abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
insect_abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insect_abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insect_abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insect_abund_fixed <- list(Response="insect_abundance", insect_abund_fixed_995[,1:4], insect_abund_fixed_975[,3:4],
                            insect_abund_fixed_95[,3:4],insect_abund_fixed_90[,3:4])
insect_abund_fixed <-data.frame(lapply(insect_abund_fixed, function(x) t(data.frame(x))))
insect_abund_fixed

#### mollusc_richness ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_mollusc_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_richness.rds")

#prob of trend
mollusc_sr_prob <- getTrendProbability(fit)
mollusc_sr_prob <- data.frame(Response="mollusc_richness", mollusc_sr_prob[,1:2])
mollusc_sr_prob

#check model
plot(fit)
mollusc_sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
mollusc_sr_loo
mollusc_sr_parento <- as.list(pareto_k_table(mollusc_sr_loo))
Count_mollusc_sr <- rbind(mollusc_sr_parento[[1]],mollusc_sr_parento[[2]],mollusc_sr_parento[[3]],mollusc_sr_parento[[4]])
colnames(Count_mollusc_sr) <- "mollusc_richness"
pp_check(fit, ndraws = 100)

#pull out fixed effects
mollusc_sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
mollusc_sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
mollusc_sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
mollusc_sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
mollusc_sr_fixed <- list(Response="mollusc_richness", mollusc_sr_fixed_995[,1:4], mollusc_sr_fixed_975[,3:4],
                        mollusc_sr_fixed_95[,3:4],mollusc_sr_fixed_90[,3:4])
mollusc_sr_fixed <-data.frame(lapply(mollusc_sr_fixed, function(x) t(data.frame(x))))
mollusc_sr_fixed

#### mollusc_abundance ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_mollusc_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_abundance.rds")

#prob of trend
mollusc_ab_prob <- getTrendProbability(fit)
mollusc_ab_prob <- data.frame(Response="mollusc_abundance", mollusc_ab_prob[,1:2])
mollusc_ab_prob

#check model
plot(fit)
mollusc_ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
mollusc_ab_loo
mollusc_ab_parento <- as.list(pareto_k_table(mollusc_ab_loo))
Count_mollusc_ab <- rbind(mollusc_ab_parento[[1]],mollusc_ab_parento[[2]],mollusc_ab_parento[[3]],mollusc_ab_parento[[4]])
colnames(Count_mollusc_ab) <- "mollusc_abundance"
pp_check(fit, ndraws = 100)

#pull out fixed effects
mollusc_abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
mollusc_abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
mollusc_abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
mollusc_abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
mollusc_abund_fixed <- list(Response="mollusc_abundance", mollusc_abund_fixed_995[,1:4], mollusc_abund_fixed_975[,3:4],
                           mollusc_abund_fixed_95[,3:4],mollusc_abund_fixed_90[,3:4])
mollusc_abund_fixed <-data.frame(lapply(mollusc_abund_fixed, function(x) t(data.frame(x))))
mollusc_abund_fixed

#### annelid_richness ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_annelid_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_richness.rds")

#prob of trend
annelid_sr_prob <- getTrendProbability(fit)
annelid_sr_prob <- data.frame(Response="annelid_richness", annelid_sr_prob[,1:2])
annelid_sr_prob

#check model
plot(fit)
annelid_sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
annelid_sr_loo
annelid_sr_parento <- as.list(pareto_k_table(annelid_sr_loo))
Count_annelid_sr <- rbind(annelid_sr_parento[[1]],annelid_sr_parento[[2]],annelid_sr_parento[[3]],annelid_sr_parento[[4]])
colnames(Count_annelid_sr) <- "annelid_richness"
pp_check(fit, ndraws = 100)

#pull out fixed effects
annelid_sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
annelid_sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
annelid_sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
annelid_sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
annelid_sr_fixed <- list(Response="annelid_richness", annelid_sr_fixed_995[,1:4], annelid_sr_fixed_975[,3:4],
                         annelid_sr_fixed_95[,3:4],annelid_sr_fixed_90[,3:4])
annelid_sr_fixed <-data.frame(lapply(annelid_sr_fixed, function(x) t(data.frame(x))))
annelid_sr_fixed

#### annelid_abundance ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_annelid_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_abundance.rds")

#prob of trend
annelid_ab_prob <- getTrendProbability(fit)
annelid_ab_prob <- data.frame(Response="annelid_abundance", annelid_ab_prob[,1:2])
annelid_ab_prob

#check model
plot(fit)
annelid_ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
annelid_ab_loo
annelid_ab_parento <- as.list(pareto_k_table(annelid_ab_loo))
Count_annelid_ab <- rbind(annelid_ab_parento[[1]],annelid_ab_parento[[2]],annelid_ab_parento[[3]],annelid_ab_parento[[4]])
colnames(Count_annelid_ab) <- "annelid_abundance"
pp_check(fit, ndraws = 100)

#pull out fixed effects
annelid_abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
annelid_abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
annelid_abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
annelid_abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
annelid_abund_fixed <- list(Response="annelid_abundance", annelid_abund_fixed_995[,1:4], annelid_abund_fixed_975[,3:4],
                            annelid_abund_fixed_95[,3:4],annelid_abund_fixed_90[,3:4])
annelid_abund_fixed <-data.frame(lapply(annelid_abund_fixed, function(x) t(data.frame(x))))
annelid_abund_fixed

#### crustacea_richness ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_crustacea_spp_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_richness.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_richness.rds")

#prob of trend
crustacea_sr_prob <- getTrendProbability(fit)
crustacea_sr_prob <- data.frame(Response="crustacea_richness", crustacea_sr_prob[,1:2])
crustacea_sr_prob

#check model
plot(fit)
crustacea_sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
crustacea_sr_loo
crustacea_sr_parento <- as.list(pareto_k_table(crustacea_sr_loo))
Count_crustacea_sr <- rbind(crustacea_sr_parento[[1]],crustacea_sr_parento[[2]],crustacea_sr_parento[[3]],crustacea_sr_parento[[4]])
colnames(Count_crustacea_sr) <- "crustacea_richness"
pp_check(fit, ndraws = 100)

#pull out fixed effects
crustacea_sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
crustacea_sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
crustacea_sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
crustacea_sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
crustacea_sr_fixed <- list(Response="crustacea_richness", crustacea_sr_fixed_995[,1:4], crustacea_sr_fixed_975[,3:4],
                         crustacea_sr_fixed_95[,3:4],crustacea_sr_fixed_90[,3:4])
crustacea_sr_fixed <-data.frame(lapply(crustacea_sr_fixed, function(x) t(data.frame(x))))
crustacea_sr_fixed

#### crustacea_abundance ####
fit <- readRDS("Sensitivity/Outputs/Metaanalysis_trends/metaanalysis_noRandom_crustacea_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Unweighted_wRandom/metaanalysis_unweighted_abundance.rds")
#fit <- readRDS("Outputs/Metaanalysis_trends/Weighted_wRandom/metaanalysis_abundance.rds")

#prob of trend
crustacea_ab_prob <- getTrendProbability(fit)
crustacea_ab_prob <- data.frame(Response="crustacea_abundance", crustacea_ab_prob[,1:2])
crustacea_ab_prob

#check model
plot(fit)
crustacea_ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
crustacea_ab_loo
crustacea_ab_parento <- as.list(pareto_k_table(crustacea_ab_loo))
Count_crustacea_ab <- rbind(crustacea_ab_parento[[1]],crustacea_ab_parento[[2]],crustacea_ab_parento[[3]],crustacea_ab_parento[[4]])
colnames(Count_crustacea_ab) <- "crustacea_abundance"
pp_check(fit, ndraws = 100)

#pull out fixed effects
crustacea_abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
crustacea_abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
crustacea_abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
crustacea_abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
crustacea_abund_fixed <- list(Response="crustacea_abundance", crustacea_abund_fixed_995[,1:4], crustacea_abund_fixed_975[,3:4],
                            crustacea_abund_fixed_95[,3:4],crustacea_abund_fixed_90[,3:4])
crustacea_abund_fixed <-data.frame(lapply(crustacea_abund_fixed, function(x) t(data.frame(x))))
crustacea_abund_fixed

#### assemble all model estimates from meta-analysis models #####

Yr_metaanaly_Ests <- rbind(sr_fixed, srr_fixed, shH_fixed, e10_fixed, abund_fixed, turn_fixed, 
                           fto_fixed, fric_fixed, feve_fixed, fdis_fixed, FRed_fixed,
                           ept_sr_fixed, ept_abund_fixed, 
                           diptera_sr_fixed, diptera_abund_fixed, 
                           insect_sr_fixed, insect_abund_fixed, 
                           mollusc_sr_fixed, mollusc_abund_fixed,
                           annelid_sr_fixed, annelid_abund_fixed,
                           crustacea_sr_fixed, crustacea_abund_fixed)
rownames(Yr_metaanaly_Ests) <- 1:23
write.csv(Yr_metaanaly_Ests, "Sensitivity/Outputs/LT_Yr_metaanaly_Ests_sensitivity.csv")

#### assemble all probabilities of increases/decreases from meta-analysis models #####

Yr_metaanaly_probs <- rbind(sr_prob, srr_prob, shH_prob, e10_prob, ab_prob, turn_prob, 
                            fto_prob, fric_prob, feve_prob, fdis_prob, FRed_prob,
                            ept_sr_prob, ept_ab_prob, 
                            diptera_sr_prob, diptera_ab_prob,
                            insect_sr_prob, insect_ab_prob,
                            mollusc_sr_prob, mollusc_ab_prob,
                            annelid_sr_prob, annelid_ab_prob,
                            crustacea_sr_prob, crustacea_ab_prob)
write.csv(Yr_metaanaly_probs, "Sensitivity/Outputs/LT_Yr_metaanaly_probabilities_sensitivity.csv")

#### assemble model counts from Parento k diagnostic values from meta-analysis models #####

Yr_metaanaly_parento <- cbind(Count_sr, Count_srr, Count_shH, Count_e10, Count_ab, Count_turn, Count_fto,
                              Count_fric, Count_feve, Count_fdis, Count_FRed,
                              Count_ept_sr, Count_ept_ab,
                              Count_diptera_sr, Count_diptera_ab,
                              Count_insect_sr, Count_insect_ab,
                              Count_mollusc_sr, Count_mollusc_ab,
                              Count_annelid_sr, Count_annelid_ab,
                              Count_crustacea_sr, Count_crustacea_ab)
rownames(Yr_metaanaly_parento) <- c("good[-Inf, 0.5]","ok[0.5, 0.7]","bad[0.7, 1]","verybad[1, Inf]")
write.csv(Yr_metaanaly_parento, "Sensitivity/Outputs/LT_Yr_meta_parento_ModelCounts_sensitivity.csv")

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
