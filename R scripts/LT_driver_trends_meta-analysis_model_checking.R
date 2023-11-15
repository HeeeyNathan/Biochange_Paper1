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
d1 <- read.csv("Data/LT_site_metadata.csv", header=T)
d1$site_id <- d1$site
d1$River_type_fact <- as.factor(d1$River_type)
d1$Heavily_modified_fact <-as.factor(d1$Heavily_modified)
head(d1)
siteData <- unique(d1[,c("site_id","study_id","Country","season","TaxonomicRes", "Year_count", "Starting_year","Ending_year", "Sampling_years","River_type", "Heavily_modified", "Heavily_modified_code", "River_type_fact", "Heavily_modified_fact")])
response_gls_pivot <- merge(siteData,response_gls_pivot,by="site_id")
head(response_gls_pivot)
# write.csv(response_gls_pivot,"Outputs/GLS_trends_all.csv")

### meta-analysis ####
getwd()
library(rstan)
library(brms)
library(loo)

#### flow ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_flow.rds")
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
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_temp.rds")
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

#### o2_dis ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_o2_dis.rds")
loo_R2(fit)

#prob of trend
o2_dis_prob <- getTrendProbability(fit)
o2_dis_prob <- data.frame(Response="o2_dis", o2_dis_prob[,1:2])
o2_dis_prob

#check model
plot(fit)
o2_dis_loo <- loo(fit, cores = getOption("mc.cores", 1))
o2_dis_loo
o2_dis_parento <- as.list(pareto_k_table(o2_dis_loo))
Count_o2_dis <- rbind(o2_dis_parento[[1]],o2_dis_parento[[2]],o2_dis_parento[[3]],o2_dis_parento[[4]])
colnames(Count_o2_dis) <- "o2_dis"
pp_check(fit, ndraws = 100)

#pull out fixed effects
o2_dis_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
o2_dis_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
o2_dis_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
o2_dis_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
o2_dis_fixed <- list(Response="o2_dis", o2_dis_fixed_995[,1:4], o2_dis_fixed_975[,3:4],
                   o2_dis_fixed_95[,3:4],o2_dis_fixed_90[,3:4])
o2_dis_fixed <-data.frame(lapply(o2_dis_fixed, function(x) t(data.frame(x))))
o2_dis_fixed

#### pH ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_pH.rds")
loo_R2(fit)

#prob of trend
pH_prob <- getTrendProbability(fit)
pH_prob <- data.frame(Response="pH", pH_prob[,1:2])
pH_prob

#check model
plot(fit)
pH_loo <- loo(fit, cores = getOption("mc.cores", 1))
pH_loo
pH_parento <- as.list(pareto_k_table(pH_loo))
Count_pH <- rbind(pH_parento[[1]],pH_parento[[2]],pH_parento[[3]],pH_parento[[4]])
colnames(Count_pH) <- "pH"
pp_check(fit, ndraws = 100)

#pull out fixed effects
pH_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
pH_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
pH_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
pH_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
pH_fixed <- list(Response="pH", pH_fixed_995[,1:4], pH_fixed_975[,3:4],
                   pH_fixed_95[,3:4],pH_fixed_90[,3:4])
pH_fixed <-data.frame(lapply(pH_fixed, function(x) t(data.frame(x))))
pH_fixed

#### NH4.N ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_NH4.N.rds")
loo_R2(fit)

#prob of trend
NH4.N_prob <- getTrendProbability(fit)
NH4.N_prob <- data.frame(Response="NH4.N", NH4.N_prob[,1:2])
NH4.N_prob

#check model
plot(fit)
NH4.N_loo <- loo(fit, cores = getOption("mc.cores", 1))
NH4.N_loo
NH4.N_parento <- as.list(pareto_k_table(NH4.N_loo))
Count_NH4.N <- rbind(NH4.N_parento[[1]],NH4.N_parento[[2]],NH4.N_parento[[3]],NH4.N_parento[[4]])
colnames(Count_NH4.N) <- "NH4.N"
pp_check(fit, ndraws = 100)

#pull out fixed effects
NH4.N_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
NH4.N_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
NH4.N_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
NH4.N_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
NH4.N_fixed <- list(Response="NH4.N", NH4.N_fixed_995[,1:4], NH4.N_fixed_975[,3:4],
                   NH4.N_fixed_95[,3:4],NH4.N_fixed_90[,3:4])
NH4.N_fixed <-data.frame(lapply(NH4.N_fixed, function(x) t(data.frame(x))))
NH4.N_fixed

#### PC_axis1 ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_PC_axis1.rds")
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

#### assemble all model estimates from meta-analysis models #####

Yr_metaanaly_Ests <- rbind(flow_fixed, temp_fixed, pH_fixed, o2_dis_fixed, NH4.N_fixed, PC_axis1_fixed)
rownames(Yr_metaanaly_Ests) <- 1:6
write.csv(Yr_metaanaly_Ests, "Outputs/LT_Yr_metaanaly_weighted_noRandom_Ests_drivers.csv")

#### assemble all probabilities of increases/decreases from meta-analysis models #####

Yr_metaanaly_probs <- rbind(flow_prob, temp_prob, pH_prob, o2_dis_prob, NH4.N_prob, PC_axis1_prob)
write.csv(Yr_metaanaly_probs, "Outputs/LT_Yr_metaanaly_weighted_noRandom_probabilities_drivers.csv")

#### assemble model counts from Parento k diagnostic values from meta-analysis models #####

Yr_metaanaly_parento <- cbind(Count_flow, Count_temp, Count_pH, Count_o2_dis, Count_NH4.N, Count_PC_axis1)
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
