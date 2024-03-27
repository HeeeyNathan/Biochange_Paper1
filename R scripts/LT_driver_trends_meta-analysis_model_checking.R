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

#### alkalinity ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_alkalinity.rds")
loo_R2(fit)

#prob of trend
alkalinity_prob <- getTrendProbability(fit)
alkalinity_prob <- data.frame(Response="alkalinity", alkalinity_prob[,1:2])
alkalinity_prob

#check model
plot(fit)
alkalinity_loo <- loo(fit, cores = getOption("mc.cores", 1))
alkalinity_loo
alkalinity_parento <- as.list(pareto_k_table(alkalinity_loo))
Count_alkalinity <- rbind(alkalinity_parento[[1]],alkalinity_parento[[2]],alkalinity_parento[[3]],alkalinity_parento[[4]])
colnames(Count_alkalinity) <- "alkalinity"
pp_check(fit, ndraws = 100)

#pull out fixed effects
alkalinity_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
alkalinity_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
alkalinity_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
alkalinity_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
alkalinity_fixed <- list(Response="alkalinity", alkalinity_fixed_995[,1:4], alkalinity_fixed_975[,3:4],
                    alkalinity_fixed_95[,3:4],alkalinity_fixed_90[,3:4])
alkalinity_fixed <-data.frame(lapply(alkalinity_fixed, function(x) t(data.frame(x))))
alkalinity_fixed

#### EC ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_EC.rds")
loo_R2(fit)

#prob of trend
EC_prob <- getTrendProbability(fit)
EC_prob <- data.frame(Response="EC", EC_prob[,1:2])
EC_prob

#check model
plot(fit)
EC_loo <- loo(fit, cores = getOption("mc.cores", 1))
EC_loo
EC_parento <- as.list(pareto_k_table(EC_loo))
Count_EC <- rbind(EC_parento[[1]],EC_parento[[2]],EC_parento[[3]],EC_parento[[4]])
colnames(Count_EC) <- "EC"
pp_check(fit, ndraws = 100)

#pull out fixed effects
EC_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
EC_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
EC_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
EC_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
EC_fixed <- list(Response="EC", EC_fixed_995[,1:4], EC_fixed_975[,3:4],
                    EC_fixed_95[,3:4],EC_fixed_90[,3:4])
EC_fixed <-data.frame(lapply(EC_fixed, function(x) t(data.frame(x))))
EC_fixed

#### NO3.N ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_NO3.N.rds")
loo_R2(fit)

#prob of trend
NO3.N_prob <- getTrendProbability(fit)
NO3.N_prob <- data.frame(Response="NO3.N", NO3.N_prob[,1:2])
NO3.N_prob

#check model
plot(fit)
NO3.N_loo <- loo(fit, cores = getOption("mc.cores", 1))
NO3.N_loo
NO3.N_parento <- as.list(pareto_k_table(NO3.N_loo))
Count_NO3.N <- rbind(NO3.N_parento[[1]],NO3.N_parento[[2]],NO3.N_parento[[3]],NO3.N_parento[[4]])
colnames(Count_NO3.N) <- "NO3.N"
pp_check(fit, ndraws = 100)

#pull out fixed effects
NO3.N_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
NO3.N_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
NO3.N_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
NO3.N_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
NO3.N_fixed <- list(Response="NO3.N", NO3.N_fixed_995[,1:4], NO3.N_fixed_975[,3:4],
                    NO3.N_fixed_95[,3:4],NO3.N_fixed_90[,3:4])
NO3.N_fixed <-data.frame(lapply(NO3.N_fixed, function(x) t(data.frame(x))))
NO3.N_fixed

#### NO2.N ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_NO2.N.rds")
loo_R2(fit)

#prob of trend
NO2.N_prob <- getTrendProbability(fit)
NO2.N_prob <- data.frame(Response="NO2.N", NO2.N_prob[,1:2])
NO2.N_prob

#check model
plot(fit)
NO2.N_loo <- loo(fit, cores = getOption("mc.cores", 1))
NO2.N_loo
NO2.N_parento <- as.list(pareto_k_table(NO2.N_loo))
Count_NO2.N <- rbind(NO2.N_parento[[1]],NO2.N_parento[[2]],NO2.N_parento[[3]],NO2.N_parento[[4]])
colnames(Count_NO2.N) <- "NO2.N"
pp_check(fit, ndraws = 100)

#pull out fixed effects
NO2.N_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
NO2.N_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
NO2.N_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
NO2.N_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
NO2.N_fixed <- list(Response="NO2.N", NO2.N_fixed_995[,1:4], NO2.N_fixed_975[,3:4],
                    NO2.N_fixed_95[,3:4],NO2.N_fixed_90[,3:4])
NO2.N_fixed <-data.frame(lapply(NO2.N_fixed, function(x) t(data.frame(x))))
NO2.N_fixed

#### mineral.N ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_mineral.N.rds")
loo_R2(fit)

#prob of trend
mineral.N_prob <- getTrendProbability(fit)
mineral.N_prob <- data.frame(Response="mineral.N", mineral.N_prob[,1:2])
mineral.N_prob

#check model
plot(fit)
mineral.N_loo <- loo(fit, cores = getOption("mc.cores", 1))
mineral.N_loo
mineral.N_parento <- as.list(pareto_k_table(mineral.N_loo))
Count_mineral.N <- rbind(mineral.N_parento[[1]],mineral.N_parento[[2]],mineral.N_parento[[3]],mineral.N_parento[[4]])
colnames(Count_mineral.N) <- "mineral.N"
pp_check(fit, ndraws = 100)

#pull out fixed effects
mineral.N_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
mineral.N_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
mineral.N_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
mineral.N_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
mineral.N_fixed <- list(Response="mineral.N", mineral.N_fixed_995[,1:4], mineral.N_fixed_975[,3:4],
                    mineral.N_fixed_95[,3:4],mineral.N_fixed_90[,3:4])
mineral.N_fixed <-data.frame(lapply(mineral.N_fixed, function(x) t(data.frame(x))))
mineral.N_fixed

#### Tot.N ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_Tot.N.rds")
loo_R2(fit)

#prob of trend
Tot.N_prob <- getTrendProbability(fit)
Tot.N_prob <- data.frame(Response="Tot.N", Tot.N_prob[,1:2])
Tot.N_prob

#check model
plot(fit)
Tot.N_loo <- loo(fit, cores = getOption("mc.cores", 1))
Tot.N_loo
Tot.N_parento <- as.list(pareto_k_table(Tot.N_loo))
Count_Tot.N <- rbind(Tot.N_parento[[1]],Tot.N_parento[[2]],Tot.N_parento[[3]],Tot.N_parento[[4]])
colnames(Count_Tot.N) <- "Tot.N"
pp_check(fit, ndraws = 100)

#pull out fixed effects
Tot.N_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
Tot.N_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
Tot.N_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
Tot.N_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
Tot.N_fixed <- list(Response="Tot.N", Tot.N_fixed_995[,1:4], Tot.N_fixed_975[,3:4],
                    Tot.N_fixed_95[,3:4],Tot.N_fixed_90[,3:4])
Tot.N_fixed <-data.frame(lapply(Tot.N_fixed, function(x) t(data.frame(x))))
Tot.N_fixed

#### PO4.P ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_PO4.P.rds")
loo_R2(fit)

#prob of trend
PO4.P_prob <- getTrendProbability(fit)
PO4.P_prob <- data.frame(Response="PO4.P", PO4.P_prob[,1:2])
PO4.P_prob

#check model
plot(fit)
PO4.P_loo <- loo(fit, cores = getOption("mc.cores", 1))
PO4.P_loo
PO4.P_parento <- as.list(pareto_k_table(PO4.P_loo))
Count_PO4.P <- rbind(PO4.P_parento[[1]],PO4.P_parento[[2]],PO4.P_parento[[3]],PO4.P_parento[[4]])
colnames(Count_PO4.P) <- "PO4.P"
pp_check(fit, ndraws = 100)

#pull out fixed effects
PO4.P_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
PO4.P_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
PO4.P_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
PO4.P_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
PO4.P_fixed <- list(Response="PO4.P", PO4.P_fixed_995[,1:4], PO4.P_fixed_975[,3:4],
                    PO4.P_fixed_95[,3:4],PO4.P_fixed_90[,3:4])
PO4.P_fixed <-data.frame(lapply(PO4.P_fixed, function(x) t(data.frame(x))))
PO4.P_fixed

#### Tot.P ####
fit <- readRDS("Outputs/Driver_metaanalysis_trends/metaanalysis_Tot.P.rds")
loo_R2(fit)

#prob of trend
Tot.P_prob <- getTrendProbability(fit)
Tot.P_prob <- data.frame(Response="Tot.P", Tot.P_prob[,1:2])
Tot.P_prob

#check model
plot(fit)
Tot.P_loo <- loo(fit, cores = getOption("mc.cores", 1))
Tot.P_loo
Tot.P_parento <- as.list(pareto_k_table(Tot.P_loo))
Count_Tot.P <- rbind(Tot.P_parento[[1]],Tot.P_parento[[2]],Tot.P_parento[[3]],Tot.P_parento[[4]])
colnames(Count_Tot.P) <- "Tot.P"
pp_check(fit, ndraws = 100)

#pull out fixed effects
Tot.P_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
Tot.P_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
Tot.P_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
Tot.P_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
Tot.P_fixed <- list(Response="Tot.P", Tot.P_fixed_995[,1:4], Tot.P_fixed_975[,3:4],
                    Tot.P_fixed_95[,3:4],Tot.P_fixed_90[,3:4])
Tot.P_fixed <-data.frame(lapply(Tot.P_fixed, function(x) t(data.frame(x))))
Tot.P_fixed

#### assemble all model estimates from meta-analysis models #####

Yr_metaanaly_Ests <- rbind(flow_fixed, temp_fixed, pH_fixed, o2_dis_fixed, NH4.N_fixed, PC_axis1_fixed,
                           alkalinity_fixed, EC_fixed, NO3.N_fixed, NO2.N_fixed, mineral.N_fixed, Tot.N_fixed, PO4.P_fixed, Tot.P_fixed)
rownames(Yr_metaanaly_Ests) <- 1:14
write.csv(Yr_metaanaly_Ests, "Outputs/LT_Yr_metaanaly_weighted_noRandom_Ests_drivers.csv")

#### assemble all probabilities of increases/decreases from meta-analysis models #####

Yr_metaanaly_probs <- rbind(flow_prob, temp_prob, pH_prob, o2_dis_prob, NH4.N_prob, PC_axis1_prob,
                            alkalinity_prob, EC_prob, NO3.N_prob, NO2.N_prob, mineral.N_prob, Tot.N_prob, PO4.P_prob, Tot.P_prob)
write.csv(Yr_metaanaly_probs, "Outputs/LT_Yr_metaanaly_weighted_noRandom_probabilities_drivers.csv")

#### assemble model counts from Parento k diagnostic values from meta-analysis models #####

Yr_metaanaly_parento <- cbind(Count_flow, Count_temp, Count_pH, Count_o2_dis, Count_NH4.N, Count_PC_axis1,
                              Count_alkalinity, Count_EC, Count_NO3.N, Count_NO2.N, Count_mineral.N, Count_Tot.N, Count_PO4.P, Count_Tot.P)
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
