# script to combine site-level trends together in a single meta-analysis

### get response for this task ######

TaskID <- read.csv("Data/LT_ResponseTrends_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]

### get site-level values for this response ####

response_gls <- readRDS("outputs/glsTrends_site_level.rds")
response_gls <- subset(response_gls, Response == myResponse)
response_gls <- subset(response_gls, !is.na(estimate))

### site metadata ######

d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T) 
siteData <- unique(d1[,c("site_id", "country")])
response_gls <- merge(siteData,response_gls,by="site_id")

### run model ####

library(rstan)
library(brms)

### decide on priors ####

prior1 = c(set_prior("normal(0,3)", class = "Intercept"))

#examine response
hist(response_gls$estimate)
summary(response_gls$estimate)

#define weights
response_gls$w <- 1/response_gls$se
summary(response_gls$w)

# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
options(mc.cores = cpus_per_task)

# fit weighted model
fit1 <- brm(estimate|se(se) ~ 1,
            data = response_gls, iter=5000, init = 0,
            chains = 4, prior = prior1,
            control = list(adapt_delta = 0.90,
                           max_treedepth = 12))

summary(fit1)

plot(fit1)
resid <- residuals(fit1)
hist(resid, breaks = 30, main = "Histogram of Residuals")
qqnorm(resid)
qqline(resid, col = "red")
plot(fitted(fit1), resid, main = "Residuals vs. Fitted values")
abline(h = 0, lty = 2, col = "red")

### save output ####
saveRDS(fit1,file=paste0("Outputs/Metaanalysis_trends/metaanalysis_noRandom_",myResponse,".rds"))
#write.csv(fit1,file=paste0("Outputs/metaanalysis_noRandom_",myResponse,".csv"))

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

