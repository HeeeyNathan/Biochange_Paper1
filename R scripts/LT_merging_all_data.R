### merge all datasets
library(dplyr)

# read in data
## taxa diversity indices
taxa_div <- read.csv("Outputs/LT_siteYr_TaxaDiversity.csv", header = T)

## func diversity indices
func_div <- read.csv("Outputs/LT_siteYr_FuncDiversity.csv", header = T)
func_div <- func_div[, c(1:11)]

## func turnover
func_turn <- read.csv("Outputs/LT_siteYr_Functurnover.csv", header = T)
func_turn <- func_turn[order(func_turn$X),]
rownames(func_turn) <- 1:335
func_turn <- func_turn[, c(2, 18:20)]

## site factors
factors <- read.csv("Data/LT_siteYr_Factors.csv", header = T)

## environmental variables
env <- read.csv("Data/LT_siteYr_EnvVariables.csv", header =  TRUE)

# join taxonomic and functional datasets
LT_siteYr_all_variables <- dplyr::left_join(factors, taxa_div, by = "site_code")
LT_siteYr_all_variables <- dplyr::left_join(LT_siteYr_all_variables, func_div, by = "site_code")
LT_siteYr_all_variables <- dplyr::left_join(LT_siteYr_all_variables, func_turn, by = "site_code")
LT_siteYr_all_variables <- dplyr::left_join(LT_siteYr_all_variables, env, by = "site_code")

#### Save output cotaining no missing values
write.csv(LT_siteYr_all_variables, "Outputs/LT_siteYr_AllData.csv", row.names = FALSE)

# fill in missing year data
library(tidyr)

# generate data for year 2018 (all values missing)
data_2018 <- data.frame(site = rep(c(unique(LT_siteYr_all_variables$site_id)), each = 1), year = rep(2018, times = 41))

# create a new data frame with all possible combinations of Site and Year
full_data <- expand.grid(site = unique(LT_siteYr_all_variables$site_id), year = unique(LT_siteYr_all_variables$year))
full_data <- rbind(full_data, data_2018) # bind full and 2018 data
full_data <- arrange(full_data, year) # arrange data by year
full_data <- arrange(full_data, site) # arrange data by site
full_data$site_code <- paste(full_data$site, full_data$year, sep = "_") # create site_code column
full_data <- select(full_data, -site, -year) # remove redundant columns

# merge the full dataframe with the original data to fill in missing values
LT_siteYr_all_variables_wNAs <- merge(full_data, LT_siteYr_all_variables,  all.x = TRUE)

#### Save output cotaining missing values
write.csv(LT_siteYr_all_variables_wNAs, "Outputs/LT_siteYr_AllData_wNAs.csv", row.names = FALSE)

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