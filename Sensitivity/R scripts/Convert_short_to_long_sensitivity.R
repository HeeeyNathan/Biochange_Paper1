# Import SHORT data
LT_taxa_SF <- read.csv("Sensitivity/Data/LT_taxalist_2010-2020_short_sensitivity.csv", sep = ",", header = T, stringsAsFactors = F)

# Add required library
library(reshape2) #used to reshape the data

# Convert from SHORT to LONG form data
attach(LT_taxa_SF)
LT_taxa_LF <- melt(LT_taxa_SF, id = c("taxon_code", "taxon_name", "taxon_id", "ept", "coleoptera", "diptera", "insect", "mollusc", "annelid", "crustacea", "group", "group2"))

# Export LONG data
write.csv(LT_taxa_LF, "Sensitivity/Data/LT_taxalist_2010-2020_long_sensitivity_updated.csv")

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