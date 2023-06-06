# Import SHORT data
LT_taxa_SF <- read.csv("Data/LT_taxalist_2010-2020_short_final.csv", sep = ",", header = T, stringsAsFactors = F)
LT_CWM_SF <- read.csv("Data/LT_2010-2020_CWMs_short.csv", sep = ",", header = T, stringsAsFactors = F)
LT_Groups_SF <- read.csv("Data/LT_taxalist_2010-2020_pie_chart_SF_perc.csv", sep = ",", header = T, stringsAsFactors = F)

# Add required library
library(reshape2) #used to reshape the data

# Convert from SHORT to LONG form data
attach(LT_taxa_SF)
LT_taxa_LF <- melt(LT_taxa_SF, id = c("taxon_code", "taxon_name", "taxon_id", "ept", "epto", "coleoptera", "diptera", "insect", "mollusc", "group", "group2"))

attach(LT_CWM_SF)
LT_CWM_LF <- melt(LT_CWM_SF, id = c("sample_id"))

attach(LT_Groups_SF)
LT_Groups_LF <- melt(LT_Groups_SF, id = c("Phylum", "SubPhylum_Class", "SuperFam_Order"))

# Export LONG data
write.csv(LT_taxa_LF, "Data/LT_taxalist_2010-2020_long.csv")
write.csv(LT_CWM_LF, "Data/LT_2010-2020_CWMs_long.csv")
write.csv(LT_Groups_LF, "Data/LT_taxalist_2010-2020_pie_chart_LF_perc.csv")

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