#libraries
library(vegan)
library(codyn)
library(mobr)
library(reshape2)
library(dplyr)

###############################################################################################################

#load data in longform
all_lf <- read.csv("Data/LT_taxalist_2010-2020_long.csv", h = T, sep = ",", stringsAsFactors = FALSE, check.names=FALSE) ##this file is on OSF
head(all_lf)

#individual-based rarefication of species richness
all_lf$ro.ab <- round(all_lf$abundance, digits = 0) #rarefication only works on integers

TD <- NULL
for(i in unique(all_lf$site_id)){
  sub <- all_lf[all_lf$site_id == i, ]
  sub.m <- dcast(sub, site_code ~ taxon_name, sum, value.var = "abundance") 	# matrix form
  sub.ta <- subset(sub.m[,c(2:length(sub.m))])          # subset matrix to remove row names
  SppRich <- specnumber(sub.ta)  							          # taxonomic richness
  Simp <- diversity(sub.ta, index = "simpson") 					# Simpson's taxonomic diversity
  Shan <- diversity(sub.ta, index = "shannon")					# Shannon's taxonomic diversity
  EvenJ <- Shan/log(SppRich) 							    	        # Pielou's evenness (J)
  E10 <- Shan/SppRich 								  	              # Shannon's evenness (E10)
  Abund <- rowSums (sub.ta) 								            # Total abundance
  S_PIE <- calc_PIE(sub.ta, ENS = TRUE)						      # effective number of common species
  DATA1_Turnover <- turnover(sub, time.var = "year", species.var = "taxon_name", abundance.var = "abundance", replicate.var = NA, metric = "total")
  Turnover <- c("NA", DATA1_Turnover$total) 					  # Turnover per yr And first yr is "NA"
  DATA1_Turnover_app <- turnover(sub, time.var = "year", species.var = "taxon_name", abundance.var = "abundance", replicate.var = NA, metric = "appearance")
  Turnover_app <- c("NA", DATA1_Turnover_app$appearance) 					  # Turnover per yr And first yr is "NA"
  DATA1_Turnover_disapp <- turnover(sub, time.var = "year", species.var = "taxon_name", abundance.var = "abundance", replicate.var = NA, metric = "disappearance")
  Turnover_disapp <- c("NA", DATA1_Turnover_disapp$disappearance) 					  # Turnover per yr And first yr is "NA"
  sub.m_r <- dcast(sub, site_code ~ taxon_name, sum, value.var = "ro.ab")    # matrix form for rarefaction with rounded richness
  sub.ta_r <- subset(sub.m_r[,c(2:length(sub.m_r))])                  	# subset matrix to remove row names for rounded sppRich    	
  rare.sppRich <- if (min(rowSums(sub.ta_r)) > 10) {
    rarefy(sub.ta_r, sample = min(rowSums(sub.ta_r)))			# rarefy based on min abundance
  } else {rarefy(sub.ta_r, sample = 10)} 					# rarefy based on abund = 10 if min is less
  TD.i <- data.frame(sub.m$site_code, SppRich, Simp, Shan, EvenJ, E10, Abund, S_PIE, Turnover, Turnover_app, Turnover_disapp, rare.sppRich)
  TD <- rbind(TD, TD.i) ; rm(TD.i, sub.m, sub.ta, sub, SppRich, Simp, Shan, 
                             EvenJ, E10, Abund, S_PIE, DATA1_Turnover, Turnover, DATA1_Turnover_app, Turnover_app, DATA1_Turnover_disapp, Turnover_disapp, sub.m_r, sub.ta_r, rare.sppRich)
} ; rm(i)

colnames(TD)[1] <- "site_code"

#### create subset for EPT
head(all_lf)
sub_lf <- subset(all_lf, ept == "Yes")
head(sub_lf)

EPT_TD <- NULL
for(i in unique(sub_lf$site_id)){
  sub <- sub_lf[sub_lf$site_id == i, ]
  sub.m <- dcast(sub, site_code ~ taxon_name, sum, value.var = "abundance") # matrix form
  if(ncol(sub.m) > 2) {sub.ta <- subset(sub.m[,c(2:length(sub.m))])} else {sub.ta <- sub.m[, 2]} # subset matrix to remove row names
  (SppRich <- specnumber(sub.ta)) # taxonomic richness
  if(length(dim(sub.ta)) >= 2) {Abund <- rowSums (sub.ta)} else {Abund <- sub.ta} # Total abundance
  EPT_TD.i <- data.frame(sub.m$site_code, SppRich, Abund)
  EPT_TD <- rbind(EPT_TD, EPT_TD.i); rm(EPT_TD.i, sub.m, sub.ta, sub, SppRich, Abund)
} ; rm(i)

colnames(EPT_TD) <- c("site_code", "ept_spp_richness", "ept_abundance")

#### create subset for Crustacea
head(all_lf)
sub_lf <- subset(all_lf, crustacea == "Yes")
head(sub_lf)

CRU_TD <- NULL
for(i in unique(sub_lf$site_id)){
  sub <- sub_lf[sub_lf$site_id == i, ] 
  sub.m <- dcast(sub, site_code ~ taxon_name, sum, value.var = "abundance") # matrix form
  if(ncol(sub.m) > 2) {sub.ta <- subset(sub.m[,c(2:length(sub.m))])} else {sub.ta <- sub.m[, 2]} # subset matrix to remove row names
  (SppRich <- specnumber(sub.ta)) # taxonomic richness
  if(length(dim(sub.ta)) >= 2) {Abund <- rowSums (sub.ta)} else {Abund <- sub.ta} # Total abundance
  CRU_TD.i <- data.frame(sub.m$site_code, SppRich, Abund)
  CRU_TD <- rbind(CRU_TD, CRU_TD.i); rm(CRU_TD.i, sub.m, sub.ta, sub, SppRich, Abund)
} ; rm(i)

colnames(CRU_TD) <- c("site_code", "crustacea_spp_richness", "crustacea_abundance")

#### create subset for Diptera
head(all_lf)
sub_lf <- subset(all_lf, diptera == "Yes")
head(sub_lf)

DIP_TD <- NULL
for(i in unique(sub_lf$site_id)){
  sub <- sub_lf[sub_lf$site_id == i, ] 
  sub.m <- dcast(sub, site_code ~ taxon_name, sum, value.var = "abundance") # matrix form
  if(ncol(sub.m) > 2) {sub.ta <- subset(sub.m[,c(2:length(sub.m))])} else {sub.ta <- sub.m[, 2]} # subset matrix to remove row names
  (SppRich <- specnumber(sub.ta)) # taxonomic richness
  if(length(dim(sub.ta)) >= 2) {Abund <- rowSums (sub.ta)} else {Abund <- sub.ta} # Total abundance
  DIP_TD.i <- data.frame(sub.m$site_code, SppRich, Abund)
  DIP_TD <- rbind(DIP_TD, DIP_TD.i); rm(DIP_TD.i, sub.m, sub.ta, sub, SppRich, Abund)
} ; rm(i)

colnames(DIP_TD) <- c("site_code", "diptera_spp_richness", "diptera_abundance")

#### create subset for Molluscs
head(all_lf)
sub_lf <- subset(all_lf, mollusc == "Yes")
head(sub_lf)

MOL_TD <- NULL
for(i in unique(sub_lf$site_id)){
  sub <- sub_lf[sub_lf$site_id == i, ] 
  sub.m <- dcast(sub, site_code ~ taxon_name, sum, value.var = "abundance") # matrix form
  if(ncol(sub.m) > 2) {sub.ta <- subset(sub.m[,c(2:length(sub.m))])} else {sub.ta <- sub.m[, 2]} # subset matrix to remove row names
  (SppRich <- specnumber(sub.ta)) # taxonomic richness
  if(length(dim(sub.ta)) >= 2) {Abund <- rowSums (sub.ta)} else {Abund <- sub.ta} # Total abundance
  MOL_TD.i <- data.frame(sub.m$site_code, SppRich, Abund)
  MOL_TD <- rbind(MOL_TD, MOL_TD.i); rm(MOL_TD.i, sub.m, sub.ta, sub, SppRich, Abund)
} ; rm(i)

colnames(MOL_TD) <- c("site_code", "mollusc_spp_richness", "mollusc_abundance")

#### create subset for Insects minus EPT
# Note: I will remove EPT from insects to reduce trend overlap
head(all_lf)
sub_lf <- subset(all_lf, insect == "Yes")
head(sub_lf)

INS_TD <- NULL
for(i in unique(sub_lf$site_id)){
  sub <- sub_lf[sub_lf$site_id == i, ] 
  sub.m <- dcast(sub, site_code ~ taxon_name, sum, value.var = "abundance") # matrix form
  if(ncol(sub.m) > 2) {sub.ta <- subset(sub.m[,c(2:length(sub.m))])} else {sub.ta <- sub.m[, 2]} # subset matrix to remove row names
  (SppRich <- specnumber(sub.ta)) # taxonomic richness
  if(length(dim(sub.ta)) >= 2) {Abund <- rowSums (sub.ta)} else {Abund <- sub.ta} # Total abundance
  INS_TD.i <- data.frame(sub.m$site_code, SppRich, Abund)
  INS_TD <- rbind(INS_TD, INS_TD.i); rm(INS_TD.i, sub.m, sub.ta, sub, SppRich, Abund)
} ; rm(i)

colnames(INS_TD) <- c("site_code", "insect_spp_richness", "insect_abundance")

#### create subset for Annelids
head(all_lf)
sub_lf <- subset(all_lf, annelid == "Yes")
head(sub_lf)

ANN_TD <- NULL
for(i in unique(sub_lf$site_id)){
  sub <- sub_lf[sub_lf$site_id == i, ] 
  sub.m <- dcast(sub, site_code ~ taxon_name, sum, value.var = "abundance") # matrix form
  if(ncol(sub.m) > 2) {sub.ta <- subset(sub.m[,c(2:length(sub.m))])} else {sub.ta <- sub.m[, 2]} # subset matrix to remove row names
  (SppRich <- specnumber(sub.ta)) # taxonomic richness
  if(length(dim(sub.ta)) >= 2) {Abund <- rowSums (sub.ta)} else {Abund <- sub.ta} # Total abundance
  ANN_TD.i <- data.frame(sub.m$site_code, SppRich, Abund)
  ANN_TD <- rbind(ANN_TD, ANN_TD.i); rm(ANN_TD.i, sub.m, sub.ta, sub, SppRich, Abund)
} ; rm(i)

colnames(ANN_TD) <- c("site_code", "annelid_spp_richness", "annelid_abundance")

#### merge all taxonomic diversity datasets
TD_ALL <- dplyr::left_join(TD, EPT_TD, by = "site_code")
TD_ALL <- dplyr::left_join(TD_ALL, CRU_TD, by = "site_code")
TD_ALL <- dplyr::left_join(TD_ALL, DIP_TD, by = "site_code")
TD_ALL <- dplyr::left_join(TD_ALL, MOL_TD, by = "site_code")
TD_ALL <- dplyr::left_join(TD_ALL, INS_TD, by = "site_code")
TD_ALL <- dplyr::left_join(TD_ALL, ANN_TD, by = "site_code")

####################################################
write.csv(TD_ALL, "Outputs/LT_siteYr_TaxaDiversity.csv", row.names = FALSE)
####################################################

# ###  Calculate relative abundance of each species in each community
# # load data in short form
# all_sf <- read.csv("Data/LT_taxalist_2010-2020_short_final_transposed.csv", h = T, sep = ",", stringsAsFactors = FALSE, check.names=FALSE) ##this file is on OSF
# rownames(all_sf) <- all_sf$sample_id
# all_sf <- all_sf[, -1]
# all_sf <- as.matrix(all_sf)
# Relative_abund <- decostand(all_sf, method = 'total')
# mean_relative_abundance <- colMeans(Relative_abund)
# 
# write.csv(Relative_abund, "Outputs/Taxa_relative_abundance.csv")
# write.csv(mean_relative_abundance, "Outputs/Mean_taxa_relative_abundance.csv")

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
