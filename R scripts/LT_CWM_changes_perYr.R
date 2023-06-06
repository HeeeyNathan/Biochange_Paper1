# library
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(cowplot)
library(tidyverse)
library(ggalluvial)

# load data
CWM <- read.csv("Data/LT_2010-2020_CWMs_short.csv")

#text to column function to split the site_code column into 2
CWM <- CWM %>%
  separate(site_code, into = c("site_id", "year"), sep = "_")

# calculate means for all traits 
CWM_Yr_means <- CWM %>% 
  summarise(across(2:64, mean), .by = year) %>%
  arrange(year)

CWM_Yr_means$year <- factor(CWM_Yr_means$year)

# Remove the year 2019 from all non-within-site analyses
CWM_Yr_means <- CWM_Yr_means[CWM_Yr_means$year != "2011",] # removes the intercept line from each response variable
CWM_Yr_means <- CWM_Yr_means[CWM_Yr_means$year != "2019",] # removes the intercept line from each response variable

# create trait groups
tachaqua <- select(CWM_Yr_means, c("year", contains("tachaqua")))
tachdisp <- select(CWM_Yr_means, c("year", contains("tachdisp")))
tachfeed <- select(CWM_Yr_means, c("year", contains("tachfeed")))
tachfood <- select(CWM_Yr_means, c("year", contains("tachfood")))
tachlifedur <- select(CWM_Yr_means, c("year", contains("tachlifedur")))
tachloco <- select(CWM_Yr_means, c("year", contains("tachloco")))
tachsize <- select(CWM_Yr_means, c("year", contains("tachsize")))
tachrepcyc <- select(CWM_Yr_means, c("year", contains("tachrepcyc")))
tachrepro <- select(CWM_Yr_means, c("year", contains("tachrepro")))
tachresist <- select(CWM_Yr_means, c("year", contains("tachresist")))
tachrespir <- select(CWM_Yr_means, c("year", contains("tachrespir")))

# covert groups to long format
tachaqua_L <- gather(tachaqua, trait, meanCWM, 2:5, factor_key=F) %>% 
  arrange(year)
tachdisp_L <- gather(tachdisp, trait, meanCWM, 2:5, factor_key=F) %>% 
  arrange(year)
tachfeed_L <- gather(tachfeed, trait, meanCWM, 2:9, factor_key=F) %>% 
  arrange(year)
tachfood_L <- gather(tachfood, trait, meanCWM, 2:10, factor_key=F) %>% 
  arrange(year)
tachlifedur_L <- gather(tachlifedur, trait, meanCWM, 2:3, factor_key=F) %>% 
  arrange(year)
tachloco_L <- gather(tachloco, trait, meanCWM, 2:9, factor_key=F) %>% 
  arrange(year)
tachsize_L <- gather(tachsize, trait, meanCWM, 2:8, factor_key=F) %>% 
  arrange(year)
tachrepcyc_L <- gather(tachrepcyc, trait, meanCWM, 2:4, factor_key=F) %>% 
  arrange(year)
tachrepro_L <- gather(tachrepro, trait, meanCWM, 2:9, factor_key=F) %>% 
  arrange(year)
tachresist_L <- gather(tachresist, trait, meanCWM, 2:6, factor_key=F) %>% 
  arrange(year)
tachrespir_L <- gather(tachrespir, trait, meanCWM, 2:6, factor_key=F) %>% 
  arrange(year)

# Stacked + percent
# tachaqua
nice_cols <- brewer.pal(4, "Pastel1")
p1 <- ggplot(tachaqua_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c("Adult", 'Egg', 'Larva', 'Nymph'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "") +
  ggtitle("Aquatic stages (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p1


# tachdisp
nice_cols <- brewer.pal(4, "Pastel1")
p2 <- ggplot(tachdisp_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c("Aerial active", 'Aerial passive', 'Aquatic active', 'Aquatic passive'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Dispersal", face = "bold") +
  ggtitle("Dispersal (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p2

# tachfeed
nice_cols <- brewer.pal(8, "Pastel1")
p3 <- ggplot(tachfeed_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c("Absorbers", 'Deposit feeders', 'Filter-feeders', 'Parasites', 'Piercers', 'Predators', 'Scrapers', 'Shredders'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Feeding habits", face = "bold") +
  ggtitle("Feeding habits (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p3

# tachfood
nice_cols <- brewer.pal(9, "Pastel1")
p4 <- ggplot(tachfood_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c("Dead animals (> 1mm)", 'Detritus (< 1mm)', 'Dead plants (>= 1mm)', 'Living macrophytes', 'Living microphytes', 'Living macroinvertebrates', 'Living microinvertebrates', 'Microorganisms', 'Vertebrates'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Food", face = "bold") +
  ggtitle("Food (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p4

# tachlifedur
nice_cols <- c("#FBB4AE", "#B3CDE3")
p5 <- ggplot(tachlifedur_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c('Up to one year', 'Longer than one year'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Lifecycle duration", face = "bold") +
  ggtitle("Lifecycle duration (Tachet et al. 2010)") +
  ylab("Percentage (%)") +
  xlab("Year")
p5

# tachloco
nice_cols <- brewer.pal(8, "Pastel1")
p6 <- ggplot(tachloco_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c('Burrower (epibenthic)', 'Crawler', 'Flier', 'Interstitial (endobenthic)', 'Permanently attached', 'Surface swimmer', 'Temporarily attached', 'Full water swimmer'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Locomotion &\nsubstrate relation", face = "bold") +
  ggtitle("Locomotion & substrate relation (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p6

# tachsize
nice_cols <- brewer.pal(7, "Pastel1")
p7 <- ggplot(tachsize_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c('> 0.25 - 0.5 cm', '> 0.5 - 1 cm', '> 1 - 2 cm', '> 2 - 4 cm', '> 4 - 8 cm', '> 8 cm', '<= 0.25 cm'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Maximum potential body\nsize", face = "bold") +
  ggtitle("Maximum potential body size (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p7

# tachrepcyc
nice_cols <- brewer.pal(3, "Pastel1")
p8 <- ggplot(tachrepcyc_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c('Multivoltine (> 1)', 'Semivoltine (< 1)', 'Univoltine (1)'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Potential no. of\nreproduction cycles\nper year", face = "bold") +
  ggtitle("Potential no. of reproduction cycles per year (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p8

# tachrepro
nice_cols <- brewer.pal(8, "Pastel1")
p9 <- ggplot(tachrepro_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c('Asexual reproduction', 'Isolated eggs, cemented', 'Clutches, cemented or fixed', 'Isolated eggs, free', 'Clutches, free', 'Ovoviviparity', 'Clutches, terrestrial', 'Clutches, in vegetation'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Reproduction", face = "bold") +
  ggtitle("Reproduction (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p9

# tachresist
nice_cols <- brewer.pal(5, "Pastel1")
p10 <- ggplot(tachresist_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c('Cocoons', 'Diapause or dormancy', 'Eggs, gemmula, statoblasts', 'Housings against desiccation', 'None'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Resistance forms", face = "bold") +
  ggtitle("Resistance forms (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p10

# tachrespir
nice_cols <- brewer.pal(5, "Pastel1")
p11 <- ggplot(tachrespir_L, aes(fill = trait, y= meanCWM , x = year)) + 
  geom_flow(aes(alluvium = trait), alpha= .5, color = "white",
            curve_type = "linear", 
            width = 0.8) +
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_fill_manual(labels=c('Gill', 'Plastron', 'Spiracle (aerial)', 'Tegument', 'Hydrostatic vesicle (aerial)'), values = nice_cols) +
  theme_minimal_grid() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.955, .925),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        legend.box.margin = margin(0, 5, 5, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8)) +
  labs(fill = "Respiration", face = "bold") +
  ggtitle("Respiration (Tachet et al. 2010)") +
  ylab("Proportion of community") +
  xlab("Year")
p11

# combined plots
(CWMs_combined <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, 
                                     align = "hv", axis = "bt", ncol = 3, nrow = 4))

# save plots
tiff(filename = "Plots/LT_Change_in_traits_check.tiff", width = 20, height = 20, units = 'in', res = 600, compression = 'lzw')
CWMs_combined
dev.off()

svg(filename = "Plots/LT_Change_in_traits.svg", width = 20, height = 20)
CWMs_combined
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