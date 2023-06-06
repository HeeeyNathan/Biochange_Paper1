# load packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggbreak)
library(patchwork)
library(ggrepel)

# read in data
df <- read.csv("Data/LT_annual_taxalist_short_wGroup.csv", header = T)

# Convert dataframe to long format
df_long <- df %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "year",
               values_to = "abundance")

# Convert the "Year" column to numeric
df_long$year <- as.numeric(gsub("X", "", df_long$year))

# Calculate relative abundance
df_long <- df_long %>%
  group_by(year) %>%
  mutate(relative_abundance = round(abundance / sum(abundance), 2))

# Remove year 2011 and 2019 from the dataframe (too few sampling sites)
df_long <- df_long %>%
  filter(year != "2011")

df_long <- df_long %>%
  filter(year != "2019")

# log transform abundance column
df_long$log_abundance <- log(df_long$abundance + 1)

# Create a new column with new groups, name it group2
df_long <- df_long %>%
  mutate(group2 = recode(group, "Ephemeroptera" = "EPT", "Plecoptera" = "EPT", "Trichoptera" = "EPT", 
                         "Lepidoptera" = "Other", "Gastropoda" = "Mollusc", "Odonata" = "Other", 
                         "Bivalvia" = "Mollusc", "Hemiptera" = "Other", "Platyhelminthes" = "Other",
                         "Oligochaeta" = "Annelid", "Hirudinea" = "Annelid", "Megaloptera" = "Other"))

df_long <- df_long %>%
  mutate(group3 = recode(group, "Diptera" = "Insect", "Ephemeroptera" = "Insect", "Coleoptera" = "Insect", 
                         "Odonata" = "Insect", "Trichoptera" = "Insect", "Hemiptera" = "Insect", "Plecoptera" = "Insect", 
                         "Lepidoptera" = "Insect", "Megaloptera" = "Insect",
                         "Gastropoda" = "Mollusc", "Bivalvia" = "Mollusc", 
                         "Oligochaeta" = "Annelid", "Hirudinea" = "Annelid"))


# Calculate the change in abundance and relative abundance
df_long <- df_long %>%
  group_by(taxonname) %>%
  arrange(year) %>%
  mutate(abund_change = round(log_abundance - lag(log_abundance, default = first(log_abundance)), 2)) %>%
  ungroup()

df_long <- df_long %>%
  group_by(taxonname) %>%
  arrange(year) %>%
  mutate(rel_abund_change = round(relative_abundance - lag(relative_abundance, default = first(relative_abundance)), 2)) %>%
  ungroup()

# Calculate the change from the first to last year for each taxon
df_long <- df_long %>%
  group_by(taxonname) %>%
  arrange(year) %>%
  mutate(abund_change_first_last = last(abundance) - first(abundance)) %>%
  ungroup()

df_long <- df_long %>%
  group_by(taxonname) %>%
  arrange(year) %>%
  mutate(log_abund_change_first_last = last(log_abundance) - first(log_abundance)) %>%
  ungroup()

# centre the change in abundance so that all values are positive and 0 is relative to 0 in the original vector used
df_long <- df_long %>%
  mutate(centred_abund_change_first_last = (abund_change_first_last-min(abund_change_first_last))/max(abund_change_first_last)) %>%
  ungroup()

df_long <- df_long %>%
  mutate(centred_log_abund_change_first_last = (log_abund_change_first_last-min(log_abund_change_first_last))/max(log_abund_change_first_last)) %>%
  ungroup()

# Create a subset with species belonging to a specific taxonomic group (based on group2)
df_long_ept <- df_long %>% filter(group2 == "EPT")
df_long_ephemeroptera <- df_long %>% filter(group == "Ephemeroptera")
df_long_plecoptera <- df_long %>% filter(group == "Plecoptera")
df_long_trichoptera <- df_long %>% filter(group == "Trichoptera")
df_long_diptera <- df_long %>% filter(group2 == "Diptera")
df_long_insect <- df_long %>% filter(group3 == "Insect")
df_long_mollusc <- df_long %>% filter(group2 == "Mollusc")
df_long_annelid <- df_long %>% filter(group2 == "Annelid")
df_long_crustacea <- df_long %>% filter(group == "Crustacea")

# ephemeroptera
# remove rare species (are likely not dominant species)
df_long_ephemeroptera <- df_long_ephemeroptera %>%
  filter(abundance >= 10)

# Filter the dataset, only keeping species present in at least 5 years
df_long_ephemeroptera <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>% # represents half of the time series
  ungroup()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_ephemeroptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_ephemeroptera, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_line(size = 0.05, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#0072B2", midpoint=median(df_long_ephemeroptera$centred_log_abund_change_first_last)) +
  # scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#009E73", midpoint=median(df_long_ephemeroptera$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "#0072B2", high = "#D55E00", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1),
                     labels = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# plecoptera
# remove rare species (are likely not dominant species)
df_long_plecoptera <- df_long_plecoptera %>%
  filter(abundance >= 10)

# Filter the dataset, only keeping species present in at least 5 years
df_long_plecoptera <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>% # represents half of the time series
  ungroup()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_plecoptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_plecoptera, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_line(size = 0.05, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#0072B2", midpoint=median(df_long_plecoptera$centred_log_abund_change_first_last)) +
  # scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#009E73", midpoint=median(df_long_plecoptera$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "#0072B2", high = "#D55E00", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1),
                     labels = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# trichoptera
# remove rare species (are likely not dominant species)
df_long_trichoptera <- df_long_trichoptera %>%
  filter(abundance >= 10)

# Filter the dataset, only keeping species present in at least 5 years
df_long_trichoptera <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>% # represents half of the time series
  ungroup()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_trichoptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_trichoptera, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_line(size = 0.05, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#0072B2", midpoint=median(df_long_trichoptera$centred_log_abund_change_first_last)) +
  # scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#009E73", midpoint=median(df_long_trichoptera$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "#0072B2", high = "#D55E00", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1),
                     labels = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# diptera
# remove rare species (are likely not dominant species)
df_long_diptera <- df_long_diptera %>%
  filter(abundance >= 10)

# Filter the dataset, only keeping species present in at least 5 years
df_long_diptera <- df_long_diptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>% # represents half of the time series
  ungroup()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_diptera %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_diptera %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_diptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_diptera, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_line(size = 0.05, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#0072B2", midpoint=median(df_long_diptera$centred_log_abund_change_first_last)) +
  # scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#009E73", midpoint=median(df_long_diptera$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "#0072B2", high = "#D55E00", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1),
                     labels = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# mollusc
# remove rare species (are likely not dominant species)
df_long_mollusc <- df_long_mollusc %>%
  filter(abundance >= 10)

# Filter the dataset, only keeping species present in at least 5 years
df_long_mollusc <- df_long_mollusc %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>% # represents half of the time series
  ungroup()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_mollusc %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_mollusc %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_mollusc.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_mollusc, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_line(size = 0.05, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#0072B2", midpoint=median(df_long_mollusc$centred_log_abund_change_first_last)) +
  # scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#009E73", midpoint=median(df_long_mollusc$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "#0072B2", high = "#D55E00", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1),
                     labels = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# annelid
# remove rare species (are likely not dominant species)
df_long_annelid <- df_long_annelid %>%
  filter(abundance >= 10)

# Filter the dataset, only keeping species present in at least 5 years
df_long_annelid <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>% # represents half of the time series
  ungroup()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_annelid.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_annelid, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_line(size = 0.05, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#0072B2", midpoint=median(df_long_annelid$centred_log_abund_change_first_last)) +
  # scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#009E73", midpoint=median(df_long_annelid$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "#0072B2", high = "#D55E00", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1),
                     labels = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# crustacea
# remove rare species (are likely not dominant species)
df_long_crustacea <- df_long_crustacea %>%
  filter(abundance >= 10)

# Filter the dataset, only keeping species present in at least 5 years
df_long_crustacea <- df_long_crustacea %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>% # represents half of the time series
  ungroup()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_crustacea %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_crustacea %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_crustacea.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_crustacea, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_line(size = 0.05, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#0072B2", midpoint=median(df_long_crustacea$centred_log_abund_change_first_last)) +
  # scale_colour_gradient2(low="#D55E00", mid="#CC79A7", high="#009E73", midpoint=median(df_long_crustacea$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "#0072B2", high = "#D55E00", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_crustacea$year), max(df_long_crustacea$year), by = 1),
                     labels = seq(min(df_long_crustacea$year), max(df_long_crustacea$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# Determine the taxonnames to display for the first and last year
first_years <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_Annelid.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_annelid, aes(x = year, y = log_abundance, color = centred_log_abund_change_first_last, group = taxonname)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group") +
  scale_colour_gradient2(low="blue", mid="goldenrod", high="red", midpoint=median(df_long_annelid$centred_log_abund_change_first_last)) +
  # scale_color_gradient(low = "blue", high = "red", aes("Absolute change\nin log(abundance)\nfrom first to last year")) +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 15, 3, 2, 8, 7, 6, 5, 4, 1, 0, 18, 9, 10, 11, 12)) +
  scale_x_continuous(breaks = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1),
                     labels = seq(min(df_long_annelid$year), max(df_long_annelid$year), by = 1)) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# Plot data points
ggplot(df_long_annelid, aes(x = year, y = log_abundance)) +
  geom_point(size = 3, shape = 18) +
  geom_smooth(method = "lm", se = FALSE, aes(group = taxonname)) +
  annotate("text", x = last_points$year, y = last_points$log_abundance, 
           label = last_points$taxonname, hjust = -0.1, vjust = 0, color = "black", size = 3) +
  # annotate("text", x = summary_df$last_year, y = summary_df$log_abundance, 
  #          label = summary_df$taxonname, hjust = -0.1, vjust = 0, color = "black", size = 3) +
  # geom_text(data = first_last_years, aes(label = taxonname), hjust = 1, vjust = -1, nudge_x = 0.3) +
  labs(x = "Year", y = "log(Abundance + 1)") +
  scale_color_discrete() +
  theme_minimal() +
  theme(legend.position = "bottom")

p + geom_text(data = first_last_years, aes(x = max(year), y = log_abundance, label = taxonname),
              hjust = -0.2, vjust = 0, color = "black", size = 3)


# Plot abundances by year with line colors grouped by taxonomic order
ggplot(df_long, aes(x = year, y = relative_abundance, color = group2, group = taxonname)) +
  geom_line() +
  geom_text_repel(aes(label = taxonname), size = 2) +
  labs(x = "Year", y = "Abundance", color = "Taxonomic Group") +
  scale_color_discrete() +
  theme_minimal() +
  theme(legend.position = "bottom")

# Filter the dataframe to get the top 20 species names with the highest relative abundance in 2020
top_species_2020 <- df_long %>%
  filter(Year == "2020") %>%
  arrange(desc(Abundance)) %>%
  head(20)

top_species_2020_relative <- df_long %>%
  filter(Year == "2020") %>%
  arrange(desc(Relative_Abundance)) %>%
  head(20)

# Plot abundances by year with line colors representing the species group
ggplot(df_long, aes(x = Year, y = log(Abundance +1), color = group2, group = taxonname)) +
  geom_line() +
  # geom_text(data = top_species_2020, aes(label = taxonname), hjust = -0.1, vjust = -0.5, size = 3) +
  geom_text_repel(data = top_species_2020, aes(label = taxonname), size = 2) +
  labs(x = "Year", y = "Abundance", color = "Taxonomic Group") +
  scale_color_discrete() +
  theme_minimal() +
  theme(legend.position = "bottom", plot.margin = margin(10, 10, 10, 10, "pt"))

# Plot realtive abundances by year with line colors representing the species group
ggplot(df_long, aes(x = Year, y = Relative_Abundance, color = group2, group = taxonname)) +
  geom_line() +
  # geom_text(data = top_species_2020_relative, aes(label = taxonname), hjust = -0.1, vjust = -0.5, size = 3) +
  geom_text_repel(data = top_species_2020_relative, aes(label = taxonname), size = 2) +
  labs(x = "Year", y = "Abundance", color = "Taxonomic Group") +
  scale_color_discrete() +
  theme_minimal() +
  scale_x_continuous(labels = function(x) as.integer(x)) +
  # scale_y_cut(breaks = c(0.05, 0.1), which=c(1, 2, 3), scales= c(1, 2, 2), space = 0.25) +
  theme(legend.position = "bottom", plot.margin = margin(10, 10, 10, 10, "pt"))

# determining winners and losers
# Sort the species by change in relative abundance
df_change_sorted <- df_long %>%
  arrange(Change)

# Get the top 10 species with the largest changes
top_largest_changes <- df_change_sorted %>%
  tail(20) %>%
  pull(taxonname)

# Get the top 20 species with the smallest changes
top_smallest_changes <- df_change_sorted %>%
  head(20) %>%
  pull(taxonname)

# Combine the two lists of species
top_species <- c(top_largest_changes, top_smallest_changes)

# Filter the dataframe to include only the top 20 species
top_species_df <- df_long %>%
  filter(taxonname %in% top_species)

# Plot realtive abundances by year with line colors representing the species group
ggplot(df_long, aes(x = Year, y = log(Abundance + 1), color = group2, group = taxonname)) +
  geom_line() +
  # geom_text(data = top_species_2020_relative, aes(label = taxonname), hjust = -0.1, vjust = -0.5, size = 3) +
  geom_text_repel(data = top_species_df, aes(label = taxonname), size = 2) +
  labs(x = "Year", y = "Abundance", color = "Taxonomic Group") +
  scale_color_discrete() +
  theme_minimal() +
  scale_x_continuous(labels = function(x) as.integer(x)) +
  # scale_y_cut(breaks = c(0.05, 0.1), which=c(1, 2, 3), scales= c(1, 2, 2), space = 0.25) +
  theme(legend.position = "bottom", plot.margin = margin(10, 10, 10, 10, "pt"))

ggplot(top_species_df, aes(x = Year, y = Relative_Abundance, color = group2, group = taxonname)) +
  geom_line() +
  # geom_text(data = top_species_2020_relative, aes(label = taxonname), hjust = -0.1, vjust = -0.5, size = 3) +
  geom_text_repel(data = top_species_df, aes(label = taxonname), size = 2) +
  labs(x = "Year", y = "Abundance", color = "Taxonomic Group") +
  scale_color_discrete() +
  theme_minimal() +
  scale_x_continuous(labels = function(x) as.integer(x)) +
  # scale_y_cut(breaks = c(0.05, 0.1), which=c(1, 2, 3), scales= c(1, 2, 2), space = 0.25) +
  theme(legend.position = "bottom", plot.margin = margin(10, 10, 10, 10, "pt"))

# Sort the species by change in relative abundance
df_change_sorted <- df_long %>%
  arrange(Change)

# Get the species with the smallest changes
top_smallest_changes <- df_change_sorted %>%
  head(10) %>%
  pull(species)

# Filter the dataframe to include only the species with the smallest changes
df_smallest_changes <- df_long %>%
  filter(taxonname %in% top_smallest_changes)%>%
  distinct(taxonname, .keep_all = TRUE)

# Plot relative abundances for species with smallest changes
ggplot(df_long, aes(x = Year, y = Relative_Abundance, color = group2, group = taxonname)) +
  geom_line() +
  geom_text_repel(data = top_species_df, aes(label = taxonname), size = 3) +
  labs(x = "Year", y = "Relative Abundance") +
  scale_color_discrete() +
  theme_minimal() +
  scale_x_continuous(labels = function(x) as.integer(x)) +
  theme(legend.position = "bottom", plot.margin = margin(10, 30, 10, 10, "pt"))

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
