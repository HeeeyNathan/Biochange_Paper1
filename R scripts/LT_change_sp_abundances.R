# load packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggbreak)
library(patchwork)
library(ggrepel)
library(broom)

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
                         "Oligochaeta" = "Annelid", "Hirudinea" = "Annelid", "Megaloptera" = "Other",
                         "Isopoda" = "Crustacea", "Decapoda" = "Crustacea", "Amphipoda" = "Crustacea", "Mysida" = "Crustacea"))

df_long <- df_long %>%
  mutate(group3 = recode(group, "Diptera" = "Insect", "Ephemeroptera" = "Insect", "Coleoptera" = "Insect", 
                         "Odonata" = "Insect", "Trichoptera" = "Insect", "Hemiptera" = "Insect", "Plecoptera" = "Insect", 
                         "Lepidoptera" = "Insect", "Megaloptera" = "Insect",
                         "Gastropoda" = "Mollusc", "Bivalvia" = "Mollusc", 
                         "Oligochaeta" = "Annelid", "Hirudinea" = "Annelid",
                         "Isopoda" = "Crustacea", "Decapoda" = "Crustacea", "Amphipoda" = "Crustacea", "Mysida" = "Crustacea"))

# Create a subset with species belonging to a specific taxonomic group (based on group2)
df_long_ept <- df_long %>% filter(group2 == "EPT")
df_long_ephemeroptera <- df_long %>% filter(group == "Ephemeroptera")
df_long_plecoptera <- df_long %>% filter(group == "Plecoptera")
df_long_trichoptera <- df_long %>% filter(group == "Trichoptera")
df_long_diptera <- df_long %>% filter(group2 == "Diptera")
df_long_insect <- df_long %>% filter(group3 == "Insect")
df_long_mollusc <- df_long %>% filter(group2 == "Mollusc")
df_long_annelid <- df_long %>% filter(group2 == "Annelid")
df_long_crustacea <- df_long %>% filter(group2 == "Crustacea")
df_long_other <- df_long %>% filter(group2 == "Other")

# ephemeroptera
# Filter the dataset, only keeping species present in at least 2 years
df_long_ephemeroptera <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_ephemeroptera <- left_join(df_long_ephemeroptera, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_ephemeroptera_filtered <- df_long_ephemeroptera %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_ephemeroptera_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_ephemeroptera_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_ephemeroptera_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_ephemeroptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_ephemeroptera_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_ephemeroptera_filtered$year), max(df_long_ephemeroptera_filtered$year), by = 1),
                     labels = seq(min(df_long_ephemeroptera_filtered$year), max(df_long_ephemeroptera_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Ephemeropteran abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# plecoptera
# Filter the dataset, only keeping species present in at least 2 years
df_long_plecoptera <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_plecoptera <- left_join(df_long_plecoptera, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_plecoptera_filtered <- df_long_plecoptera %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_plecoptera_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_plecoptera_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_plecoptera_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_plecoptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_plecoptera_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1),
                     labels = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Plecopteran abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# trichoptera
# Filter the dataset, only keeping species present in at least 2 years
df_long_trichoptera <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_trichoptera <- left_join(df_long_trichoptera, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_trichoptera_filtered <- df_long_trichoptera %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_trichoptera_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_trichoptera_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_trichoptera_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_trichoptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_trichoptera_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1),
                     labels = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Trichopteran abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# diptera
# Filter the dataset, only keeping species present in at least 2 years
df_long_diptera <- df_long_diptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_diptera %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_diptera <- left_join(df_long_diptera, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_diptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_diptera %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_diptera_filtered <- df_long_diptera %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_diptera_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_diptera_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_diptera_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_diptera.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_diptera_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_diptera_filtered$year), max(df_long_diptera_filtered$year), by = 1),
                     labels = seq(min(df_long_diptera_filtered$year), max(df_long_diptera_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Dipteran abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# insect
# Filter the dataset, only keeping species present in at least 2 years
df_long_insect <- df_long_insect %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_insect %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_insect <- left_join(df_long_insect, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_insect %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_insect %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_insect_filtered <- df_long_insect %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_insect_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_insect_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_insect_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_insect.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_insect_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1),
                     labels = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Insect abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# crustacea
# Filter the dataset, only keeping species present in at least 2 years
df_long_crustacea <- df_long_crustacea %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_crustacea %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_crustacea <- left_join(df_long_crustacea, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_crustacea %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_crustacea %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_crustacea_filtered <- df_long_crustacea %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_crustacea_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_crustacea_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_crustacea_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_crustacea.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_crustacea_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1),
                     labels = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Crustacean abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# mollusc
# Filter the dataset, only keeping species present in at least 2 years
df_long_mollusc <- df_long_mollusc %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_mollusc %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_mollusc <- left_join(df_long_mollusc, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_mollusc %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_mollusc %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_mollusc_filtered <- df_long_mollusc %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_mollusc_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_mollusc_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_mollusc_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_mollusc.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_mollusc_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1),
                     labels = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Mollusc abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# annelid
# Filter the dataset, only keeping species present in at least 2 years
df_long_annelid <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_annelid %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_annelid <- left_join(df_long_annelid, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_annelid %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_annelid %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_annelid_filtered <- df_long_annelid %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_annelid_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_annelid_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_annelid_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_annelid.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_annelid_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1),
                     labels = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Annelid abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# other
# Filter the dataset, only keeping species present in at least 2 years
df_long_other <- df_long_other %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 2) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_other %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_other <- left_join(df_long_other, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_other %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

bottom_species <- df_long_other %>%
  group_by(taxonname) %>%
  summarise(max_estimate = max(estimate)) %>%
  top_n(5, wt = -max_estimate) %>%
  arrange(desc(max_estimate)) %>% 
  ungroup()

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_other_filtered <- df_long_other %>%
  semi_join(selected_species, by = "taxonname")

unique(df_long_other_filtered$taxonname)

# Determine the taxonnames to display for the first and last year
first_years <- df_long_other_filtered %>%
  group_by(taxonname) %>%
  filter(year == min(year)) %>%
  ungroup()

last_years <- df_long_other_filtered %>%
  group_by(taxonname) %>%
  filter(n() > 1) %>%
  filter(year == max(year)) %>%
  ungroup()

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_other.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_other_filtered, aes(x = year, y = log_abundance, color = estimate, group = taxonname)) +
  geom_line(linewidth = 0.35, linetype = "dotted") +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_other_filtered$year), max(df_long_other_filtered$year), by = 1),
                     labels = seq(min(df_long_other_filtered$year), max(df_long_other_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in other invertebrate abundances through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
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
