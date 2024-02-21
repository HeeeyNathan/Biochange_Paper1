# load packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggbreak)
library(patchwork)
library(ggrepel)
library(broom)
library(vegan)

# Read in data
df_short <- read_csv("Data/LT_annual_taxalist_short_wGroup.csv", show_col_types = FALSE)

# Initial processing
df_long <- df_short %>%
  pivot_longer(cols = c(`2010`, `2011`, `2012`, `2013`, `2014`, `2015`, `2016`, `2017`, `2019`, `2020`),
               names_to = "year",
               values_to = "abundance") %>%
  mutate(year = as.numeric(year), # Convert the "year" column to numeric directly here
         relative_abundance = round(abundance / sum(abundance), 2)) %>%
  filter(!year %in% c(2011, 2019)) %>% # Combine year filters into one step
  group_by(year) %>%
  mutate(relative_abundance = round(abundance / sum(abundance), 2), # Ensure two decimal places
         relative_abundance_perc = round(relative_abundance * 100, 2), # Convert to percentage and round
         log_abundance = log(abundance + 1)) %>%
  ungroup() %>%
  mutate(
    group2 = case_when(
      group %in% c("Ephemeroptera", "Plecoptera", "Trichoptera") ~ "EPT",
      group %in% c("Diptera", "Coleoptera", "Odonata", "Hemiptera", "Lepidoptera", "Megaloptera") ~ "Insect",
      group %in% c("Gastropoda", "Bivalvia") ~ "Mollusc",
      group %in% c("Oligochaeta", "Hirudinea") ~ "Annelid",
      group %in% c("Isopoda", "Decapoda", "Amphipoda", "Mysida") ~ "Crustacea",
      TRUE ~ group
    )
  )

attach(df_long)

# ephemeroptera
df_long_ephemeroptera <- filter(df_long, group == "Ephemeroptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_ephemeroptera <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
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
  geom_line(linewidth = 0.35, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  # scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  scale_colour_gradient2(low="#f0bb00ff", mid = "white", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_ephemeroptera_filtered$year), max(df_long_ephemeroptera_filtered$year), by = 1),
                     labels = seq(min(df_long_ephemeroptera_filtered$year), max(df_long_ephemeroptera_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Ephemeropteran abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# plecoptera
df_long_plecoptera <- filter(df_long, group == "Plecoptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_plecoptera <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
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
  geom_line(linewidth = 0.35, linetype = "dotted") + # no significant changes in plecopterans were found
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  # scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  scale_colour_gradient2(low="#f0bb00ff", mid = "white", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1),
                     labels = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Plecopteran abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# trichoptera
df_long_trichoptera <- filter(df_long, group == "Trichoptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_trichoptera <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
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
  geom_line(linewidth = 0.35, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  # scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)")) +
  scale_colour_gradient2(low="#f0bb00ff", mid = "white", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)")) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1),
                     labels = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Trichopteran abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# insect
df_long_insect <- filter(df_long, group2 == "Insect")
# Filter the dataset, only keeping species present in at least 2 years
df_long_insect <- df_long_insect %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
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
  geom_line(linewidth = 0.35, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  # scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  scale_colour_gradient2(low="#f0bb00ff", mid = "white", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)")) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1),
                     labels = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Insect abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# crustacea
df_long_crustacea <- filter(df_long, group2 == "Crustacea")
# Filter the dataset, only keeping species present in at least 2 years
df_long_crustacea <- df_long_crustacea %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
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
  geom_line(linewidth = 0.35, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  # scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  scale_colour_gradient2(low="#f0bb00ff", mid = "white", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1),
                     labels = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Crustacean abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# mollusc
df_long_mollusc <- filter(df_long, group2 == "Mollusc")
# Filter the dataset, only keeping species present in at least 2 years
df_long_mollusc <- df_long_mollusc %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
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
  geom_line(linewidth = 0.35, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  # scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  scale_colour_gradient2(low="#f0bb00ff", mid = "white", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1),
                     labels = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Mollusc abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# annelid
df_long_annelid <- filter(df_long, group2 == "Annelid")
# Filter the dataset, only keeping species present in at least 2 years
df_long_annelid <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
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
  geom_line(linewidth = 0.35, linetype = "dotted") + # no significant changes in annelids were observed
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  geom_point(size = 2, aes(shape = group)) +
  geom_text_repel(data = first_years, aes(label = taxonname), hjust = -0.1, vjust = 0, nudge_y = 0.01, size = 3) +
  geom_text_repel(data = last_years, aes(label = taxonname), hjust = 1.1, vjust = 0, nudge_y = 0.075, size = 3) +
  labs(x = "Year", y = "log(Abundance + 1)", color = "Taxonomic Group", shape = "Taxonomic group", linetype = "P value") +
  # scale_colour_gradient(low="#f0bb00ff", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  scale_colour_gradient2(low="#f0bb00ff", mid = "white", high="#00ccbaff", aes("lm(log(abundance + 1) ~ year)"))+
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1),
                     labels = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(legend.position = "right", panel.border = element_rect(color = "grey80", fill = NA)) +
  ggtitle("Change in Annelid abundance through time") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# create combined table of all winners and losers
winners_losers <- bind_rows(df_long_ephemeroptera_filtered,
                            df_long_plecoptera_filtered,
                            df_long_trichoptera_filtered,
                            df_long_insect_filtered,
                            df_long_mollusc_filtered,
                            df_long_annelid_filtered,
                            df_long_crustacea_filtered)
# make better looking table
# removing duplicate rows (some species are covered by two or more datasets e.g., ephemeroptera and insect)
winners_losers <- distinct(winners_losers)
# remove unnecessary columns
winners_losers <- winners_losers %>%
  select(-year, -abundance, -relative_abundance, -relative_abundance_perc, -log_abundance, -term)
# remove species duplicates (the information is the same in every line)
winners_losers <- distinct(winners_losers, taxonname, .keep_all = TRUE)
# save the final table
write.csv(winners_losers, "Outputs/LT_2010-2020_winners_losers.csv", row.names = F)

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
