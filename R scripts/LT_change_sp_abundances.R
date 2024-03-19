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
      group %in% c("insect", "Plecoptera", "Trichoptera") ~ "EPT",
      group %in% c("Diptera", "Coleoptera", "Odonata", "Hemiptera", "Lepidoptera", "Megaloptera") ~ "Insect",
      group %in% c("Gastropoda", "Bivalvia") ~ "Mollusc",
      group %in% c("Oligochaeta", "Hirudinea") ~ "Annelid",
      group %in% c("Isopoda", "Decapoda", "Amphipoda", "Mysida") ~ "Crustacea",
      TRUE ~ group
    )
  )

# ephemeroptera
df_long_ephemeroptera <- filter(df_long, group == "Ephemeroptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_ephemeroptera <- df_long_ephemeroptera %>%
  filter(abundance >= 1) %>% # First, filter out taxa with an abundance of 0
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>% # Then, keep taxa present in at least 6 distinct years
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
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_ephemeroptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

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
tiff(filename = "Plots/Winner_Loser_ephemeroptera_new.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_ephemeroptera_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(shape = 16, aes(color = estimate), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year", 
    y = "Log(abundance + 1)", 
    title = "Change in ephemeropteran abundance through time by species",
    caption = "If available, only 12 species with the strongest estimates, negative or positive, are plotted.\nAnalysis only includes species present in at least 6 years of sampling."
  ) +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_long_ephemeroptera_filtered$year), max(df_long_ephemeroptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value"))
dev.off()

# plecoptera
df_long_plecoptera <- filter(df_long, group == "Plecoptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_plecoptera <- df_long_plecoptera %>%
  filter(abundance >= 1) %>% # First, filter out taxa with an abundance of 0
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>% # Then, keep taxa present in at least 6 distinct years (75% of time)
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
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_plecoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

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
tiff(filename = "Plots/Winner_Loser_plecoptera_new.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_plecoptera_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(shape = 16, aes(color = estimate), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate), linetype = "dashed") +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year", 
    y = "Log(abundance + 1)", 
    title = "Change in plecopteran abundance through time by species",
    caption = "If available, only 12 species with the strongest estimates, negative or positive, are plotted.\nAnalysis only includes species present in at least 6 years of sampling."
  ) +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value"))
dev.off()

# trichoptera
df_long_trichoptera <- filter(df_long, group == "Trichoptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_trichoptera <- df_long_trichoptera %>%
  filter(abundance >= 1) %>% # First, filter out taxa with an abundance of 0
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>% # Then, keep taxa present in at least 6 distinct years (75%)
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
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_trichoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

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
tiff(filename = "Plots/Winner_Loser_trichoptera_new.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_trichoptera_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(shape = 16, aes(color = estimate), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year", 
    y = "Log(abundance + 1)", 
    title = "Change in Trichopteran abundance through time by species",
    caption = "If available, only 12 species with the strongest estimates, negative or positive, are plotted.\nAnalysis only includes species present in at least 6 years of sampling."
  ) +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value"))
dev.off()

# insect
df_long_insect <- filter(df_long, group2 == "Insect")
# Filter the dataset, only keeping species present in at least 2 years
df_long_insect <- df_long_insect %>%
  filter(abundance >= 1) %>% # First, filter out taxa with an abundance of 0
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>% # Then, keep taxa present in at least 6 distinct years (75%)
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
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_insect %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

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
tiff(filename = "Plots/Winner_Loser_insect_new.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_insect_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(x = "Year", y = "Log(abundance + 1)", 
    title = "Change in Insect abundance through time by species",
    caption = "If available, only 12 species with the strongest estimates, negative or positive, are plotted.\nAnalysis only includes species present in at least 6 years of sampling.") +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value"))
dev.off()

# crustacea
df_long_crustacea <- filter(df_long, group2 == "Crustacea")
# Filter the dataset, only keeping species present in at least 2 years
df_long_crustacea <- df_long_crustacea %>%
  filter(abundance >= 1) %>% # First, filter out taxa with an abundance of 0
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>% # Then, keep taxa present in at least 6 distinct years (75%)
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
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_crustacea %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

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
tiff(filename = "Plots/Winner_Loser_crustacean_new.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_crustacea_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(x = "Year", y = "Log(abundance + 1)", 
       title = "Change in crustacea abundance through time by species",
       caption = "If available, only 12 species with the strongest estimates, negative or positive, are plotted.\nAnalysis only includes species present in at least 6 years of sampling.") +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value"))
dev.off()

# mollusc
df_long_mollusc <- filter(df_long, group2 == "Mollusc")
# Filter the dataset, only keeping species present in at least 2 years
df_long_mollusc <- df_long_mollusc %>%
  filter(abundance >= 1) %>% # First, filter out taxa with an abundance of 0
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>% # Then, keep taxa present in at least 2 distinct years
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
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_mollusc %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

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
tiff(filename = "Plots/Winner_Loser_mollusc_new.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_mollusc_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value < 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(x = "Year", y = "Log(abundance + 1)", 
       title = "Change in mollusc abundance through time by species",
       caption = "If available, only 12 species with the strongest estimates, negative or positive, are plotted.\nAnalysis only includes species present in at least 6 years of sampling.") +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value"))
dev.off()

# annelid
df_long_annelid <- filter(df_long, group2 == "Annelid")
# Filter the dataset, only keeping species present in at least 2 years
df_long_annelid <- df_long_annelid %>%
  filter(abundance >= 1) %>% # First, filter out taxa with an abundance of 0
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>% # Then, keep taxa present in at least 2 distinct years
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
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_annelid %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

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
tiff(filename = "Plots/Winner_Loser_annelid_new.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_annelid_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate), linetype = "dashed") +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(x = "Year", y = "Log(abundance + 1)", 
       title = "Change in annelid abundance through time by species",
       caption = "If available, only 12 species with the strongest estimates, negative or positive, are plotted.\nAnalysis only includes species present in at least 6 years of sampling.") +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value"))
dev.off()

# create combined table of all winners and losers
winners_losers <- bind_rows(df_long_insect_filtered,
                            df_long_plecoptera_filtered,
                            df_long_trichoptera_filtered,
                            df_long_insect_filtered,
                            df_long_crustacea_filtered,
                            df_long_mollusc_filtered,
                            df_long_annelid_filtered)
# make better looking table
# removing duplicate rows (some species are covered by two or more datasets e.g., insect and insect)
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
