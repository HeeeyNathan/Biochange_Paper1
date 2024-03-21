# load packages
library(tidyverse)
library(broom)
library(tidyr)

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
  # mutate(cYear = year - mean(year, na.rm = TRUE)) %>%
  # mutate(iYear = year - min(year) + 1) %>%
  mutate(
    group2 = case_when(
      group %in% c("Ephemeroptera") ~ "Ephemeroptera",
      group %in% c("Plecoptera") ~ "Plecoptera",
      group %in% c("Trichoptera") ~ "Trichoptera",
      group %in% c("Diptera", "Coleoptera", "Odonata", "Hemiptera", "Lepidoptera", "Megaloptera") ~ "Non-EPT Insect",
      group %in% c("Gastropoda", "Bivalvia") ~ "Mollusc",
      group %in% c("Oligochaeta", "Hirudinea") ~ "Annelid",
      group %in% c("Isopoda", "Decapoda", "Amphipoda", "Mysida") ~ "Crustacea",
      TRUE ~ group))

# Filter the dataset, only keeping species present in at least 6 years
df_long <- df_long %>%
  filter(abundance >= 1) %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Calculate trends for each species
model_outputs <- df_long %>%
  group_by(taxonname) %>%
  do(tidy(lm(log_abundance ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()
# bind trends to main dataframe
df_long <- left_join(df_long, model_outputs, by = "taxonname")

# ephemeroptera
df_long_ephemeroptera <- filter(df_long, group == "Ephemeroptera")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_ephemeroptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp <= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_ephemeroptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp >= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_ephemeroptera_filtered <- df_long_ephemeroptera %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_ephemeroptera_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_ephemeroptera_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in ephemeroptera abundance through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling. BMWP scores are annotated on each plot.\n Species with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    # axis.line = element_line(color = "gray90"), # Add x and y axis lines
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_ephemeroptera_filtered$year), max(df_long_ephemeroptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
dev.off()

# plecoptera
df_long_plecoptera <- filter(df_long, group == "Plecoptera")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_plecoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp <= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

bottom_species <- df_long_plecoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp >= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_plecoptera_filtered <- df_long_plecoptera %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_plecoptera_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_plecoptera_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_plecoptera_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(shape = 16, aes(color = estimate), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate), linetype = "dashed") +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in plecoptera abundance through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling. BMWP scores are annotated on each plot.\n Species with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
dev.off()

# trichoptera
df_long_trichoptera <- filter(df_long, group == "Trichoptera")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_trichoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp <= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

bottom_species <- df_long_trichoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp >= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_trichoptera_filtered <- df_long_trichoptera %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_trichoptera_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_trichoptera_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_trichoptera_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(shape = 16, aes(color = estimate), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in trichoptera abundance through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling. BMWP scores are annotated on each plot.\n Species with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
dev.off()

# insect
df_long_insect <- filter(df_long, group2 == "Non-EPT Insect")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_insect %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp <= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

bottom_species <- df_long_insect %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp >= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, priority, estimate, p.value)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_insect_filtered <- df_long_insect %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_insect_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_insect_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_insect_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in insect abundance through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling. BMWP scores are annotated on each plot.\n Species with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
dev.off()

# crustacea
df_long_crustacea <- filter(df_long, group2 == "Crustacea")

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

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_crustacean_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_crustacea_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in crustacea abundance through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling. BMWP scores are annotated on each plot.\n Species with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
dev.off()

# mollusc
df_long_mollusc <- filter(df_long, group2 == "Mollusc")

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

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_mollusc_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_mollusc_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in mollusc abundance through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling. BMWP scores are annotated on each plot.\n Species with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
dev.off()

# annelid
df_long_annelid <- filter(df_long, group2 == "Annelid")

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

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_annelid_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_annelid_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate), linetype = "dashed") +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in annedlid abundance through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling. BMWP scores are annotated on each plot.\n Species with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
dev.off()

# create combined table of all winners and losers
winners_losers <- bind_rows(df_long_insect,
                            df_long_ephemeroptera,
                            df_long_plecoptera,
                            df_long_trichoptera,
                            df_long_insect,
                            df_long_crustacea,
                            df_long_mollusc,
                            df_long_annelid)

winners_losers_filtered <- winners_losers %>%
  filter(taxonname %in% c("Caenis_macrura", "Heptagenia_sulphurea", 
                          "Leuctra_sp.", "Taeniopteryx_nebulosa", 
                          "Neureclipsis_bimaculata", "Hydropsyche_siltalai", 
                          "Notonecta_glauca_ssp.", "Ischnura_elegans",
                          "Asellus_aquaticus", "Gammarus_lacustris",
                          "Radix_auricularia", "Unio_tumidus_ssp."))

# control plotting order
ordering_list <- winners_losers_filtered %>%
  arrange(group2, taxonname) %>%
  pull(taxonname) %>%
  unique()

winners_losers_filtered$taxonname <- factor(winners_losers_filtered$taxonname, levels = ordering_list)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_abundance.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(winners_losers_filtered, aes(x = year, y = log_abundance)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_smooth(method = "lm", se = F, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(aes(color = estimate, shape = group2), size = 3) +  # Use different shapes for groups
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Log(abundance + 1)",
    title = "Change in abundance through time by species",
    caption = "Analysis only includes species present in at least 6 years of sampling.\nOnly select species shown.Priority given to species significantly (p \u2264 0.05) changing through time."
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(winners_losers_filtered$year), max(winners_losers_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
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
