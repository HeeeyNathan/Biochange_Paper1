# load packages
library(tidyverse)
library(broom)

# Read in data
df_long <- read_csv("Data/LT_taxalist_2010-2020_long.csv", show_col_types = FALSE)

# Initial processing
df_long <- df_long %>%
  select(-day, -month, -taxon_code, -taxon_id, -ept, -epto, -coleoptera, -diptera, -insect, -mollusc, -annelid, -crustacea) %>%
  mutate(year = as.numeric(year)) %>%
  filter(!year %in% c(2011, 2019)) %>%
  mutate(
    group2 = case_when(
      group %in% c("Ephemeroptera", "Plecoptera", "Trichoptera") ~ "EPT",
      group %in% c("Diptera", "Coleoptera", "Odonata", "Hemiptera", "Lepidoptera", "Megaloptera") ~ "Non-EPT Insect",
      group %in% c("Gastropoda", "Bivalvia") ~ "Mollusc",
      group %in% c("Oligochaeta", "Hirudinea") ~ "Annelid",
      group %in% c("Isopoda", "Decapoda", "Amphipoda", "Mysida") ~ "Crustacea",
      TRUE ~ group),
    richness = if_else(abundance > 0, 1, 0)) %>%
  rename(taxonname = taxon_name) %>%
  group_by(taxonname, year, group, group2, bmwp) %>%
  summarise(
    site_count = n_distinct(site_id),
    abundance = sum(abundance),
    .groups = "drop") %>% # Explicitly drop the grouping
  mutate(log_site_count = log(site_count + 1),
         log_abundance = log(abundance + 1))

# Filter the dataset, only keeping species present in at least 6 years
df_long_filtered <- df_long %>%
  filter(abundance >= 1) %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Calculate trends for each species
model_outputs <- df_long_filtered %>%
  group_by(taxonname) %>%
  do(tidy(lm(site_count ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

# bind trends to main dataframe
df_long_filtered <- left_join(df_long_filtered, model_outputs, by = "taxonname")

# ephemeroptera
df_long_ephemeroptera <- filter(df_long_filtered, group == "Ephemeroptera")

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

# Filter df_long_filtered to include only the rows for the selected species
df_long_ephemeroptera_filtered <- df_long_ephemeroptera %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
# tiff(filename = "Plots/Winner_Loser_ephemeroptera_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_ephemeroptera_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in ephemeroptera distribution through time by species",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_ephemeroptera_filtered$year), max(df_long_ephemeroptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
# dev.off()

# plecoptera
df_long_plecoptera <- filter(df_long_filtered, group == "Plecoptera")

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

# Filter df_long_filtered to include only the rows for the selected species
df_long_plecoptera_filtered <- df_long_plecoptera %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
# tiff(filename = "Plots/Winner_Loser_plecoptera_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_plecoptera_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(shape = 16, aes(color = estimate), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate), linetype = "dashed") +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in plecoptera distribution through time by species",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
# dev.off()

# trichoptera
df_long_trichoptera <- filter(df_long_filtered, group == "Trichoptera")

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

# Filter df_long_filtered to include only the rows for the selected species
df_long_trichoptera_filtered <- df_long_trichoptera %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
# tiff(filename = "Plots/Winner_Loser_trichoptera_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_trichoptera_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(shape = 16, aes(color = estimate), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in trichoptera distribution through time by species",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
# dev.off()

# insect
df_long_insect <- filter(df_long_filtered, group2 == "Non-EPT Insect")

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

# Filter df_long_filtered to include only the rows for the selected species
df_long_insect_filtered <- df_long_insect %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
# tiff(filename = "Plots/Winner_Loser_insect_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_insect_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in non-EPT insect distribution through time by species",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
# dev.off()

# crustacea
df_long_crustacea <- filter(df_long_filtered, group2 == "Crustacea")

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

# Filter df_long_filtered to include only the rows for the selected species
df_long_crustacea_filtered <- df_long_crustacea %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
# tiff(filename = "Plots/Winner_Loser_crustacean_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_crustacea_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in crustacea distribution through time by species",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
# dev.off()

# mollusc
df_long_mollusc <- filter(df_long_filtered, group2 == "Mollusc")

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

# Filter df_long_filtered to include only the rows for the selected species
df_long_mollusc_filtered <- df_long_mollusc %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
# tiff(filename = "Plots/Winner_Loser_mollusc_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_mollusc_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in mollusc distribution through time by species",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
# dev.off()

# annelid
df_long_annelid <- filter(df_long_filtered, group2 == "Annelid")

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

# Filter df_long_filtered to include only the rows for the selected species
df_long_annelid_filtered <- df_long_annelid %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
# tiff(filename = "Plots/Winner_Loser_annelid_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_annelid_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_point(aes(color = estimate, shape = group), size = 3) +  # Use different shapes for groups
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in annedlid distribution through time by species",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
# dev.off()

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
  filter(taxonname %in% c("Cloeon_dipterum", "Heptagenia_sulphurea", 
                          "Leuctra_sp.", "Nemoura_sp.", 
                          "Oxyethira_sp.", "Hydropsyche_pellucidula", 
                          "Notonecta_glauca_ssp.", "Calopteryx_virgo",
                          "Erpobdella_octoculata", "Helobdella_stagnalis",
                          "Radix_auricularia", "Unio_sp."))

# control plotting order
ordering_list <- winners_losers_filtered %>%
  arrange(group2, taxonname) %>%
  pull(taxonname) %>%
  unique()

winners_losers_filtered$taxonname <- factor(winners_losers_filtered$taxonname, levels = ordering_list)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(winners_losers_filtered, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_smooth(method = "lm", se = F, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(aes(color = estimate, shape = group2), size = 3) +  # Use different shapes for groups
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in the distribution of select species through time",
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
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(winners_losers_filtered$year), max(winners_losers_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2.5, fontface = "bold")
dev.off()

# Create combined table of alien species
# Filter the dataset, only keeping species present in at least 4 years
df_long_alien <- df_long %>%
  filter(abundance >= 1) %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 4) %>%
  ungroup()

# Calculate trends for each species - site counts
model_outputs_alien_count <- df_long_alien %>%
  group_by(taxonname) %>%
  do(tidy(lm(site_count ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

# bind trends to main dataframe
df_long_alien_count <- left_join(df_long_alien, model_outputs_alien_count, by = "taxonname")

# filter the count dataset, only keeping invasive species
df_long_alien_count <- df_long_alien_count %>%
  filter(taxonname %in% c("Chaetogammarus_warpachowskyi", "Chelicorophium_curvispinum", 
                          "Limnomysis_benedeni", "Obesogammarus_crassus", 
                          "Paramysis_lacustris", "Pontogammarus_robustoides", 
                          "Dreissena_polymorpha", "Lithoglyphus_naticoides", 
                          "Orconectes_limosus"))

# control plotting order
ordering_list_count <- df_long_alien_count %>%
  arrange(group2, taxonname) %>%
  pull(taxonname) %>%
  unique()

df_long_alien_count$taxonname <- factor(df_long_alien_count$taxonname, levels = ordering_list_count)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Alien_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_alien_count, aes(x = year, y = site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +  # Color line by estimate value
  geom_smooth(method = "lm", se = F, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  geom_point(aes(color = estimate, shape = group2), size = 3) +  # Use different shapes for groups
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +  # Gradient color scale for lines
  facet_wrap(~taxonname, scales = "free_y") +  # Create a facet for each species, with free y scales
  labs(
    x = "Year",
    y = "Site count",
    title = "Change in the distribution of select alien species through time",
    caption = "Analysis only includes species present in at least 4 years of sampling.\nOnly alien species shown. Priority given to species significantly (p \u2264 0.05) changing through time."
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Slant year labels at 45 degrees
    panel.grid.major = element_blank(), # Remove major gridlines
    # panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.ticks = element_line(color = "gray90"), # Add tick marks
    panel.border = element_rect(colour = "gray90", fill=NA, linewidth=0.5)
  ) +
  scale_x_continuous(breaks = seq(min(df_long_alien_count$year), max(df_long_alien_count$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), shape = guide_legend(title = "Taxonomic group"),linetype = guide_legend(title = "P Value", override.aes = list(color = "black"))) +
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
