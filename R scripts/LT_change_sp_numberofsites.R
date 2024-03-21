# load packages
library(tidyverse)
library(broom)

# Read in data
df_long <- read_csv("Data/LT_taxalist_2010-2020_long.csv", show_col_types = FALSE)

# Initial processing
df_long_summ <- df_long %>%
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
  mutate(log_site_count = log(site_count + 1))# Add the log(site_count + 1) column

# ephemeroptera
df_long_ephemeroptera <- filter(df_long_summ, group == "Ephemeroptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_ephemeroptera <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_ephemeroptera %>%
  group_by(taxonname) %>%
  do(tidy(lm(log(site_count + 1) ~ year, data = .))) %>%
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
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_ephemeroptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_ephemeroptera_filtered <- df_long_ephemeroptera %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_ephemeroptera_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_ephemeroptera_filtered, aes(x = year, y = log_site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(site count + 1)",
    title = "Change in ephemeroptera distribution through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling.\nSpecies with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time."
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
df_long_plecoptera <- filter(df_long_summ, group == "Plecoptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_plecoptera <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_plecoptera %>%
  group_by(taxonname) %>%
  do(tidy(lm(log(site_count + 1) ~ year, data = .))) %>%
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
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_plecoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_plecoptera_filtered <- df_long_plecoptera %>%
  semi_join(selected_species, by = "taxonname")

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Winner_Loser_plecoptera_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
ggplot(df_long_plecoptera_filtered, aes(x = year, y = log_site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate), linetype = "dashed") +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(site count + 1)",
    title = "Change in plecoptera distribution through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling.\nSpecies with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time."
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
  scale_x_continuous(breaks = seq(min(df_long_plecoptera_filtered$year), max(df_long_plecoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
dev.off()

# trichoptera
df_long_trichoptera <- filter(df_long_summ, group == "Trichoptera")
# Filter the dataset, only keeping species present in at least 2 years
df_long_trichoptera <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_trichoptera %>%
  group_by(taxonname) %>%
  do(tidy(lm(log(site_count + 1) ~ year, data = .))) %>%
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
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_trichoptera %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_trichoptera_filtered <- df_long_trichoptera %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_trichoptera_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Distribution_trichoptera.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_trichoptera_filtered, aes(x = year, y = log_site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(site count + 1)",
    title = "Change in trichoptera distribution through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling.\nSpecies with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time."
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
  scale_x_continuous(breaks = seq(min(df_long_trichoptera_filtered$year), max(df_long_trichoptera_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
dev.off()

# insect
df_long_insect <- filter(df_long_summ, group2 == "Non-EPT Insect")
# Filter the dataset, only keeping species present in at least 2 years
df_long_insect <- df_long_insect %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_insect %>%
  group_by(taxonname) %>%
  do(tidy(lm(log(site_count + 1) ~ year, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup()

df_long_insect <- left_join(df_long_insect, model_outputs, by = "taxonname")

# Select the five unique species with the highest amount of change from the first to the last year
top_species <- df_long_insect %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate > 0) %>%  # Keep only rows with positive estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp <= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_insect %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  # mutate(bmwp_priority = ifelse(bmwp <= 5, 1, 2)) %>% #assign priority based on bmwp score
  arrange(priority, desc(estimate)) %>%  # Arrange by priority, then by descending estimates
  slice_head(n = 6) %>%  # Select the top 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_insect_filtered <- df_long_insect %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_insect_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Distribution_insect.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_insect_filtered, aes(x = year, y = log_site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(site count + 1)",
    title = "Change in non-EPT insect distribution through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling.\nSpecies with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time."
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
  scale_x_continuous(breaks = seq(min(df_long_insect_filtered$year), max(df_long_insect_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
dev.off()

# crustacea
df_long_crustacea <- filter(df_long_summ, group2 == "Crustacea")
# Filter the dataset, only keeping species present in at least 2 years
df_long_crustacea <- df_long_crustacea %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_crustacea %>%
  group_by(taxonname) %>%
  do(tidy(lm(log(site_count + 1) ~ year, data = .))) %>%
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
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_crustacea %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_crustacea_filtered <- df_long_crustacea %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_crustacea_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Distribution_crustacea.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_crustacea_filtered, aes(x = year, y = log_site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(site count + 1)",
    title = "Change in crustacea distribution through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling.\nSpecies with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time."
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
  scale_x_continuous(breaks = seq(min(df_long_crustacea_filtered$year), max(df_long_crustacea_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
dev.off()

# mollusc
df_long_mollusc <- filter(df_long_summ, group2 == "Mollusc")
# Filter the dataset, only keeping species present in at least 2 years
df_long_mollusc <- df_long_mollusc %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_mollusc %>%
  group_by(taxonname) %>%
  do(tidy(lm(log(site_count + 1) ~ year, data = .))) %>%
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
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_mollusc %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_mollusc_filtered <- df_long_mollusc %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_mollusc_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Distribution_mollusc.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_mollusc_filtered, aes(x = year, y = log_site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(site count + 1)",
    title = "Change in mollusc distribution through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling.\nSpecies with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time."
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
  scale_x_continuous(breaks = seq(min(df_long_mollusc_filtered$year), max(df_long_mollusc_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
dev.off()

# annelid
df_long_annelid <- filter(df_long_summ, group2 == "Annelid")
# Filter the dataset, only keeping species present in at least 2 years
df_long_annelid <- df_long_annelid %>%
  group_by(taxonname) %>%
  filter(n_distinct(year) >= 6) %>%
  ungroup()

# Group the data by species and fit linear regression model for each group
model_outputs <- df_long_annelid %>%
  group_by(taxonname) %>%
  do(tidy(lm(log(site_count + 1) ~ year, data = .))) %>%
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
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, and p-values

bottom_species <- df_long_annelid %>%
  distinct(taxonname, .keep_all = TRUE) %>%  # Ensure unique species, keeping all data for now
  filter(estimate < 0) %>%  # Keep only rows with negative estimates
  mutate(priority = ifelse(p.value <= 0.05, 1, 2)) %>%  # Assign priority based on p-value
  arrange(priority, estimate) %>%  # Arrange by priority, then by ascending estimates
  slice_head(n = 6) %>%  # Select the bottom 6 species after arrangement
  select(taxonname, estimate, p.value, priority)  # Keep only species names, estimates, p-values, and priority

# Join the top and bottom species datasets
selected_species <- bind_rows(top_species, bottom_species)

# Filter df_long to include only the rows for the selected species
df_long_annelid_filtered <- df_long_annelid %>%
  semi_join(selected_species, by = "taxonname")
unique(df_long_annelid_filtered$taxonname)

# Plot abundances by year with line colors grouped by taxonomic order
tiff(filename = "Plots/Distribution_annelid.tiff", width = 20, height = 15, units = 'in', res = 600, compression = 'lzw')
ggplot(df_long_annelid_filtered, aes(x = year, y = log_site_count)) +
  geom_line(aes(color = estimate), linewidth = 0.5) +
  geom_point(shape = 16, aes(color = estimate), size = 3) +
  geom_smooth(method = "lm", se = FALSE, aes(color = estimate, linetype = ifelse(p.value <= 0.05, "p < 0.05", "p > 0.05"))) +
  scale_color_gradient(low = "#f0bb00ff", high = "#00ccbaff") +
  facet_wrap(~taxonname, scales = "free_y") +
  labs(
    x = "Year",
    y = "Log(site count + 1)",
    title = "Change in annelid distribution through time by species",
    caption = "Analysis includes species present in at least 6 years of sampling.\nSpecies with strongest estimates (positive or negative) plotted. Priority given to species significantly (p \u2264 0.05) changing through time."
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
  scale_x_continuous(breaks = seq(min(df_long_annelid_filtered$year), max(df_long_annelid_filtered$year), by = 1)) +
  guides(color = guide_legend(title = "Estimate"), linetype = guide_legend(title = "P Value")) +
  geom_text(aes(label = paste("BMWP:", bmwp), x = Inf, y = Inf), position = position_nudge(y = -0.5), hjust = 1.1, vjust = 2, check_overlap = TRUE, size = 2, fontface = "bold")
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
tiff(filename = "Plots/Winner_Loser_distribution.tiff", width = 10, height = 8, units = 'in', res = 300, compression = 'lzw')
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
