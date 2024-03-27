# Load necessary library
library(ggplot2)
library(dplyr)

# attach data
response_gls <- readRDS("Outputs/glsTrends_site_level.rds")
site_metadata <- read.csv("Data/LT_site_metadata.csv", h = T, sep = ",")

# remove unnecessary columns
colnames(site_metadata)
selected_columns <- c("site_id", "River_type", "Heavily_modified", "River")
site_metadata_clean <- site_metadata %>%
  select(all_of(selected_columns))

# merge dataframes
merged_df <- merge(response_gls, site_metadata_clean, by = "site_id", all.x = TRUE)

# Create a new column to categorize estimates as below or above zero
merged_df$EstimateCategory <- ifelse(merged_df$estimate < 0, "Negative trends", "Positive trends")

# Define the desired responses
taxo_responses <- c("abundance", "spp_richness", "E10", "shannonsH", "turnover")
func_responses <- c("FRed", "FRic", "FEve", "FDis", "F_turnover")
group_responses <- c("ept_spp_richness", "ept_abundance",
                     "insect_spp_richness", "insect_abundance", 
                     "crustacea_spp_richness", "crustacea_abundance", 
                     "mollusc_spp_richness", "mollusc_abundance",
                     "annelid_spp_richness", "annelid_abundance")

# RIVER TYPE
# Function to calculate counts and percentages for each river type
calculate_RivTyp <- function(df, response_type) {
  df %>%
    filter(Response == response_type, !is.na(estimate)) %>%
    group_by(River_type, EstimateCategory, Response) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    mutate(percentage = count / sum(count) * 100)
}

# Function to create a plot for a specific response type
plot_response_RivTyp <- function(data, title) {
  ggplot(data, aes(x = River_type, y = count, fill = EstimateCategory)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = title,
         x = "River Type", y = "Number of sites",
         fill = "GLS Estimates") +
    scale_fill_manual(values = c("Positive trends" = "#95ccba", "Negative trends" = "#f2cc84")) +
    theme_minimal() +
    facet_wrap(~Response, scales = "free_y") +
    guides(fill = guide_legend(title = "Estimate"))
}

# Create a list to store the dataframes
taxo_list <- lapply(taxo_responses, function(response) calculate_RivTyp(merged_df, response))
func_list <- lapply(func_responses, function(response) calculate_RivTyp(merged_df, response))
group_list <- lapply(group_responses, function(response) calculate_RivTyp(merged_df, response))

# Combine the dataframes
taxo_RivTyp <- bind_rows(taxo_list)
func_RivTyp <- bind_rows(func_list)
group_RivTyp <- bind_rows(group_list)

# Define the desired order of responses
taxo_order <- factor(taxo_responses, levels = taxo_responses, ordered = TRUE)
func_order <- factor(func_responses, levels = func_responses, ordered = TRUE)
group_order <- factor(group_responses, levels = group_responses, ordered = TRUE)

# Convert 'Response' to an ordered factor
taxo_RivTyp$Response <- factor(taxo_RivTyp$Response, levels = taxo_order)
func_RivTyp$Response <- factor(func_RivTyp$Response, levels = func_order)
group_RivTyp$Response <- factor(group_RivTyp$Response, levels = group_order)

# Plot each response type
tiff("Plots/Taxonomic responses by river type.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
plot_response_RivTyp(taxo_RivTyp, "Taxonomic responses by river type")
dev.off()
tiff("Plots/Functional responses by river type.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
plot_response_RivTyp(func_RivTyp, "Functional responses by river type")
dev.off()
tiff("Plots/Taxonomic group responses by river type.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
plot_response_RivTyp(group_RivTyp, "Group responses by river type")
dev.off()

# EXTENT OF MODIFICATION
# Function to calculate counts and percentages for each river type
calculate_Mod <- function(df, response_type) {
  df %>%
    filter(Response == response_type, !is.na(estimate)) %>%
    group_by(Heavily_modified, EstimateCategory, Response) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    mutate(percentage = count / sum(count) * 100)
}

# Function to create a plot for a specific response type
plot_response_Mod <- function(data, title) {
  ggplot(data, aes(x = Heavily_modified, y = count, fill = EstimateCategory)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = title,
         x = "Heavily modified", y = "Number of sites",
         fill = "GLS Estimates") +
    scale_fill_manual(values = c("Positive trends" = "#95ccba", "Negative trends" = "#f2cc84")) +
    theme_minimal() +
    facet_wrap(~Response, scales = "free_y") +
    guides(fill = guide_legend(title = "Estimate"))
}

# Create a list to store the dataframes
taxo_list <- lapply(taxo_responses, function(response) calculate_Mod(merged_df, response))
func_list <- lapply(func_responses, function(response) calculate_Mod(merged_df, response))
group_list <- lapply(group_responses, function(response) calculate_Mod(merged_df, response))

# Combine the dataframes
taxo_Mod <- bind_rows(taxo_list)
func_Mod <- bind_rows(func_list)
group_Mod <- bind_rows(group_list)

# Define the desired order of responses
taxo_order <- factor(taxo_responses, levels = taxo_responses, ordered = TRUE)
func_order <- factor(func_responses, levels = func_responses, ordered = TRUE)
group_order <- factor(group_responses, levels = group_responses, ordered = TRUE)

# Convert 'Response' to an ordered factor
taxo_Mod$Response <- factor(taxo_Mod$Response, levels = taxo_order)
func_Mod$Response <- factor(func_Mod$Response, levels = func_order)
group_Mod$Response <- factor(group_Mod$Response, levels = group_order)

# Plot each response type
tiff("Plots/Taxonomic responses by extent of modification.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
plot_response_Mod(taxo_Mod, "Taxonomic responses by extent of modification")
dev.off()
tiff("Plots/Functional responses by extent of modificatione.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
plot_response_Mod(func_Mod, "Functional responses by extent of modification")
dev.off()
tiff("Plots/Taxonomic group responses by extent of modification.tiff", width = 10, height = 10, units = 'in', res = 600, compression = 'lzw')
plot_response_Mod(group_Mod, "Taxonomic group responses by extent of modification")
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