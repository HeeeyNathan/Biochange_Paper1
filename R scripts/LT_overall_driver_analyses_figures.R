# Load required packages
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

# Read in and organise data
driver_gls <- readRDS("outputs/glsTrends_drivers_site_level.rds")
driver_gls <- driver_gls[driver_gls$driver != "(Intercept)",] # removes the intercept line from each response variable

# Change the how the driver names appear in the figures
driver_gls <- driver_gls %>%
  mutate(driver = case_when(
    driver == "sflow" ~ "Discharge",
    driver == "spH" ~ "pH",
    driver == "stemp" ~ "Temperature",
    # driver == "ssus_solid" ~ "Susp. solids",
    driver == "so2_dis" ~ "Diss. oxygen",
    # driver == "sBOD7" ~ "Biol. oxygen demand",
    driver == "sNH4.N" ~ "Ammonium",
    driver == "PC_axis1" ~ "Nutrients PCA",
    TRUE ~ driver  # Keep other values unchanged
  ))

# # Reorder the driver names so that they appear in ggplot the way we want
driver_gls$fDriver <- factor(driver_gls$driver, levels=c("Discharge", "pH", "Temperature", "Diss. oxygen", "Ammonium", "Nutrients PCA"), ordered = T) # reorder the driver names so that they appear in ggplot the way we want

# generate groups for plotting
unique(driver_gls$Response)
Group1 <- driver_gls[driver_gls$Response %in% c("abundance", "spp_richness", "E10", "shannonsH", "turnover"),]
Group2 <- driver_gls[driver_gls$Response %in% c("FRed", "FRic", "FEve", "FDis", "F_turnover"),]
Group3 <- driver_gls[driver_gls$Response %in% c("ept_spp_richness", "insect_spp_richness", "crustacea_spp_richness", "mollusc_spp_richness", "annelid_spp_richness"),]
Group4 <- driver_gls[driver_gls$Response %in% c("ept_abundance", "insect_abundance", "crustacea_abundance", "mollusc_abundance", "annelid_abundance"),]
Group5 <- driver_gls[driver_gls$Response %in% c("spp_rich_rare", "FRic.SES", "FEve.SES", "FDis.SES"),]

# Change response names to match other plots (consistency)
nrow(Group1)
fResponse <- data.frame(fResponse = factor(c(rep("Total abundance", 6), rep("Shannon evenness", 6), rep("Shannon diversity", 6), rep("Taxon richness", 6), rep("Taxon turnover", 6))))
Group1 <- cbind(Group1, fResponse)
Group1$fResponse <- factor(Group1$fResponse, levels=c("Total abundance", "Taxon richness", "Shannon evenness", "Shannon diversity", "Taxon turnover"), ordered = T) # reorder the driver names so that they appear in ggplot the way we want

# Create plot
p1 <- ggplot(Group1, aes(x = Estimate, y = fDriver)) +
  scale_color_identity() +
  geom_point(
    shape = 16, size = 4, 
    aes(color = ifelse((0 >= `2.5 %` & 0 <= `97.5 %`) | (0 <= `2.5 %` & 0 >= `97.5 %`), "gray90", ifelse(Estimate >= 0, "#95ccba", "#f2cc84")))
  ) +
  geom_errorbar(
    width = 0, linewidth = 1.5,
    aes(
      xmax = (`2.5 %`), xmin = (`97.5 %`),
      color = ifelse((0 >= `2.5 %` & 0 <= `97.5 %`) | (0 <= `2.5 %` & 0 >= `97.5 %`), "gray90", ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
    )
  ) +
  facet_wrap(~ fResponse, nrow = 5, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"),
    legend.position = "none",
    text = element_text(size = 12),  # Adjust the size of all text elements
    axis.title = element_blank(),    # Remove axis titles
    axis.text = element_text(size = 12),  # Adjust the size of tick mark labels
    strip.text = element_text(size = 12, face = "bold")  # Adjust the size and face of facet labels
  ) +
  geom_vline(xintercept = 0, lty = 3) +
  ylab("") +
  xlab("")
p1

# svg(filename = "Plots/LT_OverallDrivers_TaxoIndices.svg", width = 6, height = 10, bg = "white")
# tiff(filename = "Plots/LT_OverallDrivers_TaxoIndices.tiff", width = 6, height = 10, units = 'in', res = 600, compression = 'lzw')
# p1
# dev.off()

# Change response names to match other plots (consistency)
nrow(Group2)
Group2$Response
fResponse2 <- data.frame(fResponse = factor(c(rep("Func. turnover", 6), rep("Func. dispersion", 6), rep("Func. evenness", 6), rep("Func. redundancy", 6), rep("Func. richness", 6))))
Group2 <- cbind(Group2, fResponse2)
Group2$fResponse <- factor(Group2$fResponse, levels=c("Func. redundancy", "Func. richness", "Func. evenness", "Func. dispersion", "Func. turnover"), ordered = T)

# Create plot
p2 <- ggplot(Group2, aes(x = Estimate, y = fDriver)) +
  scale_color_identity() +
  geom_point(
    shape = 16, size = 4, 
    aes(color = ifelse((0 >= `2.5 %` & 0 <= `97.5 %`) | (0 <= `2.5 %` & 0 >= `97.5 %`), "gray90", ifelse(Estimate >= 0, "#95ccba", "#f2cc84")))
  ) +
  geom_errorbar(
    width = 0, linewidth = 1.5,
    aes(
      xmax = (`2.5 %`), xmin = (`97.5 %`),
      color = ifelse((0 >= `2.5 %` & 0 <= `97.5 %`) | (0 <= `2.5 %` & 0 >= `97.5 %`), "gray90", ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
    )
  ) +
  facet_wrap(~ fResponse, nrow = 5, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"),
    legend.position = "none",
    text = element_text(size = 12),  # Adjust the size of all text elements
    axis.title = element_blank(),    # Remove axis titles
    axis.text = element_text(size = 12),  # Adjust the size of tick mark labels
    strip.text = element_text(size = 12, face = "bold")  # Adjust the size and face of facet labels
  ) +
  geom_vline(xintercept = 0, lty = 3) +
  ylab("") +
  xlab("")
p2

# svg(filename = "Plots/LT_OverallDrivers_funcIndices.svg", width = 6, height = 10, bg = "white")
# tiff(filename = "Plots/LT_OverallDrivers_FuncIndices.tiff", width = 6, height = 10, units = 'in', res = 600, compression = 'lzw')
# p2
# dev.off()

# Change response names to match other plots (consistency)
nrow(Group3)
Group3$Response
fResponse3 <- data.frame(fResponse = factor(c(rep("Annelid richness", 6), rep("Crustacea richness", 6), rep("EPT richness", 6), rep("Insect richness", 6), rep("Mollusc richness", 6))))
Group3 <- cbind(Group3, fResponse3)
Group3$fResponse <- factor(Group3$fResponse, levels=c("EPT richness", "Insect richness", "Crustacea richness", "Mollusc richness", "Annelid richness"), ordered = T)

# Create plot
p3 <- ggplot(Group3, aes(x = Estimate, y = fDriver)) +
  scale_color_identity() +
  geom_point(
    shape = 21, size = 4, fill = "white",
    aes(color = ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
  ) +
  geom_errorbar(
    width = 0, linewidth = 1,
    aes(xmin = (`2.5 %`), xmax = (`97.5 %`), color = ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
  ) +
  facet_wrap(~ fResponse, nrow = 5, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"),
    legend.position = "none",
    text = element_text(size = 12),  # Adjust the size of all text elements
    axis.title = element_blank(),    # Remove axis titles
    axis.text = element_text(size = 12),  # Adjust the size of tick mark labels
    strip.text = element_text(size = 12, face = "bold")  # Adjust the size and face of facet labels
  ) +
  geom_vline(xintercept = 0, lty = 3) +
  ylab("") +
  xlab("")
p3

# Change response names to match other plots (consistency)
nrow(Group4)
Group4$Response
fResponse4 <- data.frame(fResponse = factor(c(rep("Annelid abundance", 6), rep("Crustacea abundance", 6), rep("EPT abundance", 6), rep("Insect abundance", 6), rep("Mollusc abundance", 6))))
Group4 <- cbind(Group4, fResponse4)
Group4$fResponse <- factor(Group4$fResponse, levels=c("EPT abundance", "Insect abundance", "Crustacea abundance", "Mollusc abundance", "Annelid abundance"), ordered = T)

# Create plot
p4 <- ggplot(Group4, aes(x = Estimate, y = fDriver)) +
  scale_color_identity() +
  geom_point(
    shape = 21, size = 4, fill = "white",
    aes(color = ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
  ) +
  geom_errorbar(
    width = 0, linewidth = 1,
    aes(xmin = (`2.5 %`), xmax = (`97.5 %`), color = ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
  ) +
  facet_wrap(~ fResponse, nrow = 5, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"),
    legend.position = "none",
    text = element_text(size = 12),  # Adjust the size of all text elements
    axis.title = element_blank(),    # Remove axis titles
    axis.text = element_text(size = 12),  # Adjust the size of tick mark labels
    strip.text = element_text(size = 12, face = "bold")  # Adjust the size and face of facet labels
  ) +
  geom_vline(xintercept = 0, lty = 3) +
  ylab("") +
  xlab("")
p4

# Change response names to match other plots (consistency)
nrow(Group5)
fResponse <- data.frame(fResponse = factor(c(rep("Standardised func. dispersion", 6), rep("Standardised func. evenness", 6), rep("Standardised func. richness", 6), rep("Rarified taxon richness", 6))))
Group5 <- cbind(Group5, fResponse)
Group5$fResponse <- factor(Group5$fResponse, levels=c("Rarified taxon richness", "Standardised func. richness", "Standardised func. evenness", "Standardised func. dispersion"), ordered = T) # reorder the driver names so that they appear in ggplot the way we want

# Create plot
p5 <- ggplot(Group5, aes(x = Estimate, y = fDriver)) +
  scale_color_identity() +
  geom_point(
    shape = 21, size = 4, fill = "white",
    aes(color = ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
  ) +
  geom_errorbar(
    width = 0, linewidth = 1,
    aes(xmin = (`2.5 %`), xmax = (`97.5 %`), color = ifelse(Estimate >= 0, "#95ccba", "#f2cc84"))
  ) +
  facet_wrap(~ fResponse, nrow = 5, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"),
    legend.position = "none",
    text = element_text(size = 12),  # Adjust the size of all text elements
    axis.title = element_blank(),    # Remove axis titles
    axis.text = element_text(size = 12),  # Adjust the size of tick mark labels
    strip.text = element_text(size = 12, face = "bold")  # Adjust the size and face of facet labels
  ) +
  geom_vline(xintercept = 0, lty = 3) +
  ylab("") +
  xlab("")
p5

# create common x and y labels
x_axis_label <- "Estimate"
# Combine ggplot 1
(combined_1 <- cowplot::plot_grid(p1, p2, align = "hv", axis = "bt", ncol = 2))
# Add common label
combined_1 <- ggdraw() +
  draw_plot(combined_1) +
  draw_plot_label(label = x_axis_label, y = 0.03, x = 0.45, fontface = "plain", size = 12)
print(combined_1)
# # save plots
# tiff(filename = "Plots/LT_Overall_Driver_Ests_TaxoFuncIndices.tiff", width = 8, height = 12, units = 'in', res = 600, compression = 'lzw')
# combined_1
# dev.off()

# Combine ggplot 2
(combined_2 <- cowplot::plot_grid(p3, p4, align = "hv", axis = "bt", ncol = 2))
# Add common label
combined_2 <- ggdraw() +
  draw_plot(combined_2) +
  draw_plot_label(label = x_axis_label, y = 0.03, x = 0.45, fontface = "plain", size = 12)
print(combined_2)
# save plots
# tiff(filename = "Plots/LT_Overall_Driver_Ests_TaxoGroups.tiff", width = 8, height = 12, units = 'in', res = 600, compression = 'lzw')
# combined_2
# dev.off()

# extra plot
# tiff(filename = "Plots/LT_Overall_Driver_Ests_Extra.tiff", width = 4, height = 9, units = 'in', res = 600, compression = 'lzw')
# p5
# dev.off()

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
