library(tidyverse)
library(viridis)
library(wesanderson)
library(treemap)

my_data <- read.csv("Data/LT_taxalist_2010-2020_pie_chart_SF_test.csv")
my_data$logAbund <- log10(my_data$Abundance + 1)
colnames(my_data)

my_data_perc <- read.csv("Data/LT_taxalist_2010-2020_pie_chart_LF_perc.csv")
my_data_perc <- my_data_perc[, -1]

# Remove the year 2019 from all non-within-site analyses
my_data_perc <- my_data_perc[my_data_perc$Year != "2019",] # removes the intercept line from each response variable

# Treemap
## 2 groups
treemap(my_data_perc[1:24,], index=c("Phylum","SubPhylum_Class"), vSize ="Abundance", type="index",
        ontsize.labels=c(15,12),                          # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","orange"),    # Color of labels
        fontface.labels=c(2,1),                           # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        #bg.labels = c("transparent"),                     # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("right", "bottom")
        ),                                                # Where to place labels in the rectangle?
        overlap.labels=0.5,                               # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                                 # If true, labels are bigger when rectangle is bigger.
        border.col=c("black","white"),                    # Color of borders of groups, of subgroups, of subsubgroups ....
        border.lwds=c(3,2)                                # Width of colors
)

treemap(my_data_perc[1:24,], index=c("SubPhylum_Class", "SuperFam_Order"), vSize="Abundance", type="index",
        ontsize.labels=c(15,12),                          # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("orange","black"),    # Color of labels
        fontface.labels=c(2,1),                           # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        #bg.labels = c("transparent"),                     # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("right", "bottom")
        ),                                                # Where to place labels in the rectangle?
        overlap.labels=0.5,                               # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                                 # If true, labels are bigger when rectangle is bigger.
        border.col=c("black","white"),                    # Color of borders of groups, of subgroups, of subsubgroups ....
        border.lwds=c(3,2)                                # Width of colors
)


## 3 groups
treemap(my_data_perc[1:24,], index=c("Phylum", "SubPhylum_Class", "SuperFam_Order"), vSize="Abundance", type="index",
        ontsize.labels=c(15,12),                          # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","orange", "black"),    # Color of labels
        fontface.labels=c(2,1),                           # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        #bg.labels = c("transparent"),                     # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("right", "bottom")
        ),                                                # Where to place labels in the rectangle?
        overlap.labels=0.5,                               # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                                 # If true, labels are bigger when rectangle is bigger.
        border.col=c("black","white"),                    # Color of borders of groups, of subgroups, of subsubgroups ....
        border.lwds=c(3,2)                                # Width of colors
)

# Stacked barchart
# Phylum ~ Year
ggplot(my_data_perc, aes(fill=Phylum, y=Abundance, x=Year)) + 
  geom_bar(position="fill", stat="identity")

# SubPhylum_Class ~ Year
ggplot(my_data_perc, aes(fill=SubPhylum_Class, y=Abundance, x=Year)) + 
  geom_bar(position="fill", stat="identity")

# SuperFam_Order ~ Year
ggplot(my_data_perc, aes(fill=SuperFam_Order, y=Abundance, x=Year)) + 
  geom_bar(position="fill", stat="identity")

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
