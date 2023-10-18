##### load functions, packages, and external sources #####
library(pacman)
library(vegan)
library(RColorBrewer)
library(dplyr)

##### read in data & format #####
d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T)
allYrs <- d1[!is.na(d1$site_id_wMissing),]

# Remove the year 2019 from all non-within-site analyses
allYrs <- allYrs[allYrs$year != "2019",]

# tranform predictor variables (if necessary)
hist(allYrs$flow)
hist(allYrs$temp)

allYrs$flow <- log10(allYrs$flow)

# Change river type variables
allYrs$river_type[allYrs$river_type == 1] <- "type1"
allYrs$river_type[allYrs$river_type == 2] <- "type2"
allYrs$river_type[allYrs$river_type == 3] <- "type3"
allYrs$river_type[allYrs$river_type == 4] <- "type4"
allYrs$river_type[allYrs$river_type == 5] <- "type5"

# make factors
allYrs$fYear <- factor(allYrs$year_wMissing)
allYrs$fmodified <- as.factor(allYrs$Heavily_modified)
allYrs$ftype <- as.factor(allYrs$river_type)
allYrs$fEQC <- as.factor(allYrs$EQC)

# Scale the covariates
# function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- lapply(df, function(x) if(is.numeric(x)){
    scale(x, center=TRUE, scale=TRUE)
  } else x)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df[,c(1, 2)], newd)
}

#apply function
sPredictors <- scaleVars(allYrs[, c(3:4, 58:72)])
sPredictors <- subset(sPredictors, select = -c(site_id, year, ssite_id, syear)) # remove ID variable
sPredictors$ID <- rownames(sPredictors)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, sPredictors, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# Rename columns using dplyr
allYrs_renamed <- allYrs %>%
  rename("pH_unchanged" = pH,
         "Suspended solids" = ssus_solid,
         "Alkalinity" = salkalinity,
         "Dissolved oxygen" = so2_dis,
         "pH" = spH,
         "Electrical conductivity" = sEC,
         "Biological oxygen demand" = sBOD7,
         "Ammonium" = sNH4.N,
         "Nitrite" =  sNO2.N,
         "Nitrate" = sNO3.N,
         "Mineralized nitrogen" = smineral.N,
         "Total nitrogen" = sTot.N,
         "Phosphate" = sPO4.P,
         "Total phosphorus" = sTot.P)

# Construct PCA to check the environmental variables and their relationships
allYrs_pca = princomp(na.omit(allYrs_renamed[, c("Alkalinity", "Electrical conductivity", 
                                                 "Nitrate", "Nitrite", "Mineralized nitrogen", "Total nitrogen",
                                                 "Phosphate", "Total phosphorus")]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
svg("Plots/Environmental_pca_new.svg", height = 15, width = 15)
# tiff("Plots/Environmental_pca.tiff", height = 15, width = 15, unit = "in", res = 600, compression = "lzw")
# Define a pastel color palette
color_palette <- c("#84bcda", "#ecc30b", "#f37748", "#283845", "#b6244f")
# Define the two point shapes you want to use (adjust as needed)
shape_palette <- c(19, 17)
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", cex.axis = 1.5,
     main = "Principal coordinate analysis of collinear environmental data", xlab = "PC1", ylab = "PC 2", cex.main = 2, cex.lab = 1.5) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
points(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], col = color_palette[allYrs_renamed$ftype], pch = shape_palette[allYrs_renamed$fmodified], cex = 1.5)  # You can customize col, pch, and cex
vec = envfit(allYrs_pca, na.omit(allYrs_renamed[, c("Alkalinity", "Electrical conductivity", 
                                                    "Nitrate", "Nitrite", "Mineralized nitrogen", "Total nitrogen",
                                                    "Phosphate", "Total phosphorus")]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "#ECB246", cex = 1.5, font = 2, pos = 4) #plots vectors (environmetal variables)
# Adding the box with text
rect_x <- min(allYrs_pca$scores[, 1]) + 9.6
rect_y <- max(allYrs_pca$scores[, 2]) + 0.2
rect_width <- 5.6
rect_height <- 0.8
rect(rect_x, rect_y, rect_x + rect_width, rect_y - rect_height, border = "black", lwd = 1)
text(rect_x + 2.8, rect_y - 0.2, "Cumulative explained variance = 90.3%", cex = 1.5, col = "black", font = 2)
text(rect_x + 2.55, rect_y - 0.4, "Variance explained by PC1 = 77.7%", cex = 1.5, col = "black", font = 2)
text(rect_x + 2.55, rect_y - 0.6, "Variance explained by PC2 = 12.6%", cex = 1.5, col = "black", font = 2)
legend("bottomleft", legend = levels(allYrs_renamed$ftype), fill = color_palette, title = "River type", cex = 1.5)
legend("topleft", legend = levels(allYrs_renamed$fmodified), pch = shape_palette, title = "Heavily modified", cex = 1.5)
dev.off()

# Check broken stick model
biplot(allYrs_pca)
allYrs_pca_scores <- (scores(allYrs_pca)[,1])
corr1 <- cor(na.omit(allYrs[, c(80, 83, 86:91)]), scores(allYrs_pca)[,1]) # correlation between original variables and principal components
round(corr1, 3)
corr2 <- cor(na.omit(allYrs[, c(80, 83, 86:91)]), scores(allYrs_pca)[,2]) # correlation between original variables and principal components
round(corr2, 3)
screeplot(allYrs_pca, bstick = TRUE, npcs = length(allYrs_pca$sdev))
(ev <- allYrs_pca$sdev^2)
n <- length (ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))}
bsm$p <- 100*bsm$p/n
bsm
svg("Plots/BrokStick_environmental_pca.svg", height = 5, width = 5)
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, main="Broken stick model", col=c("#95ccba",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch=15, col=c("#95ccba",2), bty="n")
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