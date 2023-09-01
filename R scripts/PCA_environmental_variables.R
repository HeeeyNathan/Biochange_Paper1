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
allYrs_pca = princomp(na.omit(allYrs_renamed[, c(79:91)]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", cex.axis = 1.5,
     main = "Principle coordinate", xlab = "PC1", ylab = "PC 2", cex.main = 2, cex.lab = 1.5) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
points(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], col = "#95ccba", pch = 19, cex = 1.5)  # You can customize col, pch, and cex
vec = envfit(allYrs_pca, na.omit(allYrs_renamed[, c(79:91)]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "#ECB246", cex = 1.5) #plots vectors (environmetal variables)
# Adding the box with text
text_box <- c("Cumulative explained variance = 64.2%\nPC 1 variance = 51.4%\nPC 2 variance = 12.8%")
rect_x <- min(allYrs_pca$scores[, 1]) + 0.025
rect_y <- max(allYrs_pca$scores[, 2]) - 6
rect_width <- 2.55
rect_height <- 1
rect(rect_x, rect_y, rect_x + rect_width, rect_y - rect_height, border = "black", lwd = 1)
text(rect_x + 1.28, rect_y - 0.25, "Cumulative explained variance = 64.2%", cex = 1, col = "black", font = 2)
text(rect_x + 1.175, rect_y - 0.5, "Variance explained by PC1 = 51.4%", cex = 1, col = "black", font = 2)
text(rect_x + 1.175, rect_y - 0.75, "Variance explained by PC2 = 12.8%", cex = 1, col = "black", font = 2)
