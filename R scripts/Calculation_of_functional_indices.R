##### LOAD LIBRARIES --------------------
# Load the required packages
library(ade4)
library(adespatial) # beta diversity
library(FD) # functional diversity 
library(geodist) # geographic distance matrix
library(ggplot2) # plotting
library(ks) # kernal density estimation
library(pacman) # package manager
library(plyr)
library(RColorBrewer) #colours
library(vegan) # community ecology, taxonomic diversity
library(data.table) # for df manipulation
library(BAT) 

##### LOAD ADDITIONAL FUNCTIONS --------------------
# Additional functions
source("0_FD_functions.R")
source("0_quality_funct_space_fromdist.R")

##### LOAD THE DATA --------------------
#### TAXONOMIC DATA
### Full site-by-species matrix for taxonomic analysis
comm_raw <- read.csv("Data/LT_taxalist_2010-2020_short_transposed.csv", header = TRUE) # Add "LT_taxalist_2010-2020_short_final_transposed.csv"
comm_raw <- comm_raw[, -21] # remove 'Ametropus_sp.' because it has no trait info
comm <- comm_raw[, -1] # remove some columns that are needed now
rownames(comm) <- comm_raw$sample_id

### Grouped by year
Yr.mat <- read.table("Data/LT_annual_taxalist.csv", header = TRUE, sep = ",", row.names = 1)
Yr.mat <- Yr.mat[,-20] # remove 'Ametropus_sp.' because it has no trait information

#### FUNCITONAL DATA
### Full species-by-trait matrix for functional analysis
traits_raw <- read.csv("Data/LT_traits_2010-2020_short.csv", header = TRUE) # Add "LT_taxalist_2010-2020_tachet_bio_traits.csv"
traits_raw <- traits_raw[-20, ] # remove 'Ametropus_sp.' because it has no trait info
traits <- traits_raw[, -c(1:5)] # remove some columns that are needed now
rownames(traits) <- traits_raw$taxon_name # make row names the taxon names

## Changes trait data to proportions (between 0 & 1)
# Trait category blocks
traits <- prep.fuzzy.var(traits, c(4, 4, 8, 9, 2, 8, 7, 3, 8, 5, 5)) # These numbers show the number of traits in each trait group
rowSums(traits) #rowsums should sum to the number of trait groups you included

### functional site-by-species matrix
identical(row.names(traits), colnames(comm)) # check if names are the same

#### FUNCTIONAL ALPHA DIVERSITY --------------------
### CREATING THE SUPRAREGIONAL FS (PCoA) #####
SupReg_traits <- traits[(intersect(rownames(traits), colnames(comm))),]

# Gower dissimilarity
trait_df <- ktab.list.df(list(SupReg_traits[which(rowSums(SupReg_traits) == 11),]))
tr.dist <- dist.ktab(trait_df, type = c("F"))

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist, nbdim = 15) 
qual_fs$meanSD # 0.006 at 10D

# Supraregional FS (PCoA)
supreg.pco <- dudi.pco(tr.dist, scan = F, nf = 11)
summary(supreg.pco)

cumsum(supreg.pco$eig)[6]/sum(supreg.pco$eig)*100 # Variance explained by FS (PCoA); 54.31%

# Validating the PCoA
biplot(supreg.pco)
summary(supreg.pco)
screeplot(supreg.pco, bstick = TRUE, npcs = length(supreg.pco$eig))
(ev <- supreg.pco$eig^2)
n <- length (ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))}
bsm$p <- 100*bsm$p/n
bsm
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, main="Broken stick model", col=c("blue",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch=15, col=c("blue",2), bty="n")

# Spearman rank correlations between original trait categories and functional space axes
cor.res <- round(cor(SupReg_traits[which(rowSums(SupReg_traits) == 11),], supreg.pco$li, method = "spearman"), 2) #11 trait groups used
# write.csv(cor.res, "Outputs/FuncDiv_correlations.csv")

## FIGURE 3A: PROBABILITY DENSITY USING KERNEL DENSITY ESTIMATION #####
# Spearman rank correlations between original trait categories and functional space axes, which provides the contribution of each trait on the PCoA axes
cor.res
supreg.pco$li

# functional space
# Use dudi.pco for PCA for being consistent with scaling of scores in analysis
pco12 <- supreg.pco$li[,1:2]

# Kernel Density Estimation
Ho12 <- Hpi(x = pco12) # optimal bandwidth estimation
esto12 <- kde(x = pco12, H = Ho12, compute.cont = TRUE)     # kernel density estimation
plot(esto12)

# Set contour probabilities for drawing contour levels
clo12 <- contourLevels(esto12, prob = c(0.5, 0.05, 0.001), approx = TRUE)

# Setting arrows based on correlations (see cor.res -- 0.5 was the threshold)
# pc1:All trait categories selected on axis1 -> tachfeed_scr (-0.59), 
#                                               tachfood_lima (-0.51), 
#                                               tachfood_lmac (0.56), 
#                                               tachlifedur_less1 (-0.74), 
#                                               tachlifedur_gr1 (0.74), 
#                                               tachsize_0.5to1 (-0.61) 
# pc2:All trait categories selected on axis2->  tachdisp_aquapass (0.55), 
#                                               tachdisp_aeract (-0.61), 
#                                               tachloco_crw (-0.67), 
#                                               tachrepcyc_gr1 (0.51), 
#                                               tachrepro_ovo (0.55)
# pc3:All trait categories selected on axis3->  tachaqua_adult (0.7), 
#                                               tachfeed_pre (-0.51)
# pc4:All trait categories selected on axis4->  tachloco_flr (0.51), 
#                                               tachloco_wsw (0.55),
#                                               tachrespir_gil (-0.57),
#                                               tachrespir_spi (0.66)
# pc6:All trait categories selected on axis6->  tachresist_did (0.62), 
#                                               tachresist_non (-0.63)

cor.res # select traits with highest correlation (e.g. >= 0.5, <= -0.5)
SupReg_traits_sel12 <- subset(SupReg_traits, select = c(tachfeed_scr, 
                                                      tachfood_lima, 
                                                      tachfood_lmac, 
                                                      tachlifedur_less1, 
                                                      tachlifedur_gr1, 
                                                      tachsize_0.5to1,
                                                      
                                                      tachdisp_aquapass, 
                                                      tachdisp_aeract, 
                                                      tachloco_crw, 
                                                      tachrepcyc_gr1, 
                                                      tachrepro_ovo))

fit12 <- envfit(pco12, SupReg_traits_sel12)  # use envfit to draw arrows, can be also done using trait loadings

# Plot Kernel Density Estimations
plot(esto12, cont = seq(1, 100, by = 2.25), display  = "filled.contour", add = FALSE, ylab = "PC2", xlab = "PC1", 
     cex.axis = 0.75, ylim = c(-0.4, 0.6), xlim = c(-0.4, 0.6) , las = 1)

## FIGURE 1: RAIT PROBABILITY DENSITY PLOT #####
# pdf(file="Plots/Supraregional_FS_density_vectors.pdf",onefile = T, width = 6, height = 6)
# svg(file="Plots/Supraregional_FS_density_vectors.svg", onefile = T, width = 12, height = 12)
# tiff(file="Plots/Supraregional_FS_density_vectors.tiff", width = 12, height = 12, units = 'in', res = 600, compression = 'lzw')
par(mfrow = c(1,1), cex.axis = 1.85, cex.lab = 2, cex.main = 2, mar = c(5,5,4,1))
plot(esto12, cont = seq(1, 100, by = 2.25), display = "filled.contour", add = FALSE, 
     ylab = paste("PCo-2"), 
     xlab = paste("PCo-1"),
     main = paste("Supraregional FS"),
     #    ylab = paste("PC2","(",round(supreg.pco$eig[2]/sum(supreg.pco$eig), 3),"%)"), 
     #    xlab = paste("PC1","(",round(supreg.pco$eig[1]/sum(supreg.pco$eig), 3),"%)"), 
     cex.axis = 0.75, ylim=c(-1, 1), xlim = c(-1, 1), las = 1)
abline(h = 0, lty = 2, col = grey(0.5, alpha=0.2))
abline(v = 0, lty = 2, col = grey(0.5, alpha=0.2))
plot(esto12, abs.cont = clo12[1], labels = c(0.5), labcex = 0.75, add = TRUE, lwd = 0.75, col = "grey30")
plot(esto12, abs.cont = clo12[2], labels = c(0.95), labcex = 0.75, add = TRUE, lwd = 0.5, col = "grey60")
plot(esto12, abs.cont = clo12[3], labels = c(0.99), labcex = 0.75, add = TRUE, lwd = 0.5, col = "grey60")
points(pco12[,], pch = 16, cex = 0.25, col = "black")
plot(fit12, cex = 0.90, col = 1)
# dev.off()

## FIGURE 2: POSITION OF TAXONOMIC GROUPS WITHIN THE FS #####
library(adegraphics)
supra.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(comm))),]

# pdf(file = "Plots/Supraregional_FS_taxa_grouping.pdf", onefile = T, width = 6, height = 6)
svg(file = "Plots/Supraregional_FS_taxa_grouping.svg", onefile = T, width = 7, height = 7)
# tiff(file="Plots/Supraregional_FS_taxa_grouping.tiff", width = 12, height = 12, units = 'in', res = 600, compression = 'lzw')
par(mfrow = c(1,1), mar = c(5,5,4,1))
# Define a vector of pastel colors
pastel_colors <- c("#FFC0CB", "#FFD700", "#87CEEB", "#98FB98", "#FFA07A", "#9370DB", "#F0E68C", "#FF69B4", "#00CED1", "#B0E0E6", "#DDA0DD", "#20B2AA", "#FFE4B5", "#00FF7F", "#7B68EE", "#AFEEEE", "#F08080", "#40E0D0")
s.class(supreg.pco$li, fac = as.factor(supra.gr$group2[which(rowSums(traits)==11)]), plines.col = 1:18, col = pastel_colors)
# s.class(supreg.pco$li, fac = as.factor(supra.gr$group2[which(rowSums(traits)==11)]), plines.col = 1:18, col = T)
dev.off()

## FIGURE 3: REPRESENTATION OF THE SUPRAREGIONAL FS #####
# pdf(file = "Plots/Supraregional_FS_convexhull.pdf",onefile = T, width = 12, height = 16) 
# svg(file = "Plots/Supraregional_FS_convexhull.svg",onefile = T, width = 12, height = 12) 
# tiff(file="Plots/Supraregional_FS_convexhull.tiff", width = 12, height = 12, units = 'in', res = 600, compression = 'lzw')
par(mfrow = c(1,1), cex.axis = 1.85, cex.lab = 2, cex.main = 2, mar = c(5,5,3,1))
plot(range(supreg.pco$li[1]), range(supreg.pco$li[2]),type = "n", main="Supraregional FS",xlab = "PCo-1", cex.axis = 1, ylab = "PCo-1")
points(supreg.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
plot_chull2D(supreg.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
colMeans(supreg.pco$li[,c(1,2)])->cent_ov
points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
colMeans(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(comm))),c(1,2)])->cent_r
points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)
# dev.off()

# combined plots
svg(file = "Plots/funcitonal_space_largerlabs.svg",onefile = T, width = 21, height = 7) 
par(mfrow = c(1,3), cex.axis = 1.85, cex.lab = 2, cex.main = 2, mar = c(5,5,3,1))
# Supraregional functional space
plot(range(supreg.pco$li[1]), range(supreg.pco$li[2]), type = "n", cex.axis = 1.5, cex.lab = 1.5,
     xlab = paste("PCo-1","(",round(supreg.pco$eig[1]/sum(supreg.pco$eig),3)*100,"%)"),
     ylab = paste("PCo-2","(",round(supreg.pco$eig[2]/sum(supreg.pco$eig),3)*100,"%)"))
points(supreg.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
plot_chull2D(supreg.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
colMeans(supreg.pco$li[,c(1,2)])->cent_ov
points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
colMeans(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(comm))),c(1,2)])->cent_r
points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)
#Probability denisty
plot(esto12, cont = seq(1, 100, by = 2.25), display = "filled.contour", add = FALSE, 
     xlab = paste("PCo-1","(",round(supreg.pco$eig[1]/sum(supreg.pco$eig), 3)*100,"%)"),
     ylab = "",
     cex.lab = 1.5,
     cex.axis = 1.5, ylim=c(-1, 1), xlim = c(-1, 1), las = 1)
abline(h = 0, lty = 2, col = grey(0.5, alpha=0.2))
abline(v = 0, lty = 2, col = grey(0.5, alpha=0.2))
plot(esto12, abs.cont = clo12[1], labels = c(0.5), labcex = 0.75, add = TRUE, lwd = 0.75, col = "grey30")
plot(esto12, abs.cont = clo12[2], labels = c(0.95), labcex = 0.75, add = TRUE, lwd = 0.5, col = "grey60")
plot(esto12, abs.cont = clo12[3], labels = c(0.99), labcex = 0.75, add = TRUE, lwd = 0.25, col = "grey70")
points(pco12[,], pch = 16, cex = 0.25, col = "grey30")
plot(fit12, cex = 1.25, col = 1)
dev.off()

## FIGURE 4: YEARLY TAXON POOL REPRESENTATION WITHIN THE SUPRAREGIONAL FSs #####
# pdf(file = "Plots/Change_in_FS_through_time.pdf",onefile = T, width = 12, height = 16) 
# svg(file = "Plots/Change_in_FS_through_time.svg",onefile = T, width = 12, height = 16)
# tiff(file="Plots/Change_in_FS_through_time.tiff", width = 12, height = 16, units = 'in', res = 600, compression = 'lzw')
par(mfrow = c(4,3), cex.axis = 1.85, cex.lab = 2, cex.main = 2, mar = c(5,5,3,1))

c("royalblue4", "royalblue4", "royalblue4",
  "royalblue4", "royalblue4", "royalblue4", 
  "royalblue4", "royalblue4", "royalblue4",
  "royalblue4") -> col.reg
c("2010", "2011", "2012", 
  "2013", "2014", "2015", 
  "2016", "2017", "2019", 
  "2020") -> lab.reg

# probability denisty
plot(esto12, cont = seq(1, 100, by = 1.79), display = "filled.contour", add = FALSE, 
     ylab = paste("PCo-2"), 
     xlab = paste("PCo-1"),
     main = paste("Supraregional FS"),
     #    ylab = paste("PC2","(",round(supreg.pco$eig[2]/sum(supreg.pco$eig), 3),"%)"), 
     #    xlab = paste("PC1","(",round(supreg.pco$eig[1]/sum(supreg.pco$eig), 3),"%)"), 
     cex.axis = 0.75, ylim=c(-1, 1), xlim = c(-1, 1), las = 1)
abline(h = 0, lty = 2, col = grey(0.5, alpha=0.2))
abline(v = 0, lty = 2, col = grey(0.5, alpha=0.2))
plot(esto12, abs.cont = clo12[1], labels = c(0.5), labcex = 0.75, add = TRUE, lwd = 0.75, col = "grey30")
plot(esto12, abs.cont = clo12[2], labels = c(0.95), labcex = 0.75, add = TRUE, lwd = 0.5, col = "grey60")
plot(esto12, abs.cont = clo12[3], labels = c(0.99), labcex = 0.75, add = TRUE, lwd = 0.5, col = "grey60")
points(pco12[,], pch = 16, cex = 0.25, col = "black")
plot(fit12, cex = 0.90, col = 1)

# Supraregional functional space
plot(range(supreg.pco$li[1]), range(supreg.pco$li[2]),type = "n", main="Supraregional FS",xlab = "PCo-1", cex.axis = 1, ylab = "PCo-1")
points(supreg.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
plot_chull2D(supreg.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
colMeans(supreg.pco$li[,c(1,2)])->cent_ov
points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
colMeans(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(comm))),c(1,2)])->cent_r
points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)

# Yearly funcitonal spaces
for (i in 1:10)
  
{
  plot(range(supreg.pco$li[1]), range(supreg.pco$li[2]),type = "n", main=lab.reg[i], cex.axis = 1, xlab= "PCo-1", ylab="PCo-1")
  points(supreg.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
  plot_chull2D(supreg.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
  plot_chull2D(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(Yr.mat[which(Yr.mat[i,]>0)]))),c(1,2)],col_ch=adjustcolor(col.reg[i], alpha.f = 0.35),border_ch=adjustcolor(col.reg[i], alpha.f = 0.4))
  colMeans(supreg.pco$li[,c(1,2)])-> cent_ov
  points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
  colMeans(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(Yr.mat[which(Yr.mat[i,]>0)]))),c(1,2)])-> cent_r
  points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)
}

# dev.off()

# Building annual FSs (PCoA)
# Years
rownames(Yr.mat)

# Matrix to store correlations between continental and regional FS axes
cor.reg<-data.frame(matrix(NA, 10, 10))
rownames(cor.reg)<-c( "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2019", "2020")
colnames(cor.reg)<-paste("Axis", 1:10)

### Plotting temporal-scale FS
AllYrs.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat))),]
Yr2010.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[1,]>0)]))),]
Yr2011.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[2,]>0)]))),]
Yr2012.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[3,]>0)]))),]
Yr2013.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[4,]>0)]))),]
Yr2014.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[5,]>0)]))),]
Yr2015.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[6,]>0)]))),]
Yr2016.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[7,]>0)]))),]
Yr2017.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[8,]>0)]))),]
Yr2019.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[9,]>0)]))),]
Yr2020.gr <- data.frame(group2 = traits_raw$group2, traits)[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[10,]>0)]))),]

# Selecting taxa from the regional pool
Yr2010.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[1,]>0)]))),]
Yr2011.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[2,]>0)]))),]
Yr2012.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[3,]>0)]))),]
Yr2013.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[4,]>0)]))),]
Yr2014.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[5,]>0)]))),]
Yr2015.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[6,]>0)]))),]
Yr2016.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[7,]>0)]))),]
Yr2017.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[8,]>0)]))),]
Yr2019.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[9,]>0)]))),]
Yr2020.tr<-traits[(intersect(rownames(traits), colnames(Yr.mat[which(Yr.mat[10,]>0)]))),]

### 2010
# calculating annual Gower dissimilarity
traits_2010 <- ktab.list.df(list(Yr2010.tr[which(rowSums(Yr2010.tr)==11),]))
tr.dist_2010 <- dist.ktab(traits_2010, type= c("F"))

# FS (PCoA)
Yr2010.pco <- dudi.pco(tr.dist_2010, scan = F, nf = 10)
Yr2010.tax <- Yr.mat[1,(intersect(rownames(Yr2010.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2010, nbdim=15) 
qual_fs$meanSD

cumsum(Yr2010.pco$eig)[6]/sum(Yr2010.pco$eig)*100 # Variance explained by FS (PCoA); 55.81%

# Correlation betwewen regional and continental-scale FS axes
Yr2010.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2010.pco$li))),]
cor.reg[1,]<-round(diag(cor(Yr2010.pco$li , Yr2010.cfs, method="pearson")),3)

### 2011
# calculating annual Gower dissimilarity
traits_2011 <- ktab.list.df(list(Yr2011.tr[which(rowSums(Yr2011.tr)==11),]))
tr.dist_2011 <- dist.ktab(traits_2011, type= c("F"))

# FS (PCoA)
Yr2011.pco <- dudi.pco(tr.dist_2011, scan = F, nf = 10)
Yr2011.tax <- Yr.mat[1,(intersect(rownames(Yr2011.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2011, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2011.pco$eig)[6]/sum(Yr2011.pco$eig)*100 # Variance explained by FS (PCoA); 62.74%

# Correlation betwewen regional and continental-scale FS axes
Yr2011.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2011.pco$li))),]
cor.reg[2,]<-round(diag(cor(Yr2011.pco$li , Yr2011.cfs, method="pearson")),3)

### 2012
# calculating annual Gower dissimilarity
traits_2012 <- ktab.list.df(list(Yr2012.tr[which(rowSums(Yr2012.tr)==11),]))
tr.dist_2012 <- dist.ktab(traits_2012, type= c("F"))

# FS (PCoA)
Yr2012.pco <- dudi.pco(tr.dist_2012, scan = F, nf = 10)
Yr2012.tax <- Yr.mat[1,(intersect(rownames(Yr2012.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2012, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2012.pco$eig)[6]/sum(Yr2012.pco$eig)*100 # Variance explained by FS (PCoA); 53.86%

# Correlation betwewen regional and continental-scale FS axes
Yr2012.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2012.pco$li))),]
cor.reg[3,]<-round(diag(cor(Yr2012.pco$li , Yr2012.cfs, method="pearson")),3)

### 2013
# calculating annual Gower dissimilarity
traits_2013 <- ktab.list.df(list(Yr2013.tr[which(rowSums(Yr2013.tr)==11),]))
tr.dist_2013 <- dist.ktab(traits_2013, type= c("F"))

# FS (PCoA)
Yr2013.pco <- dudi.pco(tr.dist_2013, scan = F, nf = 10)
Yr2013.tax <- Yr.mat[1,(intersect(rownames(Yr2013.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2013, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2013.pco$eig)[6]/sum(Yr2013.pco$eig)*100 # Variance explained by FS (PCoA); 54.20%

# Correlation betwewen regional and continental-scale FS axes
Yr2013.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2013.pco$li))),]
cor.reg[4,]<-round(diag(cor(Yr2013.pco$li , Yr2013.cfs, method="pearson")),3)

### 2014
# calculating annual Gower dissimilarity
traits_2014 <- ktab.list.df(list(Yr2014.tr[which(rowSums(Yr2014.tr)==11),]))
tr.dist_2014 <- dist.ktab(traits_2014, type= c("F"))

# FS (PCoA)
Yr2014.pco <- dudi.pco(tr.dist_2014, scan = F, nf = 10)
Yr2014.tax <- Yr.mat[1,(intersect(rownames(Yr2014.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2014, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2014.pco$eig)[6]/sum(Yr2014.pco$eig)*100 # Variance explained by FS (PCoA); 54.99%

# Correlation betwewen regional and continental-scale FS axes
Yr2014.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2014.pco$li))),]
cor.reg[5,]<-round(diag(cor(Yr2014.pco$li , Yr2014.cfs, method="pearson")),3)

### 2015
# calculating annual Gower dissimilarity
traits_2015 <- ktab.list.df(list(Yr2015.tr[which(rowSums(Yr2015.tr)==11),]))
tr.dist_2015 <- dist.ktab(traits_2015, type= c("F"))

# FS (PCoA)
Yr2015.pco <- dudi.pco(tr.dist_2015, scan = F, nf = 10)
Yr2015.tax <- Yr.mat[1,(intersect(rownames(Yr2015.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2015, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2015.pco$eig)[6]/sum(Yr2015.pco$eig)*100 # Variance explained by FS (PCoA); 53.25%

# Correlation betwewen regional and continental-scale FS axes
Yr2015.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2015.pco$li))),]
cor.reg[6,]<-round(diag(cor(Yr2015.pco$li , Yr2015.cfs, method="pearson")),3)

### 2016
# calculating annual Gower dissimilarity
traits_2016 <- ktab.list.df(list(Yr2016.tr[which(rowSums(Yr2016.tr)==11),]))
tr.dist_2016 <- dist.ktab(traits_2016, type= c("F"))

# FS (PCoA)
Yr2016.pco <- dudi.pco(tr.dist_2016, scan = F, nf = 10)
Yr2016.tax <- Yr.mat[1,(intersect(rownames(Yr2016.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2016, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2016.pco$eig)[6]/sum(Yr2016.pco$eig)*100 # Variance explained by FS (PCoA); 54.76%

# Correlation betwewen regional and continental-scale FS axes
Yr2016.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2016.pco$li))),]
cor.reg[7,]<-round(diag(cor(Yr2016.pco$li , Yr2016.cfs, method="pearson")),3)

### 2017
# calculating annual Gower dissimilarity
traits_2017 <- ktab.list.df(list(Yr2017.tr[which(rowSums(Yr2017.tr)==11),]))
tr.dist_2017 <- dist.ktab(traits_2017, type= c("F"))

# FS (PCoA)
Yr2017.pco <- dudi.pco(tr.dist_2017, scan = F, nf = 10)
Yr2017.tax <- Yr.mat[1,(intersect(rownames(Yr2017.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2017, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2017.pco$eig)[6]/sum(Yr2017.pco$eig)*100 # Variance explained by FS (PCoA); 54.98%

# Correlation betwewen regional and continental-scale FS axes
Yr2017.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2017.pco$li))),]
cor.reg[8,]<-round(diag(cor(Yr2017.pco$li , Yr2017.cfs, method="pearson")),3)

### 2019
# calculating annual Gower dissimilarity
traits_2019 <- ktab.list.df(list(Yr2019.tr[which(rowSums(Yr2019.tr)==11),]))
tr.dist_2019 <- dist.ktab(traits_2019, type= c("F"))

# FS (PCoA)
Yr2019.pco <- dudi.pco(tr.dist_2019, scan = F, nf = 10)
Yr2019.tax <- Yr.mat[1,(intersect(rownames(Yr2019.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2019, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2019.pco$eig)[6]/sum(Yr2019.pco$eig)*100 # Variance explained by FS (PCoA); 71.36%

# Correlation betwewen regional and continental-scale FS axes
Yr2019.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2019.pco$li))),]
cor.reg[9,]<-round(diag(cor(Yr2019.pco$li , Yr2019.cfs, method="pearson")),3)

### 2020
# calculating annual Gower dissimilarity
traits_2020 <- ktab.list.df(list(Yr2020.tr[which(rowSums(Yr2020.tr)==11),]))
tr.dist_2020 <- dist.ktab(traits_2020, type= c("F"))

# FS (PCoA)
Yr2020.pco <- dudi.pco(tr.dist_2020, scan = F, nf = 10)
Yr2020.tax <- Yr.mat[1,(intersect(rownames(Yr2020.pco$li), colnames(Yr.mat)))]

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist_2020, nbdim=15) 
qual_fs$meanSD
cumsum(Yr2020.pco$eig)[6]/sum(Yr2020.pco$eig)*100 # Variance explained by FS (PCoA); 55.09%

# Correlation betwewen regional and continental-scale FS axes
Yr2020.cfs<- supreg.pco$li[(intersect(rownames(supreg.pco$li), rownames(Yr2020.pco$li))),]
cor.reg[10,]<-round(diag(cor(Yr2020.pco$li , Yr2020.cfs, method="pearson")),3)

### Plotting annual FS
# All Years
#pdf(file="continental_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(supreg.pco$li, fac = factor(AllYrs.gr$group2[which(rowSums(traits)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2010
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2010.pco$li, fac = factor(Yr2010.gr$group2[which(rowSums(Yr2010.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2011
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2011.pco$li, fac = factor(Yr2011.gr$group2[which(rowSums(Yr2011.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2012
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2012.pco$li, fac = factor(Yr2012.gr$group2[which(rowSums(Yr2012.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2013
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2013.pco$li, fac = factor(Yr2013.gr$group2[which(rowSums(Yr2013.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2014
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2014.pco$li, fac = factor(Yr2014.gr$group2[which(rowSums(Yr2014.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2015
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2015.pco$li, fac = factor(Yr2015.gr$group2[which(rowSums(Yr2015.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2016
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2016.pco$li, fac = factor(Yr2016.gr$group2[which(rowSums(Yr2016.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2017
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2017.pco$li, fac = factor(Yr2017.gr$group2[which(rowSums(Yr2017.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2019
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2019.pco$li, fac = factor(Yr2019.gr$group2[which(rowSums(Yr2019.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

# 2020
#pdf(file="rif_fs.pdf",onefile=T,width=4,height=4)  
par(mfrow=c(1,1),cex.axis=1.85,cex.lab=2,cex.main=2,mar=c(5,5,4,1))
s.class(Yr2020.pco$li, fac = factor(Yr2020.gr$group2[which(rowSums(Yr2020.tr)==11)]), plines.col = 1:18, col = TRUE)
#dev.off()

### COMPUTING FD METRICS (FRIC, FDIS, FEVE)#####
# Using supraregional-scale FS
chull_3d(supreg.pco$li, m = 6, prec = "Qt") -> ch_st # estimates the convex hull of a Functional Space
fric_3d(comm, supreg.pco$li, m = 6, prec = "QJ", fric.3d.max = ch_st) -> FRic_co
feve_k(supreg.pco$li, comm, m = 6) -> FEve_co
fdisp_k_sub(tr.dist, comm, tax_sub = colnames(comm), m = 6)$FDis -> FDis_co

# Substituting NAs by 0 (only if necessary)
FRic_co[which(is.na(FRic_co) == T)] <- 0
FEve_co[which(is.na(FEve_co) == T)] <- 0

### Add these indices to the a func_div data frame
func_div <- data.frame(FRic_co, FEve_co, FDis_co)
colnames(func_div) <- c("FRic", "FEve", "FDis")
func_div

## Computing FRed
comm
tr.dist
library(adiv)
# Calculated from the adiv package (Redundancy = 1-Uniqueness)
Redundancy <- uniqueness(comm, tr.dist, tol = 1e-08, abundance = TRUE) # creates 3 dataframs: kbar, V, and Red
# Isolate redundancy metrics
func_div$FRed <- Redundancy$red$R

### NULL MODELS OF FD METRICS #####
### shuffling null model: keeping the functional space the same, but shuffling taxa within the FS
# FRic
## inputs
comm
supreg.pco$li
ch_st # the maximum size of the convex hull

## model
FRic.shuff <- function(x){
  rownames(x) <- sample(rownames(x), length(rownames(x)), replace = F)
  x <- x[order(rownames(x)),]
  fric_3d(comm, x, m = 6, prec = "QJ", fric.3d.max = ch_st)
}
set.seed(1) # make results repeatable
FRic.obs.null.output.shuff <- cbind(fric_3d(comm, supreg.pco$li, m = 6, prec = "QJ", fric.3d.max = ch_st), 
                                    replicate(3, FRic.shuff(supreg.pco$li))) # change this number to 999 after you have determined the code runs

FRic.ses.value <- (FRic.obs.null.output.shuff[,1] - 
                     apply(FRic.obs.null.output.shuff, 1, mean))/
  apply(FRic.obs.null.output.shuff,1,sd)

qFRic <- NaN*FRic.obs.null.output.shuff[,1]

for(i in seq(qFRic)){
  qFRic[i] <- sum(FRic.obs.null.output.shuff[,1][i] > FRic.obs.null.output.shuff[i,]) / length(FRic.obs.null.output.shuff[i,])
}

sigFRic <- qFRic < 0.05 | qFRic > 0.95 # test if outside distribution

FRic_output <- as.data.frame(cbind(FRic.ses.value, sigFRic))
colnames(FRic_output) <- c("FRic.SES", "FRic.SES.sig")
boxplot(FRic.obs.null.output.shuff[,1], FRic_output[,1]) # obsevred vs standardised FRic

null_outputs <- FRic_output
# write.csv(FRic_output, "FRic_output.csv")

# FEve
## inputs
comm
supreg.pco$li

## model
FEve.shuff <- function(x){
  rownames(x) <- sample(rownames(x), length(rownames(x)), replace = F)
  x <- x[order(rownames(x)),]
  feve_k(x, comm, m = 6)
}
set.seed(1) # make results repeatable
FEve.obs.null.output.shuff <- cbind(feve_k(supreg.pco$li, comm, m = 6), 
                                    replicate(3, FEve.shuff(supreg.pco$li))) # change this number to 999 after you have determined the code runs

FEve.ses.value <- (FEve.obs.null.output.shuff[,1] - 
                     apply(FEve.obs.null.output.shuff, 1, mean))/
  apply(FEve.obs.null.output.shuff,1,sd)

qFEve <- NaN*FEve.obs.null.output.shuff[,1]

for(i in seq(qFEve)){
  qFEve[i] <- sum(FEve.obs.null.output.shuff[,1][i] > FEve.obs.null.output.shuff[i,]) / length(FEve.obs.null.output.shuff[i,])
}

sigFEve <- qFEve < 0.05 | qFEve > 0.95 # test if outside distribution

FEve_output <- as.data.frame(cbind(FEve.ses.value, sigFEve))
colnames(FEve_output) <- c("FEve.SES", "FEve.SES.sig")
boxplot(FEve.obs.null.output.shuff[,1], FEve_output[,1])#obsevred vs standardised FRic

null_outputs <- cbind(FRic_output, FEve_output)
# write.csv(FEve_output, "FEve_output.csv")

# FDisp
## inputs
comm
tr.dist

## model - changed species identities by shuffling the site-by-taxa matrix
FDis.shuff <- function(x){
  colnames(x) <- sample(colnames(x), length(colnames(x)), replace = F)
  x <- x[, order(names(x))]
  colnames(x) <- colnames(comm)
  fdisp_k(tr.dist, x, m = 6)$FDis
}
set.seed(1) # make results repeatable
FDis.obs.null.output.shuff <- cbind(fdisp_k(tr.dist, comm, m = 6)$FDis, 
                                    replicate(3, FDis.shuff(comm))) # change this number to 999 after you have determined the code runs

FDis.ses.value <- (FDis.obs.null.output.shuff[,1] - 
                     apply(FDis.obs.null.output.shuff, 1, mean))/
  apply(FDis.obs.null.output.shuff,1,sd)

qFDis <- NaN*FDis.obs.null.output.shuff[,1]

for(i in seq(qFDis)){
  qFDis[i] <- sum(FDis.obs.null.output.shuff[,1][i] > FDis.obs.null.output.shuff[i,]) / length(FDis.obs.null.output.shuff[i,])
}

sigFDis <- qFDis < 0.05 | qFDis > 0.95 # test if outside distribution

FDis_output <- as.data.frame(cbind(FDis.ses.value, sigFDis))
colnames(FDis_output) <- c("FDis.SES", "FDis.SES.sig")
boxplot(FDis.obs.null.output.shuff[,1], FDis_output[,1])# obsevred vs standardised FDis

# write.csv(FDis_output, "FDis_output.csv")

# Read in computed null model values
null.outputs <- read.csv("Outputs/null_model_outputs_wNAs.csv", h = TRUE, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)

# Substituting NAs by 0 (only if necessary)
null.outputs$FRic.SES[which(is.na(null.outputs$FRic.SES) == T)] <- 0
null.outputs$FEve.SES[which(is.na(null.outputs$FEve.SES) == T)] <- 0

# write.csv(null.outputs, "null_model_outputs.csv")

#### EXTRA FUNCTIONAL DIVERSITY INDICES & CWMS --------------------
### inputs
comm
traits
### check
identical(row.names(traits), colnames(comm)) # check if names are the same
CWM <- dbFD(traits, as.matrix(comm), calc.FRic = FALSE, calc.FDiv = FALSE, calc.CWM = TRUE, scale.RaoQ = TRUE, CWM.type = "all")
# write.csv(CWM, "Outputs/dbFD_output_wCWM.csv")

#### COMBINE ALL FUNCTIONAL INDICES
### Inputs
func_div
null.outputs
CWM

### Combine all data
func_div_all <- cbind(func_div, null.outputs, CWM$CWM)
## organise dataframe
# create site_code column
library(stringr)
sample_id_split <- as.data.frame(str_split_fixed(rownames(func_div_all), "_", 3))
site_code <- as.data.frame(paste(sample_id_split$V1, sample_id_split$V2, sep = "_"))
colnames(site_code)[1] <- "site_code"
func_div_all <- cbind(func_div_all, site_code)
# rearrange columns and change row names
func_div_all <- func_div_all[, c(74, 1:73)]
rownames(func_div_all) <- 1:335
func_div_all
write.csv(func_div_all, "Outputs/LT_siteYr_FuncDiversity.csv", row.names = FALSE) # save output

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

