# attach data
response_gls <- readRDS("Outputs/glsTrends_site_level_drivers.rds")
head(response_gls)
unique(response_gls$Response)
length(unique(response_gls$site_id))

Ests <- read.csv("Outputs/LT_Yr_metaanaly_weighted_noRandom_Ests_drivers.csv")
Ests <- Ests[, -1]

d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

library(vegan)

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

# # Change flow to discharge
# names(allYrs)[names(allYrs) == "flow"] <- "discharge"

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

# Construct PCA to check the environmental variables and their relationships
allYrs_pca = princomp(na.omit(allYrs[, c("salkalinity", "sEC", 
                                         "sNO3.N", "sNO2.N", "smineral.N", "sTot.N",
                                         "sPO4.P", "sTot.P")]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
# Define a pastel color palette
color_palette <- c("#84bcda", "#ecc30b", "#f37748", "#283845", "#b6244f")
# Define the two point shapes you want to use (adjust as needed)
shape_palette <- c(19, 17)
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", cex.axis = 1.5,
     main = "Principal coordinate analysis of collinear environmental data", xlab = "PC1", ylab = "PC 2", cex.main = 2, cex.lab = 1.5) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
points(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], col = color_palette[allYrs$ftype], pch = shape_palette[allYrs$fmodified], cex = 1.5)  # You can customize col, pch, and cex
vec = envfit(allYrs_pca, na.omit(allYrs[, c("salkalinity", "sEC", 
                                            "sNO3.N", "sNO2.N", "smineral.N", "sTot.N",
                                            "sPO4.P", "sTot.P")]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "#ECB246", cex = 1.5, font = 2, pos = 4) #plots vectors (environmetal variables)

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
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, main="Broken stick model", col=c("#95ccba",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch=15, col=c("#95ccba",2), bty="n")

# isolate pc axis for use as covariate in model
pc1_scores <- as.data.frame(allYrs_pca$scores[, 1]) # reverse signs to help with interpretation later
colnames(pc1_scores) <- c("PC_axis1")
pc1_scores$ID <- rownames(pc1_scores)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, pc1_scores, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# Create plot
tiff(filename = "Plots/LT_slopeDistributions_Drivers_updated.tiff", width = 12, height = 6, units = 'in', res = 600, compression = 'lzw')
# svg(filename = "Plots/LT_slopeDistributions_TaxoFuncIndices.svg", width = 12, height = 10, bg = "white")

# set plotting parametres
par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,2))

#### Flow #####
# log10
flow <- subset(response_gls, Response == "flow")
flow <- flow$estimate[!is.na(flow$estimate)]
mean(allYrs$flow, na.rm = T)
ave_flow <- mean(allYrs$flow, na.rm = T)
percChange_perYr <- (10^flow-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
flow_Est <- subset(Ests, Response == "flow")
stand_flow <- lapply(flow_Est[,2:11], function(x) (10^x-1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (3/4*(b-a)+a)
points(x=stand_flow$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_flow$Q2.5, stand_flow$Q2.5, stand_flow$Q97.5, stand_flow$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_flow$Q5, stand_flow$Q5, stand_flow$Q95, stand_flow$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_flow$Q10, stand_flow$Q10, stand_flow$Q90, stand_flow$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-18, y=(4/4*(b-a)+a), legend=expression(paste("a, Discharge (",m^3, sep = ".", s^-1,")")), bty="n", cex=1.3)
legend(x=5, y=(4/4*(b-a)+a), legend=expression(paste("-0.44%",y^-1, sep = "")), bty="n", cex=1.3)

#### Temperature #####
# none
temp <- subset(response_gls, Response == "temp")
temp <- temp$estimate[!is.na(temp$estimate)]
mean(allYrs$temp, na.rm = T)
ave_temp <- mean(allYrs$temp, na.rm = T)
percChange_perYr <- (temp/ave_temp)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
temp_Est <- subset(Ests, Response == "temp")
stand_temp <- lapply(temp_Est[,2:11],"*",100/ave_temp)
yy <- (1.7/4*(b-a)+a)
points(x=stand_temp$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_temp$Q2.5, stand_temp$Q2.5, stand_temp$Q97.5, stand_temp$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_temp$Q5, stand_temp$Q5, stand_temp$Q95, stand_temp$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_temp$Q10, stand_temp$Q10, stand_temp$Q90, stand_temp$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(2.6/4*(b-a)+a), legend=("b, Temperature (°C)"), bty="n", cex=1.3)
legend(x=5, y=(2.6/4*(b-a)+a), legend=expression(paste("+0.45% ",y^-1, sep = "")), bty="n", cex=1.3)

#### pH #####
# none
pH <- subset(response_gls, Response == "pH")
pH <- pH$estimate[!is.na(pH$estimate)]
mean(allYrs$pH, na.rm = T)
ave_pH <- mean(allYrs$pH, na.rm = T)
percChange_perYr <- (pH/ave_pH)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*0
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
pH_Est <- subset(Ests, Response == "pH")
stand_pH <- lapply(pH_Est[,2:11],"*",100/ave_pH)
yy <- (0.4/4*(b-a)+a)
points(x=stand_pH$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_pH$Q2.5, stand_pH$Q2.5, stand_pH$Q97.5, stand_pH$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_pH$Q5, stand_pH$Q5, stand_pH$Q95, stand_pH$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_pH$Q10, stand_pH$Q10, stand_pH$Q90, stand_pH$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(1.3/4*(b-a)+a), legend=("c, pH"), bty="n", cex=1.3)
legend(x=5, y=(1.3/4*(b-a)+a), legend=expression(paste("+0.18% ",y^-1, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

#### o2_dis #####
# log10
o2_dis <- subset(response_gls, Response == "o2_dis")
o2_dis <- o2_dis$estimate[!is.na(o2_dis$estimate)]
mean(allYrs$o2_dis, na.rm = T)
ave_o2_dis <- mean(allYrs$o2_dis, na.rm = T)
percChange_perYr <- (o2_dis/ave_o2_dis)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,15),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
o2_dis_Est <- subset(Ests, Response == "o2_dis")
stand_o2_dis <- lapply(o2_dis_Est[,2:11],"*",100/ave_o2_dis) # this must be the same as how you calculated the percChange_perYr
yy <- (3/4*(b-a)+a)
points(x=stand_o2_dis$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_o2_dis$Q2.5, stand_o2_dis$Q2.5, stand_o2_dis$Q97.5, stand_o2_dis$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_o2_dis$Q5, stand_o2_dis$Q5, stand_o2_dis$Q95, stand_o2_dis$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_o2_dis$Q10, stand_o2_dis$Q10, stand_o2_dis$Q90, stand_o2_dis$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-26, y=(4/4*(b-a)+a), legend=expression(paste("d, Dissolved oxygen (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=2, y=(4/4*(b-a)+a), legend=expression(paste("+1.24%",y^-1, sep = "")), bty="n", cex=1.3)

#### NH4.N #####
# none
NH4.N <- subset(response_gls, Response == "NH4.N")
NH4.N <- NH4.N$estimate[!is.na(NH4.N$estimate)]
mean(allYrs$NH4.N, na.rm = T)
ave_NH4.N <- mean(allYrs$NH4.N, na.rm = T)
percChange_perYr <- (10^NH4.N-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
NH4.N_Est <- subset(Ests, Response == "NH4.N")
stand_NH4.N <- lapply(NH4.N_Est[,2:11], function(x) (10^x-1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (1.7/4*(b-a)+a)
points(x=stand_NH4.N$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_NH4.N$Q2.5, stand_NH4.N$Q2.5, stand_NH4.N$Q97.5, stand_NH4.N$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_NH4.N$Q5, stand_NH4.N$Q5, stand_NH4.N$Q95, stand_NH4.N$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_NH4.N$Q10, stand_NH4.N$Q10, stand_NH4.N$Q90, stand_NH4.N$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-26, y=(2.6/4*(b-a)+a), legend=expression(paste("e, Ammonium (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=2, y=(2.6/4*(b-a)+a), legend=expression(paste("-2.86% ",y^-1, sep = "")), bty="n", cex=1.3)

#### PC_axis1 #####
# log10 + 3.5
PC_axis1 <- subset(response_gls, Response == "PC_axis1")
PC_axis1 <- PC_axis1$estimate[!is.na(PC_axis1$estimate)]
mean(allYrs$PC_axis1, na.rm = T)
ave_PC_axis1 <- mean(allYrs$PC_axis1, na.rm = T)
percChange_perYr <- (10^PC_axis1-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*0
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
PC_axis1_Est <- subset(Ests, Response == "PC_axis1")
stand_PC_axis1 <- lapply(PC_axis1_Est[,2:11], function(x) (10^x-1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.4/4*(b-a)+a)
points(x=stand_PC_axis1$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_PC_axis1$Q2.5, stand_PC_axis1$Q2.5, stand_PC_axis1$Q97.5, stand_PC_axis1$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_PC_axis1$Q5, stand_PC_axis1$Q5, stand_PC_axis1$Q95, stand_PC_axis1$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_PC_axis1$Q10, stand_PC_axis1$Q10, stand_PC_axis1$Q90, stand_PC_axis1$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-26, y=(1.3/4*(b-a)+a), legend=("f, Nutrients PCA axis 1"), bty="n", cex=1.3)
legend(x=2, y=(1.3/4*(b-a)+a), legend=expression(paste("1.50% ",y^-1, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

dev.off()
########################################

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
