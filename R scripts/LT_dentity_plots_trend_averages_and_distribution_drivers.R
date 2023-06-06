# attach data
response_gls <- readRDS("Outputs/glsTrends_site_level_drivers.rds")
head(response_gls)
unique(response_gls$Response)
length(unique(response_gls$site_id))

Ests <- read.csv("Outputs/LT_Yr_metaanaly_weighted_noRandom_Ests_drivers.csv")
Ests <- Ests[, -1]

d1 <- read.csv("Data/LT_siteYr_DivIndices_EnvVariables_final.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

library(vegan)

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
sPredictors <- scaleVars(allYrs[, c(3, 4, 67:81)])
sPredictors <- subset(sPredictors, select = -c(site_id, year, ssite_id, syear)) # remove ID variable
sPredictors$ID <- rownames(sPredictors)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, sPredictors, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

# define new variable for nutrients
# Construct PCA to check the environmental variables and their relationships
allYrs_pca = princomp(na.omit(allYrs[, c(94:106)]), scores = TRUE) #computes PCA - gives us all the PCA results
summary(allYrs_pca)
plot(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], type = "n", xlab = "Axis.1", ylab = "Axis.2", las = 1) #empty PCA Plot
abline(v = 0, lty = 3) #plots origin line verticle
abline(h = 0, lty = 3) #plots origin line horizontal
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$site_code, col = RivTypCols) #adds text for samples
text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$site_id, col = "red", cex = 0.6) #adds text for samples
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$ftype, col = "darkblue") #adds text for samples
#text(allYrs_pca$scores[, 1], allYrs_pca$scores[, 2], allYrs$fmodified, col = "forestgreen") #adds text for samples
vec = envfit(allYrs_pca, na.omit(allYrs[, c(94:106)]), choices = c(1, 2)) #computes the environmental variable vectors
plot(vec, col = "blue", cex = 1) #plots vectors (environmetal variables)
pc1_scores <- as.data.frame(allYrs_pca$scores[, 1]*-1) # reverse signs to help with interpretation later
pc2_scores <- as.data.frame(allYrs_pca$scores[, 2]) # reverse signs to help with interpretation later
colnames(pc1_scores) <- c("PC_axis1")
colnames(pc2_scores) <- c("PC_axis2")
pc1_scores$ID <- rownames(pc1_scores)
pc2_scores$ID <- rownames(pc2_scores)
allYrs$ID <- rownames(allYrs)
allYrs <- dplyr::left_join(allYrs, pc1_scores, by = "ID")
allYrs <- dplyr::left_join(allYrs, pc2_scores, by = "ID")
allYrs <- subset(allYrs, select = -c(ID)) # remove ID variable

tiff(filename = "Plots/LT_slopeDistributions_drivers.tiff", width = 6, height = 8, units = 'in', res = 600, compression = 'lzw')

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,1))

##################### Extra metrics ##############################
#### Flow #####
# none
flow <- subset(response_gls, Response == "flow")
flow <- flow$estimate[!is.na(flow$estimate)]
mean(allYrs$flow, na.rm = T)
ave_flow <- mean(allYrs$flow, na.rm = T)
percChange_perYr <- (10^flow-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
flow_Est <- subset(Ests, Response == "flow")
stand_flow <- lapply(flow_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (3.2/4*(b-a)+a)
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
legend(x=-16, y=(4/4*(b-a)+a), legend=expression(paste("a, Flow (",m^3, sep = "", s^-1,")")), bty="n", cex=1.3)
legend(x=5, y=(4/4*(b-a)+a), legend=expression(paste("-0.66%",y^-1, sep = "")), bty="n", cex=1.3)

#### Temperature #####
# none
temp <- subset(response_gls, Response == "temp")
temp <- temp$estimate[!is.na(temp$estimate)]
mean(allYrs$temp, na.rm = T)
ave_temp <- mean(allYrs$temp, na.rm = T)
percChange_perYr <- (temp/ave_temp)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-2
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
yy <- (2.2/4*(b-a)+a)
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
legend(x=-16, y=(3/4*(b-a)+a), legend=("b, Temperature (°C)"), bty="n", cex=1.3)
legend(x=5, y=(3/4*(b-a)+a), legend=expression(paste("+0.44% ",y^-1, sep = "")), bty="n", cex=1.3)

#### PC_axis1 #####
# none
PC_axis1 <- subset(response_gls, Response == "PC_axis1")
PC_axis1 <- PC_axis1$estimate[!is.na(PC_axis1$estimate)]
mean(allYrs$PC_axis1, na.rm = T)
ave_PC_axis1 <- mean(allYrs$PC_axis1, na.rm = T)
percChange_perYr <- ((PC_axis1*2)/sqrt(ave_PC_axis1+3.2))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
PC_axis1_Est <- subset(Ests, Response == "PC_axis1")
stand_PC_axis1 <- lapply(PC_axis1_Est[,2:11], function(x) ((x*2) / sqrt(ave_PC_axis1+3.2))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (1.2/4*(b-a)+a)
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
legend(x=-16, y=(2/4*(b-a)+a), legend=("c, Princ. comp. axis 1"), bty="n", cex=1.3)
legend(x=5, y=(2/4*(b-a)+a), legend=expression(paste("+0.41% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Standardised func. dispersion #####
# none
PC_axis2 <- subset(response_gls, Response == "PC_axis2")
PC_axis2 <- PC_axis2$estimate[!is.na(PC_axis2$estimate)]
mean(allYrs$PC_axis2, na.rm = T)
ave_PC_axis2 <- mean(allYrs$PC_axis2, na.rm = T)
percChange_perYr<-(PC_axis2)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
PC_axis2_Est <- subset(Ests, Response == "PC_axis2")
stand_PC_axis2 <- lapply(PC_axis2_Est[,2:11],"*",100/ave_PC_axis2)
yy <- (0.2/4*(b-a)+a)
points(x=stand_PC_axis2$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_PC_axis2$Q2.5, stand_PC_axis2$Q2.5, stand_PC_axis2$Q97.5, stand_PC_axis2$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_PC_axis2$Q5, stand_PC_axis2$Q5, stand_PC_axis2$Q95, stand_PC_axis2$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_PC_axis2$Q10, stand_PC_axis2$Q10, stand_PC_axis2$Q90, stand_PC_axis2$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-16, y=(1/4*(b-a)+a), legend=("d, Princ. comp. axis 2"), bty="n", cex=1.3)
legend(x=5, y=(1/4*(b-a)+a), legend=expression(paste("+1.44% ",y^-1, sep = "")), bty="n", cex=1.3)

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
