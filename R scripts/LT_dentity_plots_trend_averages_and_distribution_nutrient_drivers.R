# attach data
response_gls <- readRDS("Outputs/glsTrends_site_level_drivers.rds")
head(response_gls)
unique(response_gls$Response)
length(unique(response_gls$site_id))

Ests <- read.csv("Outputs/LT_Yr_metaanaly_weighted_noRandom_Ests_drivers.csv")
Ests <- Ests[, -1]

d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

tiff(filename = "Plots/LT_slopeDistributions_Nutrients.tiff", width = 12, height = 10, units = 'in', res = 600, compression = 'lzw')
# svg(filename = "Plots/LT_slopeDistributions_TaxoFuncIndices.svg", width = 12, height = 10, bg = "white")

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,2))

##################### Colinear nutirents ##############################
#### Alkalinity #####
# none
alkalinity <- subset(response_gls, Response == "alkalinity")
alkalinity <- alkalinity$estimate[!is.na(alkalinity$estimate)]
mean(allYrs$alkalinity, na.rm = T)
ave_alkalinity <- mean(allYrs$alkalinity, na.rm = T)
percChange_perYr <- (alkalinity/ave_alkalinity)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
alkalinity_Est <- subset(Ests, Response == "alkalinity")
stand_alkalinity <- lapply(alkalinity_Est[,2:11],"*",100/ave_alkalinity) # this must be the same as how you calculated the percChange_perYr
yy <- (3.2/4*(b-a)+a)
points(x=stand_alkalinity$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_alkalinity$Q2.5, stand_alkalinity$Q2.5, stand_alkalinity$Q97.5, stand_alkalinity$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_alkalinity$Q5, stand_alkalinity$Q5, stand_alkalinity$Q95, stand_alkalinity$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_alkalinity$Q10, stand_alkalinity$Q10, stand_alkalinity$Q90, stand_alkalinity$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-22.5, y=(4/4*(b-a)+a), legend=expression(paste("a, Alkalinity (",mmol, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(4/4*(b-a)+a), legend=expression(paste("+0.85% ",y^-1, sep = "")), bty="n", cex=1.3)

#### EC #####
# log10(x + 0.1)
EC <- subset(response_gls, Response == "EC")
EC <- EC$estimate[!is.na(EC$estimate)]
mean(allYrs$EC, na.rm = T)
ave_EC <- mean(allYrs$EC, na.rm = T)
percChange_perYr <- (10^EC-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
EC_Est <- subset(Ests, Response == "EC")
stand_EC <- lapply(EC_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (2.2/4*(b-a)+a)
points(x=stand_EC$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_EC$Q2.5, stand_EC$Q2.5, stand_EC$Q97.5, stand_EC$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_EC$Q5, stand_EC$Q5, stand_EC$Q95, stand_EC$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_EC$Q10, stand_EC$Q10, stand_EC$Q90, stand_EC$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-22.5, y=(3/4*(b-a)+a), legend=expression(paste("b, Electrical conductivity (",µS, sep = ".", cm^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(3/4*(b-a)+a), legend=expression(paste("+0.21% ",y^-1, sep = "")), bty="n", cex=1.3)

#### NO3.N #####
# log10(x + 0.1)
NO3.N <- subset(response_gls, Response == "NO3.N")
NO3.N <- NO3.N$estimate[!is.na(NO3.N$estimate)]
mean(allYrs$NO3.N, na.rm = T)
ave_NO3.N <- mean(allYrs$NO3.N, na.rm = T)
percChange_perYr <- (10^NO3.N-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
NO3.N_Est <- subset(Ests, Response == "NO3.N")
stand_NO3.N <- lapply(NO3.N_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (1.2/4*(b-a)+a)
points(x=stand_NO3.N$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_NO3.N$Q2.5, stand_NO3.N$Q2.5, stand_NO3.N$Q97.5, stand_NO3.N$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_NO3.N$Q5, stand_NO3.N$Q5, stand_NO3.N$Q95, stand_NO3.N$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_NO3.N$Q10, stand_NO3.N$Q10, stand_NO3.N$Q90, stand_NO3.N$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-22.5, y=(2/4*(b-a)+a), legend=expression(paste("c, Nitrate (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(2/4*(b-a)+a), legend=expression(paste("+4.40% ",y^-1, sep = "")), bty="n", cex=1.3)

#### NO2.N #####
# log10(x + 0.1)
NO2.N <- subset(response_gls, Response == "NO2.N")
NO2.N <- NO2.N$estimate[!is.na(NO2.N$estimate)]
mean(allYrs$NO2.N, na.rm = T)
ave_NO2.N <- mean(allYrs$NO2.N, na.rm = T)
percChange_perYr <- (10^NO2.N-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
NO2.N_Est <- subset(Ests, Response == "NO2.N")
stand_NO2.N <- lapply(NO2.N_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.2/4*(b-a)+a)
points(x=stand_NO2.N$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_NO2.N$Q2.5, stand_NO2.N$Q2.5, stand_NO2.N$Q97.5, stand_NO2.N$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_NO2.N$Q5, stand_NO2.N$Q5, stand_NO2.N$Q95, stand_NO2.N$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_NO2.N$Q10, stand_NO2.N$Q10, stand_NO2.N$Q90, stand_NO2.N$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-22.5, y=(1/4*(b-a)+a), legend=expression(paste("d, Nitrite (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(1/4*(b-a)+a), legend=expression(paste("-0.06% ",y^-1, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

#### mineral.N #####
# log10(x + 0.1)
mineral.N <- subset(response_gls, Response == "mineral.N")
mineral.N <- mineral.N$estimate[!is.na(mineral.N$estimate)]
mean(allYrs$mineral.N, na.rm = T)
ave_mineral.N <- mean(allYrs$mineral.N, na.rm = T)
percChange_perYr <- (10^mineral.N-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
mineral.N_Est <- subset(Ests, Response == "mineral.N")
stand_mineral.N <- lapply(mineral.N_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (3.2/4*(b-a)+a)
points(x=stand_mineral.N$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_mineral.N$Q2.5, stand_mineral.N$Q2.5, stand_mineral.N$Q97.5, stand_mineral.N$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_mineral.N$Q5, stand_mineral.N$Q5, stand_mineral.N$Q95, stand_mineral.N$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_mineral.N$Q10, stand_mineral.N$Q10, stand_mineral.N$Q90, stand_mineral.N$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-22.5, y=(4/4*(b-a)+a), legend=expression(paste("e, Mineralized Nitrogen (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(4/4*(b-a)+a), legend=expression(paste("+3.53% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Tot.N #####
# log10(x + 0.1)
Tot.N <- subset(response_gls, Response == "Tot.N")
Tot.N <- Tot.N$estimate[!is.na(Tot.N$estimate)]
mean(allYrs$Tot.N, na.rm = T)
ave_Tot.N <- mean(allYrs$Tot.N, na.rm = T)
percChange_perYr <- (10^Tot.N-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
Tot.N_Est <- subset(Ests, Response == "Tot.N")
stand_Tot.N <- lapply(Tot.N_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (2.2/4*(b-a)+a)
points(x=stand_Tot.N$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_Tot.N$Q2.5, stand_Tot.N$Q2.5, stand_Tot.N$Q97.5, stand_Tot.N$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_Tot.N$Q5, stand_Tot.N$Q5, stand_Tot.N$Q95, stand_Tot.N$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_Tot.N$Q10, stand_Tot.N$Q10, stand_Tot.N$Q90, stand_Tot.N$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-22.5, y=(3/4*(b-a)+a), legend=expression(paste("f, Total Nitrogen (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(3/4*(b-a)+a), legend=expression(paste("+1.95% ",y^-1, sep = "")), bty="n", cex=1.3)

#### PO4.P #####
# log10(x + 0.1)
PO4.P <- subset(response_gls, Response == "PO4.P")
PO4.P <- PO4.P$estimate[!is.na(PO4.P$estimate)]
mean(allYrs$PO4.P, na.rm = T)
ave_PO4.P <- mean(allYrs$PO4.P, na.rm = T)
percChange_perYr <- (10^PO4.P-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
PO4.P_Est <- subset(Ests, Response == "PO4.P")
stand_PO4.P <- lapply(PO4.P_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (1.2/4*(b-a)+a)
points(x=stand_PO4.P$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_PO4.P$Q2.5, stand_PO4.P$Q2.5, stand_PO4.P$Q97.5, stand_PO4.P$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_PO4.P$Q5, stand_PO4.P$Q5, stand_PO4.P$Q95, stand_PO4.P$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_PO4.P$Q10, stand_PO4.P$Q10, stand_PO4.P$Q90, stand_PO4.P$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-22.5, y=(2/4*(b-a)+a), legend=expression(paste("g, Phosphate (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(2/4*(b-a)+a), legend=expression(paste("-0.36% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Tot.P #####
# log10(x + 0.1)
Tot.P <- subset(response_gls, Response == "Tot.P")
Tot.P <- Tot.P$estimate[!is.na(Tot.P$estimate)]
mean(allYrs$Tot.P, na.rm = T)
ave_Tot.P <- mean(allYrs$Tot.P, na.rm = T)
percChange_perYr <- (10^Tot.P-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Tot.P", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
Tot.P_Est <- subset(Ests, Response == "Tot.P")
stand_Tot.P <- lapply(Tot.P_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.2/4*(b-a)+a)
points(x=stand_Tot.P$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_Tot.P$Q2.5, stand_Tot.P$Q2.5, stand_Tot.P$Q97.5, stand_Tot.P$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_Tot.P$Q5, stand_Tot.P$Q5, stand_Tot.P$Q95, stand_Tot.P$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_Tot.P$Q10, stand_Tot.P$Q10, stand_Tot.P$Q90, stand_Tot.P$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-22.5, y=(1/4*(b-a)+a), legend=expression(paste("h, Total Phosphorus (",mg, sep = ".", l^-1,")")), bty="n", cex=1.3)
legend(x=7.5, y=(1/4*(b-a)+a), legend=expression(paste("-0.55% ",y^-1, sep = "")), bty="n", cex=1.3)

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
