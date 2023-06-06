# attach data
response_gls <- readRDS("Outputs/glsTrends_site_level.rds")
head(response_gls)
unique(response_gls$Response)
length(unique(response_gls$site_id))

Ests <- read.csv("Outputs/LT_Yr_metaanaly_Ests.csv")
Ests <- Ests[, -1]

d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

tiff(filename = "Plots/LT_slopeDistributions_Extra.tiff", width = 6, height = 8, units = 'in', res = 600, compression = 'lzw')

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,1))

##################### Extra metrics ##############################
#### Rarified spp richness #####
# sqrt(x)
rare_SppRich <- subset(response_gls, Response == "spp_rich_rare")
rare_SppRich <- rare_SppRich$estimate[!is.na(rare_SppRich$estimate)]
mean(allYrs$spp_rich_rare, na.rm = T)
ave_rare_SppRich <- mean(allYrs$spp_rich_rare, na.rm = T)
percChange_perYr <- ((rare_SppRich*2)/sqrt(ave_rare_SppRich))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-100,100),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
rare_SppRich_Est <- subset(Ests, Response == "spp_rich_rare")
stand_rare_SppRich <- lapply(rare_SppRich_Est[,2:11], function(x) ((x*2) / sqrt(ave_rare_SppRich))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (3.2/4*(b-a)+a)
points(x=stand_rare_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_rare_SppRich$Q2.5, stand_rare_SppRich$Q2.5, stand_rare_SppRich$Q97.5, stand_rare_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_rare_SppRich$Q5, stand_rare_SppRich$Q5, stand_rare_SppRich$Q95, stand_rare_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_rare_SppRich$Q10, stand_rare_SppRich$Q10, stand_rare_SppRich$Q90, stand_rare_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-115, y=(4/4*(b-a)+a), legend=("a, Rarified taxon richness"), bty="n", cex=1.3)
legend(x=42, y=(4/4*(b-a)+a), legend=expression(paste("+1.56% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Standardised func. richness #####
# none
FRic_SES <- subset(response_gls, Response == "FRic.SES")
FRic_SES <- FRic_SES$estimate[!is.na(FRic_SES$estimate)]
mean(allYrs$FRic.SES, na.rm = T)
ave_FRic_SES <- mean(allYrs$FRic.SES, na.rm = T)*-1
percChange_perYr <- (FRic_SES/ave_FRic_SES)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-100,100),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
FRic_SES_Est <- subset(Ests, Response == "FRic.SES")
stand_FRic_SES <- lapply(FRic_SES_Est[,2:11],"*",100/ave_FRic_SES)
yy <- (2.2/4*(b-a)+a)
points(x=stand_FRic_SES$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FRic_SES$Q2.5, stand_FRic_SES$Q2.5, stand_FRic_SES$Q97.5, stand_FRic_SES$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FRic_SES$Q5, stand_FRic_SES$Q5, stand_FRic_SES$Q95, stand_FRic_SES$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FRic_SES$Q10, stand_FRic_SES$Q10, stand_FRic_SES$Q90, stand_FRic_SES$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-115, y=(3/4*(b-a)+a), legend=("b, Standardised func. richness"), bty="n", cex=1.3)
legend(x=42, y=(3/4*(b-a)+a), legend=expression(paste("+1.95% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Standardised func. evenness #####
# none
FEve_SES <- subset(response_gls, Response == "FEve.SES")
FEve_SES <- FEve_SES$estimate[!is.na(FEve_SES$estimate)]
mean(allYrs$FEve.SES)
ave_FEve_SES <- mean(allYrs$FEve.SES)*-1
percChange_perYr<-(FEve_SES/ave_FEve_SES)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-100,100),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
FEve_SES_Est <- subset(Ests, Response == "FEve.SES")
stand_FEve_SES <- lapply(FEve_SES_Est[,2:11],"*",100/ave_FEve_SES)
yy <- (1.2/4*(b-a)+a)
points(x=stand_FEve_SES$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FEve_SES$Q2.5, stand_FEve_SES$Q2.5, stand_FEve_SES$Q97.5, stand_FEve_SES$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FEve_SES$Q5, stand_FEve_SES$Q5, stand_FEve_SES$Q95, stand_FEve_SES$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FEve_SES$Q10, stand_FEve_SES$Q10, stand_FEve_SES$Q90, stand_FEve_SES$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-115, y=(2/4*(b-a)+a), legend=("c, Standardised func. evenness"), bty="n", cex=1.3)
legend(x=42, y=(2/4*(b-a)+a), legend=expression(paste("-1.71% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Standardised func. dispersion #####
# none
FDis_SES <- subset(response_gls, Response == "FDis.SES")
FDis_SES <- FDis_SES$estimate[!is.na(FDis_SES$estimate)]
mean(allYrs$FDis.SES)
ave_FDis_SES <- mean(allYrs$FDis.SES)*-1
percChange_perYr<-(FDis_SES/ave_FDis_SES)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-100,100),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
FDis_SES_Est <- subset(Ests, Response == "FDis.SES")
stand_FDis_SES <- lapply(FDis_SES_Est[,2:11],"*",100/ave_FDis_SES)
yy <- (0.2/4*(b-a)+a)
points(x=stand_FDis_SES$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FDis_SES$Q2.5, stand_FDis_SES$Q2.5, stand_FDis_SES$Q97.5, stand_FDis_SES$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FDis_SES$Q5, stand_FDis_SES$Q5, stand_FDis_SES$Q95, stand_FDis_SES$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FDis_SES$Q10, stand_FDis_SES$Q10, stand_FDis_SES$Q90, stand_FDis_SES$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-115, y=(1/4*(b-a)+a), legend=("d, Standardised func. dispersion"), bty="n", cex=1.3)
legend(x=42, y=(1/4*(b-a)+a), legend=expression(paste("-1.31% ",y^-1, sep = "")), bty="n", cex=1.3)

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
