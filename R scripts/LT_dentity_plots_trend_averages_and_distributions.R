# attach data
response_gls <- readRDS("Outputs/glsTrends_site_level.rds")
head(response_gls)
unique(response_gls$Response)
length(unique(response_gls$site_id))

Ests <- read.csv("Outputs/LT_Yr_metaanaly_Ests.csv")
Ests <- Ests[, -1]

d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

tiff(filename = "Plots/LT_slopeDistributions_TaxoFuncIndices.tiff", width = 12, height = 10, units = 'in', res = 600, compression = 'lzw')

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,2))

#####################Taxonomic metrics ##############################
#### Abundance #####
# log10(x + 1)
abund <- subset(response_gls, Response == "abundance")
abund <- abund$estimate[!is.na(abund$estimate)]
mean(allYrs$abundance)
ave_Abund <- mean(allYrs$abundance)
percChange_perYr <- (10^abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-4
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
abund_Est <- subset(Ests, Response == "abundance")
stand_abund <- lapply(abund_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (3.35/4*(b-a)+a)
points(x=stand_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_abund$Q2.5, stand_abund$Q2.5, stand_abund$Q97.5, stand_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_abund$Q5, stand_abund$Q5, stand_abund$Q95, stand_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_abund$Q10, stand_abund$Q10, stand_abund$Q90, stand_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(4/4*(b-a)+a), legend=("a, Abundance"), bty="n", cex=1.3)
legend(x=25, y=(4/4*(b-a)+a), legend=expression(paste("+3.53% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Spp Richness #####
# sqrt(x)
SppRich <- subset(response_gls, Response == "spp_richness")
SppRich <- SppRich$estimate[!is.na(SppRich$estimate)]
mean(allYrs$spp_richness)
ave_SppRich <- mean(allYrs$spp_richness)
percChange_perYr <- ((SppRich*2)/sqrt(ave_SppRich))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-3
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
SppRich_Est <- subset(Ests, Response == "spp_richness")
stand_SppRich <- lapply(SppRich_Est[,2:11], function(x) ((x*2) / sqrt(ave_SppRich))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (2.55/4*(b-a)+a)
points(x=stand_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_SppRich$Q2.5, stand_SppRich$Q2.5, stand_SppRich$Q97.5, stand_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRich$Q5, stand_SppRich$Q5, stand_SppRich$Q95, stand_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRich$Q10, stand_SppRich$Q10, stand_SppRich$Q90, stand_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-50, y=(3.2/4*(b-a)+a), legend=("b, Taxon richness"), bty="n", cex=1.3)
legend(x=25, y=(3.2/4*(b-a)+a), legend=expression(paste("+3.28% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Shannon's Evenness #####
# none
E10 <- subset(response_gls, Response == "E10")
E10 <- E10$estimate[!is.na(E10$estimate)]
mean(allYrs$E10)
ave_E10 <- mean(allYrs$E10)
percChange_perYr <- (E10/ave_E10)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-2
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
E10_Est <- subset(Ests, Response == "E10")
stand_E10 <- lapply(E10_Est[,2:11],"*",100/ave_E10)
yy <- (1.75/4*(b-a)+a)
points(x=stand_E10$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_E10$Q2.5, stand_E10$Q2.5, stand_E10$Q97.5, stand_E10$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_E10$Q5, stand_E10$Q5, stand_E10$Q95, stand_E10$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_E10$Q10, stand_E10$Q10, stand_E10$Q90, stand_E10$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(2.4/4*(b-a)+a), legend=("c, Evenness"), bty="n", cex=1.3)
legend(x=25, y=(2.4/4*(b-a)+a), legend=expression(paste("-1.90% ",y^-1, sep = "")), bty="n", cex=1.3)

#### ShannonsH #####
# none
ShannonsH <- subset(response_gls, Response == "shannonsH")
ShannonsH <- ShannonsH$estimate[!is.na(ShannonsH$estimate)]
mean(allYrs$shannonsH)
ave_ShannonsH <- mean(allYrs$shannonsH)
percChange_perYr <- (ShannonsH/ave_ShannonsH)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
ShannonsH_Est <- subset(Ests, Response == "shannonsH")
stand_ShannonsH <- lapply(ShannonsH_Est[,2:11],"*",100/ave_ShannonsH)
yy <- (0.95/4*(b-a)+a)
points(x=stand_ShannonsH$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_ShannonsH$Q2.5, stand_ShannonsH$Q2.5, stand_ShannonsH$Q97.5, stand_ShannonsH$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_ShannonsH$Q5, stand_ShannonsH$Q5, stand_ShannonsH$Q95, stand_ShannonsH$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_ShannonsH$Q10, stand_ShannonsH$Q10, stand_ShannonsH$Q90, stand_ShannonsH$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(1.6/4*(b-a)+a), legend=("d, Shannon's H"), bty="n", cex=1.3)
legend(x=25, y=(1.6/4*(b-a)+a), legend=expression(paste("+1.71% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Turnover #####
# x^2
turnover <- subset(response_gls, Response == "turnover")
turnover <- turnover$estimate[!is.na(turnover$estimate)]
mean(allYrs$turnover, na.rm = T)
ave_turnover <- mean(allYrs$turnover, na.rm = T)
percChange_perYr<-((turnover/3)/ave_turnover^3)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
turnover_Est <- subset(Ests, Response == "turnover")
stand_turnover <- lapply(turnover_Est[,2:11], function(x) {((x/3) / ave_turnover^3)*100})
yy <- (0.15/4*(b-a)+a)
points(x=stand_turnover$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_turnover$Q2.5, stand_turnover$Q2.5, stand_turnover$Q97.5, stand_turnover$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_turnover$Q5, stand_turnover$Q5, stand_turnover$Q95, stand_turnover$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_turnover$Q10, stand_turnover$Q10, stand_turnover$Q90, stand_turnover$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(0.8/4*(b-a)+a), legend=("e, Turnover"), bty="n", cex=1.3)
legend(x=25, y=(0.8/4*(b-a)+a), legend=expression(paste("+0.77% ",y^-1, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

#####################Functional metrics ##############################
#### Functional Redundancy #####
# none
FRed <- subset(response_gls, Response == "FRed")
FRed <- FRed$estimate[!is.na(FRed$estimate)]
mean(allYrs$FRed)
ave_FRed <- mean(allYrs$FRed)
percChange_perYr<-(FRed/ave_FRed)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-4
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
FRed_Est <- subset(Ests, Response == "FRed")
stand_FRed <- lapply(FRed_Est[,2:11],"*",100/ave_FRed)
yy <- (3.35/4*(b-a)+a)
points(x=stand_FRed$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FRed$Q2.5, stand_FRed$Q2.5, stand_FRed$Q97.5, stand_FRed$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FRed$Q5, stand_FRed$Q5, stand_FRed$Q95, stand_FRed$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FRed$Q10, stand_FRed$Q10, stand_FRed$Q90, stand_FRed$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(4/4*(b-a)+a), legend=("f, Func. redundancy"), bty="n", cex=1.3)
legend(x=25, y=(4/4*(b-a)+a), legend=expression(paste("+0.27% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Functional Richness #####
# log10(x + 0.28)
FRic <- subset(response_gls, Response == "FRic")
FRic <- FRic$estimate[!is.na(FRic$estimate)]
mean(allYrs$FRic)
ave_FRic = mean(allYrs$FRic)
percChange_perYr<-(10^FRic-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-3
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
FRic_Est <- subset(Ests, Response == "FRic")
stand_FRic <- lapply(FRic_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (2.55/4*(b-a)+a)
points(x=stand_FRic$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FRic$Q2.5, stand_FRic$Q2.5, stand_FRic$Q97.5, stand_FRic$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FRic$Q5, stand_FRic$Q5, stand_FRic$Q95, stand_FRic$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FRic$Q10, stand_FRic$Q10, stand_FRic$Q90, stand_FRic$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(3.2/4*(b-a)+a), legend=("g, Func. richness"), bty="n", cex=1.3)
legend(x=25, y=(3.2/4*(b-a)+a), legend=expression(paste("+0.80% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Functional evenness #####
# x^2
FEve <- subset(response_gls, Response == "FEve")
FEve <- FEve$estimate[!is.na(FEve$estimate)]
mean(allYrs$FEve)
ave_FEve <- mean(allYrs$FEve)
percChange_perYr<-((FEve/3)/ave_FEve^3)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="FEve", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
FEve_Est <- subset(Ests, Response == "FEve")
stand_FEve <- lapply(FEve_Est[,2:11], function(x) {((x/3) / ave_FEve^3)*100})
yy <- (1.75/4*(b-a)+a)
points(x=stand_FEve$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FEve$Q2.5, stand_FEve$Q2.5, stand_FEve$Q97.5, stand_FEve$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FEve$Q5, stand_FEve$Q5, stand_FEve$Q95, stand_FEve$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FEve$Q10, stand_FEve$Q10, stand_FEve$Q90, stand_FEve$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(2.4/4*(b-a)+a), legend=("h, Func. evenness"), bty="n", cex=1.3)
legend(x=25, y=(2.4/4*(b-a)+a), legend=expression(paste("+0.08% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Functional dispersion #####
# x^3
FDis <- subset(response_gls, Response == "FDis")
FDis <- FDis$estimate[!is.na(FDis$estimate)]
mean(allYrs$FDis)
ave_FDis <- mean(allYrs$FDis)
percChange_perYr<-((FDis/3)/ave_FDis^3)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
FDis_Est <- subset(Ests, Response == "FDis")
stand_FDis <- lapply(FDis_Est[,2:11], function(x) {((x/3) / ave_FDis^3)*100})
yy <- (0.95/4*(b-a)+a)
points(x=stand_FDis$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FDis$Q2.5, stand_FDis$Q2.5, stand_FDis$Q97.5, stand_FDis$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FDis$Q5, stand_FDis$Q5, stand_FDis$Q95, stand_FDis$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FDis$Q10, stand_FDis$Q10, stand_FDis$Q90, stand_FDis$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(1.6/4*(b-a)+a), legend=("i, Func. dispersion"), bty="n", cex=1.3)
legend(x=25, y=(1.6/4*(b-a)+a), legend=expression(paste("+0.36% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Func turnover #####
# log10(x + 0.01)
F_to <- subset(response_gls, Response == "F_turnover")
F_to <- F_to$estimate[!is.na(F_to$estimate)]
mean(allYrs$F_turnover, na.rm = T)
ave_F_to <- mean(allYrs$F_turnover, na.rm = T)
percChange_perYr<-(10^F_to-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =3)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
F_to_Est <- subset(Ests, Response == "F_turnover")
stand_F_to <- lapply(F_to_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.15/4*(b-a)+a)
points(x=stand_F_to$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_F_to$Q2.5, stand_F_to$Q2.5, stand_F_to$Q97.5, stand_F_to$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_F_to$Q5, stand_F_to$Q5, stand_F_to$Q95, stand_F_to$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_F_to$Q10, stand_F_to$Q10, stand_F_to$Q90, stand_F_to$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(0.8/4*(b-a)+a), legend=("j, Func. turnover"), bty="n", cex=1.3)
legend(x=25, y=(0.8/4*(b-a)+a), legend=expression(paste("-1.10% ",y^-1, sep = "")), bty="n", cex=1.3)

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
