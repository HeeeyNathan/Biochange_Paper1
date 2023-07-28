# attach data
response_gls <- readRDS("Sensitivity/Outputs/glsTrends_site_level_sensitivity.rds")
head(response_gls)
unique(response_gls$Response)
length(unique(response_gls$site_id))

Ests <- read.csv("Sensitivity/Outputs/LT_Yr_metaanaly_Ests_sensitivity.csv")
Ests <- Ests[, -1]

allYrs <- read.csv("Sensitivity/Data/LT_siteYr_AllData_sensitivity.csv", header=T) 

tiff(filename = "Sensitivity/Plots/LT_slopeDistributions_sensitivity.tiff", width = 12, height = 10, units = 'in', res = 600, compression = 'lzw')

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
legend(x=25, y=(4/4*(b-a)+a), legend=expression(paste("+3.49% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=25, y=(3.2/4*(b-a)+a), legend=expression(paste("+2.37% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=25, y=(2.4/4*(b-a)+a), legend=expression(paste("-1.59% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=25, y=(1.6/4*(b-a)+a), legend=expression(paste("+1.70% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=25, y=(0.8/4*(b-a)+a), legend=expression(paste("+0.90% ",y^-1, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

##################### Richness metrics ##############################
#### EPT Spp Richness #####
# sqrt(x)
ept_SppRich <- subset(response_gls, Response == "ept_spp_richness")
ept_SppRich <- ept_SppRich$estimate[!is.na(ept_SppRich$estimate)]
mean(allYrs$ept_spp_richness, na.rm = T)
ave_ept_SppRich <- mean(allYrs$ept_spp_richness, na.rm = T)
percChange_perYr <- ((ept_SppRich*2)/sqrt(ave_ept_SppRich))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-4
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,40),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
ept_SppRich_Est <- subset(Ests, Response == "ept_richness")
stand_ept_SppRich <- lapply(ept_SppRich_Est[,2:11], function(x) ((x*2) / sqrt(ave_ept_SppRich))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (3.35/4*(b-a)+a)
points(x=stand_ept_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_ept_SppRich$Q2.5, stand_ept_SppRich$Q2.5, stand_ept_SppRich$Q97.5, stand_ept_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_ept_SppRich$Q5, stand_ept_SppRich$Q5, stand_ept_SppRich$Q95, stand_ept_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_ept_SppRich$Q10, stand_ept_SppRich$Q10, stand_ept_SppRich$Q90, stand_ept_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-31.25, y=(4/4*(b-a)+a), legend=("f, EPT richness"), bty="n", cex=1.3)
legend(x=16.8, y=(4/4*(b-a)+a), legend=expression(paste("+7.99% ",y^-1, sep = "")), bty="n", cex=1.3)

#### insect SPP rich #####
# sqrt(x)
insect_SppRich <- subset(response_gls, Response == "insect_spp_richness")
insect_SppRich <- insect_SppRich$estimate[!is.na(insect_SppRich$estimate)]
mean(allYrs$insect_spp_richness, na.rm = T)
ave_insect_SppRich <- mean(allYrs$insect_spp_richness, na.rm = T)
percChange_perYr <- ((insect_SppRich*2)/sqrt(ave_insect_SppRich))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-3
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,40),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
insect_SppRich_Est <- subset(Ests, Response == "insect_richness")
stand_insect_SppRich <- lapply(insect_SppRich_Est[,2:11], function(x) ((x*2) / sqrt(ave_insect_SppRich))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (2.55/4*(b-a)+a)
points(x=stand_insect_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_insect_SppRich$Q2.5, stand_insect_SppRich$Q2.5, stand_insect_SppRich$Q97.5, stand_insect_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_insect_SppRich$Q5, stand_insect_SppRich$Q5, stand_insect_SppRich$Q95, stand_insect_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_insect_SppRich$Q10, stand_insect_SppRich$Q10, stand_insect_SppRich$Q90, stand_insect_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-31.25, y=(3.2/4*(b-a)+a), legend=("g, Insect richness"), bty="n", cex=1.3)
legend(x=16.8, y=(3.2/4*(b-a)+a), legend=expression(paste("+2.68% ",y^-1, sep = "")), bty="n", cex=1.3)

#### crustacea richness #####
# sqrt(x)
crustacea_SppRich <- subset(response_gls, Response == "crustacea_spp_richness")
crustacea_SppRich <- crustacea_SppRich$estimate[!is.na(crustacea_SppRich$estimate)]
mean(allYrs$crustacea_spp_richness)
ave_crustacea_SppRich <- mean(allYrs$crustacea_spp_richness)
percChange_perYr <- ((crustacea_SppRich*2)/sqrt(ave_crustacea_SppRich))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-2
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,40),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
crustacea_SppRich_Est <- subset(Ests, Response == "crustacea_richness")
stand_crustacea_SppRich <- lapply(crustacea_SppRich_Est[,2:11], function(x) ((x*2) / sqrt(ave_crustacea_SppRich))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (1.75/4*(b-a)+a)
points(x=stand_crustacea_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_crustacea_SppRich$Q2.5, stand_crustacea_SppRich$Q2.5, stand_crustacea_SppRich$Q97.5, stand_crustacea_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_crustacea_SppRich$Q5, stand_crustacea_SppRich$Q5, stand_crustacea_SppRich$Q95, stand_crustacea_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_crustacea_SppRich$Q10, stand_crustacea_SppRich$Q10, stand_crustacea_SppRich$Q90, stand_crustacea_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-31.25, y=(2.4/4*(b-a)+a), legend=("h, Crustacea richness"), bty="n", cex=1.3)
legend(x=16.8, y=(2.4/4*(b-a)+a), legend=expression(paste("+4.61% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Mollusc richness #####
# sqrt(x)
mollusc_SppRich <- subset(response_gls, Response == "mollusc_spp_richness")
mollusc_SppRich <- mollusc_SppRich$estimate[!is.na(mollusc_SppRich$estimate)]
mean(allYrs$mollusc_spp_richness)
ave_mollusc_SppRich <- mean(allYrs$mollusc_spp_richness)
percChange_perYr <- ((mollusc_SppRich*2)/sqrt(ave_mollusc_SppRich))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,40),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
mollusc_SppRich_Est <- subset(Ests, Response == "mollusc_richness")
stand_mollusc_SppRich <- lapply(mollusc_SppRich_Est[,2:11], function(x) ((x*2) / sqrt(ave_mollusc_SppRich))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.95/4*(b-a)+a)
points(x=stand_mollusc_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_mollusc_SppRich$Q2.5, stand_mollusc_SppRich$Q2.5, stand_mollusc_SppRich$Q97.5, stand_mollusc_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_mollusc_SppRich$Q5, stand_mollusc_SppRich$Q5, stand_mollusc_SppRich$Q95, stand_mollusc_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_mollusc_SppRich$Q10, stand_mollusc_SppRich$Q10, stand_mollusc_SppRich$Q90, stand_mollusc_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-31.25, y=(1.6/4*(b-a)+a), legend=("i, Mollusc richness"), bty="n", cex=1.3)
legend(x=16.8, y=(1.6/4*(b-a)+a), legend=expression(paste("+5.08% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Annelid richness #####
# sqrt(x)
annelid_SppRich <- subset(response_gls, Response == "annelid_spp_richness")
annelid_SppRich <- annelid_SppRich$estimate[!is.na(annelid_SppRich$estimate)]
mean(allYrs$annelid_spp_richness)
ave_annelid_SppRich <- mean(allYrs$annelid_spp_richness)
percChange_perYr <- ((annelid_SppRich*2)/sqrt(ave_annelid_SppRich))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-25,40),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
annelid_SppRich_Est <- subset(Ests, Response == "annelid_richness")
stand_annelid_SppRich <- lapply(annelid_SppRich_Est[,2:11], function(x) ((x*2) / sqrt(ave_annelid_SppRich))*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.15/4*(b-a)+a)
points(x=stand_annelid_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_annelid_SppRich$Q2.5, stand_annelid_SppRich$Q2.5, stand_annelid_SppRich$Q97.5, stand_annelid_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_annelid_SppRich$Q5, stand_annelid_SppRich$Q5, stand_annelid_SppRich$Q95, stand_annelid_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_annelid_SppRich$Q10, stand_annelid_SppRich$Q10, stand_annelid_SppRich$Q90, stand_annelid_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-31.25, y=(0.8/4*(b-a)+a), legend=("j, Annelid richness"), bty="n", cex=1.3)
legend(x=16.8, y=(0.8/4*(b-a)+a), legend=expression(paste("+3.52% ",y^-1, sep = "")), bty="n", cex=1.3)

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
