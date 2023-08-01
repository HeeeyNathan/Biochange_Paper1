# attach data
response_gls <- readRDS("Outputs/glsTrends_site_level.rds")
head(response_gls)
unique(response_gls$Response)
length(unique(response_gls$site_id))

Ests <- read.csv("Outputs/LT_Yr_metaanaly_Ests.csv")
Ests <- Ests[, -1]

d1 <- read.csv("Data/LT_siteYr_AllData_wNAs_modified.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

tiff(filename = "Plots/LT_slopeDistributions_TaxoGroups.tiff", width = 12, height = 10, units = 'in', res = 600, compression = 'lzw')

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,2))

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
legend(x=-31.25, y=(4/4*(b-a)+a), legend=("a, EPT richness"), bty="n", cex=1.3)
legend(x=16.8, y=(4/4*(b-a)+a), legend=expression(paste("+2.80% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=-31.25, y=(3.2/4*(b-a)+a), legend=("b, Insect richness"), bty="n", cex=1.3)
legend(x=16.8, y=(3.2/4*(b-a)+a), legend=expression(paste("+2.33% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=-31.25, y=(2.4/4*(b-a)+a), legend=("c, Crustacea richness"), bty="n", cex=1.3)
legend(x=16.8, y=(2.4/4*(b-a)+a), legend=expression(paste("+4.30% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=-31.25, y=(1.6/4*(b-a)+a), legend=("d, Mollusc richness"), bty="n", cex=1.3)
legend(x=16.8, y=(1.6/4*(b-a)+a), legend=expression(paste("+6.54% ",y^-1, sep = "")), bty="n", cex=1.3)

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
legend(x=-31.25, y=(0.8/4*(b-a)+a), legend=("e, Annelid richness"), bty="n", cex=1.3)
legend(x=16.8, y=(0.8/4*(b-a)+a), legend=expression(paste("+3.55% ",y^-1, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

##################### Abundance metrics ##############################
#### EPT abundance #####
# log10(x + 1)
ept_abund <- subset(response_gls, Response == "ept_abundance")
ept_abund <- ept_abund$estimate[!is.na(ept_abund$estimate)]
mean(allYrs$ept_abundance)
ave_ept_abund <- mean(allYrs$ept_abundance)
percChange_perYr <- (10^ept_abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-4
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
ept_abund_Est <- subset(Ests, Response == "ept_abundance")
stand_ept_abund <- lapply(ept_abund_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (3.35/4*(b-a)+a)
points(x=stand_ept_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_ept_abund$Q2.5, stand_ept_abund$Q2.5, stand_ept_abund$Q97.5, stand_ept_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_ept_abund$Q5, stand_ept_abund$Q5, stand_ept_abund$Q95, stand_ept_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_ept_abund$Q10, stand_ept_abund$Q10, stand_ept_abund$Q90, stand_ept_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(4/4*(b-a)+a), legend=("f, EPT abundance"), bty="n", cex=1.3)
legend(x=25, y=(4/4*(b-a)+a), legend=expression(paste("+2.94% ",y^-1, sep = "")), bty="n", cex=1.3)

#### insect abundance #####
# log10(x + 1)
insect_abund <- subset(response_gls, Response == "insect_abundance")
insect_abund <- insect_abund$estimate[!is.na(insect_abund$estimate)]
mean(allYrs$insect_abundance)
ave_insect_abund <- mean(allYrs$insect_abundance)
percChange_perYr <- (10^insect_abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
insect_abund_Est <- subset(Ests, Response == "insect_abundance")
stand_insect_abund <- lapply(insect_abund_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (2.55/4*(b-a)+a)
points(x=stand_insect_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_insect_abund$Q2.5, stand_insect_abund$Q2.5, stand_insect_abund$Q97.5, stand_insect_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_insect_abund$Q5, stand_insect_abund$Q5, stand_insect_abund$Q95, stand_insect_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_insect_abund$Q10, stand_insect_abund$Q10, stand_insect_abund$Q90, stand_insect_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(3.2/4*(b-a)+a), legend=("g, Insect abundance"), bty="n", cex=1.3)
legend(x=25, y=(3.2/4*(b-a)+a), legend=expression(paste("+1.01% ",y^-1, sep = "")), bty="n", cex=1.3)

#### crustacea abundance #####
# log(x + 1)
crustacea_abund <- subset(response_gls, Response == "crustacea_abundance")
crustacea_abund <- crustacea_abund$estimate[!is.na(crustacea_abund$estimate)]
mean(allYrs$crustacea_abundance)
ave_crustacea_abund <- mean(allYrs$crustacea_abundance)
percChange_perYr <- (10^crustacea_abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="crustacea_abund", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
crustacea_abund_Est <- subset(Ests, Response == "crustacea_abundance")
stand_crustacea_abund <- lapply(crustacea_abund_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (1.75/4*(b-a)+a)
points(x=stand_crustacea_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_crustacea_abund$Q2.5, stand_crustacea_abund$Q2.5, stand_crustacea_abund$Q97.5, stand_crustacea_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_crustacea_abund$Q5, stand_crustacea_abund$Q5, stand_crustacea_abund$Q95, stand_crustacea_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_crustacea_abund$Q10, stand_crustacea_abund$Q10, stand_crustacea_abund$Q90, stand_crustacea_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(2.4/4*(b-a)+a), legend=("h, Crustacea abundance"), bty="n", cex=1.3)
legend(x=25, y=(2.4/4*(b-a)+a), legend=expression(paste("+12.40% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Mollusc abundance #####
# log10(x + 1)
mollusc_abund <- subset(response_gls, Response == "mollusc_abundance")
mollusc_abund <- mollusc_abund$estimate[!is.na(mollusc_abund$estimate)]
mean(allYrs$mollusc_abundance)
ave_mollusc_abund <- mean(allYrs$mollusc_abundance)
percChange_perYr <- (10^mollusc_abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-40,60),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "#95ccba", border = "#95ccba", lwd =3)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "#f2cc84", border = "#f2cc84", lwd =2)
mollusc_abund_Est <- subset(Ests, Response == "mollusc_abundance")
stand_mollusc_abund <- lapply(mollusc_abund_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.95/4*(b-a)+a)
points(x=stand_mollusc_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_mollusc_abund$Q2.5, stand_mollusc_abund$Q2.5, stand_mollusc_abund$Q97.5, stand_mollusc_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_mollusc_abund$Q5, stand_mollusc_abund$Q5, stand_mollusc_abund$Q95, stand_mollusc_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_mollusc_abund$Q10, stand_mollusc_abund$Q10, stand_mollusc_abund$Q90, stand_mollusc_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(1.6/4*(b-a)+a), legend=("i, Mollusc abundance"), bty="n", cex=1.3)
legend(x=25, y=(1.6/4*(b-a)+a), legend=expression(paste("+10.03% ",y^-1, sep = "")), bty="n", cex=1.3)

#### Annelid abundance #####
#log10(x + 1)
annelid_abund <- subset(response_gls, Response == "annelid_abundance")
annelid_abund <- annelid_abund$estimate[!is.na(annelid_abund$estimate)]
mean(allYrs$annelid_abundance)
ave_annelid_abund <- mean(allYrs$annelid_abundance)
percChange_perYr <- (10^annelid_abund-1)*100
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
annelid_abund_Est <- subset(Ests, Response == "annelid_abundance")
stand_annelid_abund <- lapply(annelid_abund_Est[,2:11], function(x) (10^x - 1)*100) # this must be the same as how you calculated the percChange_perYr
yy <- (0.15/4*(b-a)+a)
points(x=stand_annelid_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_annelid_abund$Q2.5, stand_annelid_abund$Q2.5, stand_annelid_abund$Q97.5, stand_annelid_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_annelid_abund$Q5, stand_annelid_abund$Q5, stand_annelid_abund$Q95, stand_annelid_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_annelid_abund$Q10, stand_annelid_abund$Q10, stand_annelid_abund$Q90, stand_annelid_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-50, y=(0.8/4*(b-a)+a), legend=("j, Annelid abundance"), bty="n", cex=1.3)
legend(x=25, y=(0.8/4*(b-a)+a), legend=expression(paste("+7.91% ",y^-1, sep = "")), bty="n", cex=1.3)

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
