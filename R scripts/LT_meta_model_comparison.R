# attach data
ESTs <- read.csv("Outputs/LT_Yr_Ests_Comparison.csv")
head(ESTs)

tiff(filename = "Plots/LT_EstComparison.tiff", width = 12, height = 10, units = 'in', res = 600, compression = 'lzw')

####layout
layout(mat = matrix(c(1:4), 
                    nrow = 1, 
                    ncol = 4),
       heights = c(1,1,1,1),
       widths = c(1,1,1,1))

par(mar=c(4,0.4,0.4,0.8))
#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

sub <- subset(ESTs, Model == "weighted_random")
yy <- c(13:0)
est <-sub$percEstimate
min(sub$percQ0.5)
max(sub$percQ99.5)
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-16.6,17.2), ylim=c(0,13.5),cex=2, bty="n")
#polygon(x=c(-100,-100,0,0),
#        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
segments(x0=0,y0=0,x1=0,y1=13.3,lty=2, lwd=2,col="grey60")
axis(2, at=yy, lwd = 0, lwd.ticks = 0, labels=c("Taxon richness", "Rarefied taxon richness", "Shannon's H", "Evenness", "Abundance", "Turnover",
                                                "Functional turnover", "Functional richness", "Functional evenness", "Functional dispersion", 
                                                "Stand. functional richness", "Stand. functional dispersion", "Stand. functional eveness", 
                                                "Functional redundancy"), las=1,cex.axis=1.3)
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
legend("top", legend=("a, Weighted Meta-analysis wRand"), bty="n", cex=1.3)
mm <- cbind(yy,est)
mm_neg <- as.data.frame(mm[ which(est < 0),])
mm_pos <- as.data.frame(mm[ which(est > 0),])
#points(yy ~ est, pch="l",cex=2.2,col="white")
points(mm_neg$yy ~ mm_neg$est, pch="l",cex=2,col="firebrick2")
points(mm_pos$yy ~ mm_pos$est, pch="l",cex=2,col="dodgerblue")
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)

polygon(x=c(if (sub$percQ2.5[1]<0) {sub$percQ2.5[1]} else {0}, if (sub$percQ2.5[1]<0) {sub$percQ2.5[1]} else {0}, if (sub$percQ97.5[1]<0) {sub$percQ97.5[1]} else {0}, if (sub$percQ97.5[1]<0) {sub$percQ97.5[1]} else {0}),y=(yyy1+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[1]<0) {sub$percQ5[1]} else {0}, if (sub$percQ5[1]<0) {sub$percQ5[1]} else {0}, if (sub$percQ95[1]<0) {sub$percQ95[1]} else {0}, if (sub$percQ95[1]<0) {sub$percQ95[1]} else {0}),y=(yyy2+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[1]<0) {sub$percQ10[1]} else {0}, if (sub$percQ10[1]<0) {sub$percQ10[1]} else {0}, if (sub$percQ90[1]<0) {sub$percQ90[1]} else {0}, if (sub$percQ90[1]<0) {sub$percQ90[1]} else {0}),y=(yyy3+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[1]>0) {sub$percQ2.5[1]} else {0}, if (sub$percQ2.5[1]>0) {sub$percQ2.5[1]} else {0}, if (sub$percQ97.5[1]>0) {sub$percQ97.5[1]} else {0}, if (sub$percQ97.5[1]>0) {sub$percQ97.5[1]} else {0}),y=(yyy1+12), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[1]>0) {sub$percQ5[1]} else {0}, if (sub$percQ5[1]>0) {sub$percQ5[1]} else {0}, if (sub$percQ95[1]>0) {sub$percQ95[1]} else {0}, if (sub$percQ95[1]>0) {sub$percQ95[1]} else {0}),y=(yyy2+12), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[1]>0) {sub$percQ10[1]} else {0}, if (sub$percQ10[1]>0) {sub$percQ10[1]} else {0}, if (sub$percQ90[1]>0) {sub$percQ90[1]} else {0}, if (sub$percQ90[1]>0) {sub$percQ90[1]} else {0}),y=(yyy3+12), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[2]<0) {sub$percQ2.5[2]} else {0}, if (sub$percQ2.5[2]<0) {sub$percQ2.5[2]} else {0}, if (sub$percQ97.5[2]<0) {sub$percQ97.5[2]} else {0}, if (sub$percQ97.5[2]<0) {sub$percQ97.5[2]} else {0}),y=(yyy1+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[2]<0) {sub$percQ5[2]} else {0}, if (sub$percQ5[2]<0) {sub$percQ5[2]} else {0}, if (sub$percQ95[2]<0) {sub$percQ95[2]} else {0}, if (sub$percQ95[2]<0) {sub$percQ95[2]} else {0}),y=(yyy2+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[2]<0) {sub$percQ10[2]} else {0}, if (sub$percQ10[2]<0) {sub$percQ10[2]} else {0}, if (sub$percQ90[2]<0) {sub$percQ90[2]} else {0}, if (sub$percQ90[2]<0) {sub$percQ90[2]} else {0}),y=(yyy3+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[2]>0) {sub$percQ2.5[2]} else {0}, if (sub$percQ2.5[2]>0) {sub$percQ2.5[2]} else {0}, if (sub$percQ97.5[2]>0) {sub$percQ97.5[2]} else {0}, if (sub$percQ97.5[2]>0) {sub$percQ97.5[2]} else {0}),y=(yyy1+11), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[2]>0) {sub$percQ5[2]} else {0}, if (sub$percQ5[2]>0) {sub$percQ5[2]} else {0}, if (sub$percQ95[2]>0) {sub$percQ95[2]} else {0}, if (sub$percQ95[2]>0) {sub$percQ95[2]} else {0}),y=(yyy2+11), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[2]>0) {sub$percQ10[2]} else {0}, if (sub$percQ10[2]>0) {sub$percQ10[2]} else {0}, if (sub$percQ90[2]>0) {sub$percQ90[2]} else {0}, if (sub$percQ90[2]>0) {sub$percQ90[2]} else {0}),y=(yyy3+11), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[3]<0) {sub$percQ2.5[3]} else {0}, if (sub$percQ2.5[3]<0) {sub$percQ2.5[3]} else {0}, if (sub$percQ97.5[3]<0) {sub$percQ97.5[3]} else {0}, if (sub$percQ97.5[3]<0) {sub$percQ97.5[3]} else {0}),y=(yyy1+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[3]<0) {sub$percQ5[3]} else {0}, if (sub$percQ5[3]<0) {sub$percQ5[3]} else {0}, if (sub$percQ95[3]<0) {sub$percQ95[3]} else {0}, if (sub$percQ95[3]<0) {sub$percQ95[3]} else {0}),y=(yyy2+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[3]<0) {sub$percQ10[3]} else {0}, if (sub$percQ10[3]<0) {sub$percQ10[3]} else {0}, if (sub$percQ90[3]<0) {sub$percQ90[3]} else {0}, if (sub$percQ90[3]<0) {sub$percQ90[3]} else {0}),y=(yyy3+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[3]>0) {sub$percQ2.5[3]} else {0}, if (sub$percQ2.5[3]>0) {sub$percQ2.5[3]} else {0}, if (sub$percQ97.5[3]>0) {sub$percQ97.5[3]} else {0}, if (sub$percQ97.5[3]>0) {sub$percQ97.5[3]} else {0}),y=(yyy1+10), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[3]>0) {sub$percQ5[3]} else {0}, if (sub$percQ5[3]>0) {sub$percQ5[3]} else {0}, if (sub$percQ95[3]>0) {sub$percQ95[3]} else {0}, if (sub$percQ95[3]>0) {sub$percQ95[3]} else {0}),y=(yyy2+10), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[3]>0) {sub$percQ10[3]} else {0}, if (sub$percQ10[3]>0) {sub$percQ10[3]} else {0}, if (sub$percQ90[3]>0) {sub$percQ90[3]} else {0}, if (sub$percQ90[3]>0) {sub$percQ90[3]} else {0}),y=(yyy3+10), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[4]<0) {sub$percQ2.5[4]} else {0}, if (sub$percQ2.5[4]<0) {sub$percQ2.5[4]} else {0}, if (sub$percQ97.5[4]<0) {sub$percQ97.5[4]} else {0}, if (sub$percQ97.5[4]<0) {sub$percQ97.5[4]} else {0}),y=(yyy1+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[4]<0) {sub$percQ5[4]} else {0}, if (sub$percQ5[4]<0) {sub$percQ5[4]} else {0}, if (sub$percQ95[4]<0) {sub$percQ95[4]} else {0}, if (sub$percQ95[4]<0) {sub$percQ95[4]} else {0}),y=(yyy2+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[4]<0) {sub$percQ10[4]} else {0}, if (sub$percQ10[4]<0) {sub$percQ10[4]} else {0}, if (sub$percQ90[4]<0) {sub$percQ90[4]} else {0}, if (sub$percQ90[4]<0) {sub$percQ90[4]} else {0}),y=(yyy3+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[4]>0) {sub$percQ2.5[4]} else {0}, if (sub$percQ2.5[4]>0) {sub$percQ2.5[4]} else {0}, if (sub$percQ97.5[4]>0) {sub$percQ97.5[4]} else {0}, if (sub$percQ97.5[4]>0) {sub$percQ97.5[4]} else {0}),y=(yyy1+9), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[4]>0) {sub$percQ5[4]} else {0}, if (sub$percQ5[4]>0) {sub$percQ5[4]} else {0}, if (sub$percQ95[4]>0) {sub$percQ95[4]} else {0}, if (sub$percQ95[4]>0) {sub$percQ95[4]} else {0}),y=(yyy2+9), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[4]>0) {sub$percQ10[4]} else {0}, if (sub$percQ10[4]>0) {sub$percQ10[4]} else {0}, if (sub$percQ90[4]>0) {sub$percQ90[4]} else {0}, if (sub$percQ90[4]>0) {sub$percQ90[4]} else {0}),y=(yyy3+9), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[5]<0) {sub$percQ2.5[5]} else {0}, if (sub$percQ2.5[5]<0) {sub$percQ2.5[5]} else {0}, if (sub$percQ97.5[5]<0) {sub$percQ97.5[5]} else {0}, if (sub$percQ97.5[5]<0) {sub$percQ97.5[5]} else {0}),y=(yyy1+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[5]<0) {sub$percQ5[5]} else {0}, if (sub$percQ5[5]<0) {sub$percQ5[5]} else {0}, if (sub$percQ95[5]<0) {sub$percQ95[5]} else {0}, if (sub$percQ95[5]<0) {sub$percQ95[5]} else {0}),y=(yyy2+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[5]<0) {sub$percQ10[5]} else {0}, if (sub$percQ10[5]<0) {sub$percQ10[5]} else {0}, if (sub$percQ90[5]<0) {sub$percQ90[5]} else {0}, if (sub$percQ90[5]<0) {sub$percQ90[5]} else {0}),y=(yyy3+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[5]>0) {sub$percQ2.5[5]} else {0}, if (sub$percQ2.5[5]>0) {sub$percQ2.5[5]} else {0}, if (sub$percQ97.5[5]>0) {sub$percQ97.5[5]} else {0}, if (sub$percQ97.5[5]>0) {sub$percQ97.5[5]} else {0}),y=(yyy1+8), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[5]>0) {sub$percQ5[5]} else {0}, if (sub$percQ5[5]>0) {sub$percQ5[5]} else {0}, if (sub$percQ95[5]>0) {sub$percQ95[5]} else {0}, if (sub$percQ95[5]>0) {sub$percQ95[5]} else {0}),y=(yyy2+8), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[5]>0) {sub$percQ10[5]} else {0}, if (sub$percQ10[5]>0) {sub$percQ10[5]} else {0}, if (sub$percQ90[5]>0) {sub$percQ90[5]} else {0}, if (sub$percQ90[5]>0) {sub$percQ90[5]} else {0}),y=(yyy3+8), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[6]<0) {sub$percQ2.5[6]} else {0}, if (sub$percQ2.5[6]<0) {sub$percQ2.5[6]} else {0}, if (sub$percQ97.5[6]<0) {sub$percQ97.5[6]} else {0}, if (sub$percQ97.5[6]<0) {sub$percQ97.5[6]} else {0}),y=(yyy1+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[6]<0) {sub$percQ5[6]} else {0}, if (sub$percQ5[6]<0) {sub$percQ5[6]} else {0}, if (sub$percQ95[6]<0) {sub$percQ95[6]} else {0}, if (sub$percQ95[6]<0) {sub$percQ95[6]} else {0}),y=(yyy2+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[6]<0) {sub$percQ10[6]} else {0}, if (sub$percQ10[6]<0) {sub$percQ10[6]} else {0}, if (sub$percQ90[6]<0) {sub$percQ90[6]} else {0}, if (sub$percQ90[6]<0) {sub$percQ90[6]} else {0}),y=(yyy3+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[6]>0) {sub$percQ2.5[6]} else {0}, if (sub$percQ2.5[6]>0) {sub$percQ2.5[6]} else {0}, if (sub$percQ97.5[6]>0) {sub$percQ97.5[6]} else {0}, if (sub$percQ97.5[6]>0) {sub$percQ97.5[6]} else {0}),y=(yyy1+7), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[6]>0) {sub$percQ5[6]} else {0}, if (sub$percQ5[6]>0) {sub$percQ5[6]} else {0}, if (sub$percQ95[6]>0) {sub$percQ95[6]} else {0}, if (sub$percQ95[6]>0) {sub$percQ95[6]} else {0}),y=(yyy2+7), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[6]>0) {sub$percQ10[6]} else {0}, if (sub$percQ10[6]>0) {sub$percQ10[6]} else {0}, if (sub$percQ90[6]>0) {sub$percQ90[6]} else {0}, if (sub$percQ90[6]>0) {sub$percQ90[6]} else {0}),y=(yyy3+7), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[7]<0) {sub$percQ2.5[7]} else {0}, if (sub$percQ2.5[7]<0) {sub$percQ2.5[7]} else {0}, if (sub$percQ97.5[7]<0) {sub$percQ97.5[7]} else {0}, if (sub$percQ97.5[7]<0) {sub$percQ97.5[7]} else {0}),y=(yyy1+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[7]<0) {sub$percQ5[7]} else {0}, if (sub$percQ5[7]<0) {sub$percQ5[7]} else {0}, if (sub$percQ95[7]<0) {sub$percQ95[7]} else {0}, if (sub$percQ95[7]<0) {sub$percQ95[7]} else {0}),y=(yyy2+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[7]<0) {sub$percQ10[7]} else {0}, if (sub$percQ10[7]<0) {sub$percQ10[7]} else {0}, if (sub$percQ90[7]<0) {sub$percQ90[7]} else {0}, if (sub$percQ90[7]<0) {sub$percQ90[7]} else {0}),y=(yyy3+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[7]>0) {sub$percQ2.5[7]} else {0}, if (sub$percQ2.5[7]>0) {sub$percQ2.5[7]} else {0}, if (sub$percQ97.5[7]>0) {sub$percQ97.5[7]} else {0}, if (sub$percQ97.5[7]>0) {sub$percQ97.5[7]} else {0}),y=(yyy1+6), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[7]>0) {sub$percQ5[7]} else {0}, if (sub$percQ5[7]>0) {sub$percQ5[7]} else {0}, if (sub$percQ95[7]>0) {sub$percQ95[7]} else {0}, if (sub$percQ95[7]>0) {sub$percQ95[7]} else {0}),y=(yyy2+6), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[7]>0) {sub$percQ10[7]} else {0}, if (sub$percQ10[7]>0) {sub$percQ10[7]} else {0}, if (sub$percQ90[7]>0) {sub$percQ90[7]} else {0}, if (sub$percQ90[7]>0) {sub$percQ90[7]} else {0}),y=(yyy3+6), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[8]<0) {sub$percQ2.5[8]} else {0}, if (sub$percQ2.5[8]<0) {sub$percQ2.5[8]} else {0}, if (sub$percQ97.5[8]<0) {sub$percQ97.5[8]} else {0}, if (sub$percQ97.5[8]<0) {sub$percQ97.5[8]} else {0}),y=(yyy1+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[8]<0) {sub$percQ5[8]} else {0}, if (sub$percQ5[8]<0) {sub$percQ5[8]} else {0}, if (sub$percQ95[8]<0) {sub$percQ95[8]} else {0}, if (sub$percQ95[8]<0) {sub$percQ95[8]} else {0}),y=(yyy2+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[8]<0) {sub$percQ10[8]} else {0}, if (sub$percQ10[8]<0) {sub$percQ10[8]} else {0}, if (sub$percQ90[8]<0) {sub$percQ90[8]} else {0}, if (sub$percQ90[8]<0) {sub$percQ90[8]} else {0}),y=(yyy3+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[8]>0) {sub$percQ2.5[8]} else {0}, if (sub$percQ2.5[8]>0) {sub$percQ2.5[8]} else {0}, if (sub$percQ97.5[8]>0) {sub$percQ97.5[8]} else {0}, if (sub$percQ97.5[8]>0) {sub$percQ97.5[8]} else {0}),y=(yyy1+5), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[8]>0) {sub$percQ5[8]} else {0}, if (sub$percQ5[8]>0) {sub$percQ5[8]} else {0}, if (sub$percQ95[8]>0) {sub$percQ95[8]} else {0}, if (sub$percQ95[8]>0) {sub$percQ95[8]} else {0}),y=(yyy2+5), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[8]>0) {sub$percQ10[8]} else {0}, if (sub$percQ10[8]>0) {sub$percQ10[8]} else {0}, if (sub$percQ90[8]>0) {sub$percQ90[8]} else {0}, if (sub$percQ90[8]>0) {sub$percQ90[8]} else {0}),y=(yyy3+5), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[9]<0) {sub$percQ2.5[9]} else {0}, if (sub$percQ2.5[9]<0) {sub$percQ2.5[9]} else {0}, if (sub$percQ97.5[9]<0) {sub$percQ97.5[9]} else {0}, if (sub$percQ97.5[9]<0) {sub$percQ97.5[9]} else {0}),y=(yyy1+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[9]<0) {sub$percQ5[9]} else {0}, if (sub$percQ5[9]<0) {sub$percQ5[9]} else {0}, if (sub$percQ95[9]<0) {sub$percQ95[9]} else {0}, if (sub$percQ95[9]<0) {sub$percQ95[9]} else {0}),y=(yyy2+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[9]<0) {sub$percQ10[9]} else {0}, if (sub$percQ10[9]<0) {sub$percQ10[9]} else {0}, if (sub$percQ90[9]<0) {sub$percQ90[9]} else {0}, if (sub$percQ90[9]<0) {sub$percQ90[9]} else {0}),y=(yyy3+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[9]>0) {sub$percQ2.5[9]} else {0}, if (sub$percQ2.5[9]>0) {sub$percQ2.5[9]} else {0}, if (sub$percQ97.5[9]>0) {sub$percQ97.5[9]} else {0}, if (sub$percQ97.5[9]>0) {sub$percQ97.5[9]} else {0}),y=(yyy1+4), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[9]>0) {sub$percQ5[9]} else {0}, if (sub$percQ5[9]>0) {sub$percQ5[9]} else {0}, if (sub$percQ95[9]>0) {sub$percQ95[9]} else {0}, if (sub$percQ95[9]>0) {sub$percQ95[9]} else {0}),y=(yyy2+4), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[9]>0) {sub$percQ10[9]} else {0}, if (sub$percQ10[9]>0) {sub$percQ10[9]} else {0}, if (sub$percQ90[9]>0) {sub$percQ90[9]} else {0}, if (sub$percQ90[9]>0) {sub$percQ90[9]} else {0}),y=(yyy3+4), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[10]<0) {sub$percQ2.5[10]} else {0}, if (sub$percQ2.5[10]<0) {sub$percQ2.5[10]} else {0}, if (sub$percQ97.5[10]<0) {sub$percQ97.5[10]} else {0}, if (sub$percQ97.5[10]<0) {sub$percQ97.5[10]} else {0}),y=(yyy1+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[10]<0) {sub$percQ5[10]} else {0}, if (sub$percQ5[10]<0) {sub$percQ5[10]} else {0}, if (sub$percQ95[10]<0) {sub$percQ95[10]} else {0}, if (sub$percQ95[10]<0) {sub$percQ95[10]} else {0}),y=(yyy2+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[10]<0) {sub$percQ10[10]} else {0}, if (sub$percQ10[10]<0) {sub$percQ10[10]} else {0}, if (sub$percQ90[10]<0) {sub$percQ90[10]} else {0}, if (sub$percQ90[10]<0) {sub$percQ90[10]} else {0}),y=(yyy3+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[10]>0) {sub$percQ2.5[10]} else {0}, if (sub$percQ2.5[10]>0) {sub$percQ2.5[10]} else {0}, if (sub$percQ97.5[10]>0) {sub$percQ97.5[10]} else {0}, if (sub$percQ97.5[10]>0) {sub$percQ97.5[10]} else {0}),y=(yyy1+3), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[10]>0) {sub$percQ5[10]} else {0}, if (sub$percQ5[10]>0) {sub$percQ5[10]} else {0}, if (sub$percQ95[10]>0) {sub$percQ95[10]} else {0}, if (sub$percQ95[10]>0) {sub$percQ95[10]} else {0}),y=(yyy2+3), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[10]>0) {sub$percQ10[10]} else {0}, if (sub$percQ10[10]>0) {sub$percQ10[10]} else {0}, if (sub$percQ90[10]>0) {sub$percQ90[10]} else {0}, if (sub$percQ90[10]>0) {sub$percQ90[10]} else {0}),y=(yyy3+3), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[11]<0) {sub$percQ2.5[11]} else {0}, if (sub$percQ2.5[11]<0) {sub$percQ2.5[11]} else {0}, if (sub$percQ97.5[11]<0) {sub$percQ97.5[11]} else {0}, if (sub$percQ97.5[11]<0) {sub$percQ97.5[11]} else {0}),y=(yyy1+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[11]<0) {sub$percQ5[11]} else {0}, if (sub$percQ5[11]<0) {sub$percQ5[11]} else {0}, if (sub$percQ95[11]<0) {sub$percQ95[11]} else {0}, if (sub$percQ95[11]<0) {sub$percQ95[11]} else {0}),y=(yyy2+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[11]<0) {sub$percQ10[11]} else {0}, if (sub$percQ10[11]<0) {sub$percQ10[11]} else {0}, if (sub$percQ90[11]<0) {sub$percQ90[11]} else {0}, if (sub$percQ90[11]<0) {sub$percQ90[11]} else {0}),y=(yyy3+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[11]>0) {sub$percQ2.5[11]} else {0}, if (sub$percQ2.5[11]>0) {sub$percQ2.5[11]} else {0}, if (sub$percQ97.5[11]>0) {sub$percQ97.5[11]} else {0}, if (sub$percQ97.5[11]>0) {sub$percQ97.5[11]} else {0}),y=(yyy1+2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[11]>0) {sub$percQ5[11]} else {0}, if (sub$percQ5[11]>0) {sub$percQ5[11]} else {0}, if (sub$percQ95[11]>0) {sub$percQ95[11]} else {0}, if (sub$percQ95[11]>0) {sub$percQ95[11]} else {0}),y=(yyy2+2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[11]>0) {sub$percQ10[11]} else {0}, if (sub$percQ10[11]>0) {sub$percQ10[11]} else {0}, if (sub$percQ90[11]>0) {sub$percQ90[11]} else {0}, if (sub$percQ90[11]>0) {sub$percQ90[11]} else {0}),y=(yyy3+2), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[12]<0) {sub$percQ2.5[12]} else {0}, if (sub$percQ2.5[12]<0) {sub$percQ2.5[12]} else {0}, if (sub$percQ97.5[12]<0) {sub$percQ97.5[12]} else {0}, if (sub$percQ97.5[12]<0) {sub$percQ97.5[12]} else {0}),y=(yyy1+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[12]<0) {sub$percQ5[12]} else {0}, if (sub$percQ5[12]<0) {sub$percQ5[12]} else {0}, if (sub$percQ95[12]<0) {sub$percQ95[12]} else {0}, if (sub$percQ95[12]<0) {sub$percQ95[12]} else {0}),y=(yyy2+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[12]<0) {sub$percQ10[12]} else {0}, if (sub$percQ10[12]<0) {sub$percQ10[12]} else {0}, if (sub$percQ90[12]<0) {sub$percQ90[12]} else {0}, if (sub$percQ90[12]<0) {sub$percQ90[12]} else {0}),y=(yyy3+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[12]>0) {sub$percQ2.5[12]} else {0}, if (sub$percQ2.5[12]>0) {sub$percQ2.5[12]} else {0}, if (sub$percQ97.5[12]>0) {sub$percQ97.5[12]} else {0}, if (sub$percQ97.5[12]>0) {sub$percQ97.5[12]} else {0}),y=(yyy1+1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[12]>0) {sub$percQ5[12]} else {0}, if (sub$percQ5[12]>0) {sub$percQ5[12]} else {0}, if (sub$percQ95[12]>0) {sub$percQ95[12]} else {0}, if (sub$percQ95[12]>0) {sub$percQ95[12]} else {0}),y=(yyy2+1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[12]>0) {sub$percQ10[12]} else {0}, if (sub$percQ10[12]>0) {sub$percQ10[12]} else {0}, if (sub$percQ90[12]>0) {sub$percQ90[12]} else {0}, if (sub$percQ90[12]>0) {sub$percQ90[12]} else {0}),y=(yyy3+1), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[13]<0) {sub$percQ2.5[13]} else {0}, if (sub$percQ2.5[13]<0) {sub$percQ2.5[13]} else {0}, if (sub$percQ97.5[13]<0) {sub$percQ97.5[13]} else {0}, if (sub$percQ97.5[13]<0) {sub$percQ97.5[13]} else {0}),y=(yyy1+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[13]<0) {sub$percQ5[13]} else {0}, if (sub$percQ5[13]<0) {sub$percQ5[13]} else {0}, if (sub$percQ95[13]<0) {sub$percQ95[13]} else {0}, if (sub$percQ95[13]<0) {sub$percQ95[13]} else {0}),y=(yyy2+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[13]<0) {sub$percQ10[13]} else {0}, if (sub$percQ10[13]<0) {sub$percQ10[13]} else {0}, if (sub$percQ90[13]<0) {sub$percQ90[13]} else {0}, if (sub$percQ90[13]<0) {sub$percQ90[13]} else {0}),y=(yyy3+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[13]>0) {sub$percQ2.5[13]} else {0}, if (sub$percQ2.5[13]>0) {sub$percQ2.5[13]} else {0}, if (sub$percQ97.5[13]>0) {sub$percQ97.5[13]} else {0}, if (sub$percQ97.5[13]>0) {sub$percQ97.5[13]} else {0}),y=(yyy1+0), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[13]>0) {sub$percQ5[13]} else {0}, if (sub$percQ5[13]>0) {sub$percQ5[13]} else {0}, if (sub$percQ95[13]>0) {sub$percQ95[13]} else {0}, if (sub$percQ95[13]>0) {sub$percQ95[13]} else {0}),y=(yyy2+0), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[13]>0) {sub$percQ10[13]} else {0}, if (sub$percQ10[13]>0) {sub$percQ10[13]} else {0}, if (sub$percQ90[13]>0) {sub$percQ90[13]} else {0}, if (sub$percQ90[13]>0) {sub$percQ90[13]} else {0}),y=(yyy3+0), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[14]<0) {sub$percQ2.5[14]} else {0}, if (sub$percQ2.5[14]<0) {sub$percQ2.5[14]} else {0}, if (sub$percQ97.5[14]<0) {sub$percQ97.5[14]} else {0}, if (sub$percQ97.5[14]<0) {sub$percQ97.5[14]} else {0}),y=(yyy1-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[14]<0) {sub$percQ5[14]} else {0}, if (sub$percQ5[14]<0) {sub$percQ5[14]} else {0}, if (sub$percQ95[14]<0) {sub$percQ95[14]} else {0}, if (sub$percQ95[14]<0) {sub$percQ95[14]} else {0}),y=(yyy2-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[14]<0) {sub$percQ10[14]} else {0}, if (sub$percQ10[14]<0) {sub$percQ10[14]} else {0}, if (sub$percQ90[14]<0) {sub$percQ90[14]} else {0}, if (sub$percQ90[14]<0) {sub$percQ90[14]} else {0}),y=(yyy3-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[14]>0) {sub$percQ2.5[14]} else {0}, if (sub$percQ2.5[14]>0) {sub$percQ2.5[14]} else {0}, if (sub$percQ97.5[14]>0) {sub$percQ97.5[14]} else {0}, if (sub$percQ97.5[14]>0) {sub$percQ97.5[14]} else {0}),y=(yyy1-1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[14]>0) {sub$percQ5[14]} else {0}, if (sub$percQ5[14]>0) {sub$percQ5[14]} else {0}, if (sub$percQ95[14]>0) {sub$percQ95[14]} else {0}, if (sub$percQ95[14]>0) {sub$percQ95[14]} else {0}),y=(yyy2-1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[14]>0) {sub$percQ10[14]} else {0}, if (sub$percQ10[14]>0) {sub$percQ10[14]} else {0}, if (sub$percQ90[14]>0) {sub$percQ90[14]} else {0}, if (sub$percQ90[14]>0) {sub$percQ90[14]} else {0}),y=(yyy3-1), col = "dodgerblue", border = NA)

abline(h=7.5,lty=2, col="grey60")

###########
sub <- subset(ESTs, Model == "unweighted_random")
yy <- c(13:0)
est <-sub$percEstimate
min(sub$percQ0.5)
max(sub$percQ99.5)
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-17.6,20.1), ylim=c(0,13.5),cex=2, bty="n")
#polygon(x=c(-100,-100,0,0),
#        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
segments(x0=0,y0=0,x1=0,y1=13.3,lty=2, lwd=2,col="grey60")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
legend("top", legend=("b, Unweighted Meta-analysis wRand"), bty="n", cex=1.3)
mm <- cbind(yy,est)
mm_neg <- as.data.frame(mm[ which(est < 0),])
mm_pos <- as.data.frame(mm[ which(est > 0),])
#points(yy ~ est, pch="l",cex=2.2,col="white")
points(mm_neg$yy ~ mm_neg$est, pch="l",cex=2,col="firebrick2")
points(mm_pos$yy ~ mm_pos$est, pch="l",cex=2,col="dodgerblue")
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)

polygon(x=c(if (sub$percQ2.5[1]<0) {sub$percQ2.5[1]} else {0}, if (sub$percQ2.5[1]<0) {sub$percQ2.5[1]} else {0}, if (sub$percQ97.5[1]<0) {sub$percQ97.5[1]} else {0}, if (sub$percQ97.5[1]<0) {sub$percQ97.5[1]} else {0}),y=(yyy1+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[1]<0) {sub$percQ5[1]} else {0}, if (sub$percQ5[1]<0) {sub$percQ5[1]} else {0}, if (sub$percQ95[1]<0) {sub$percQ95[1]} else {0}, if (sub$percQ95[1]<0) {sub$percQ95[1]} else {0}),y=(yyy2+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[1]<0) {sub$percQ10[1]} else {0}, if (sub$percQ10[1]<0) {sub$percQ10[1]} else {0}, if (sub$percQ90[1]<0) {sub$percQ90[1]} else {0}, if (sub$percQ90[1]<0) {sub$percQ90[1]} else {0}),y=(yyy3+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[1]>0) {sub$percQ2.5[1]} else {0}, if (sub$percQ2.5[1]>0) {sub$percQ2.5[1]} else {0}, if (sub$percQ97.5[1]>0) {sub$percQ97.5[1]} else {0}, if (sub$percQ97.5[1]>0) {sub$percQ97.5[1]} else {0}),y=(yyy1+12), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[1]>0) {sub$percQ5[1]} else {0}, if (sub$percQ5[1]>0) {sub$percQ5[1]} else {0}, if (sub$percQ95[1]>0) {sub$percQ95[1]} else {0}, if (sub$percQ95[1]>0) {sub$percQ95[1]} else {0}),y=(yyy2+12), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[1]>0) {sub$percQ10[1]} else {0}, if (sub$percQ10[1]>0) {sub$percQ10[1]} else {0}, if (sub$percQ90[1]>0) {sub$percQ90[1]} else {0}, if (sub$percQ90[1]>0) {sub$percQ90[1]} else {0}),y=(yyy3+12), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[2]<0) {sub$percQ2.5[2]} else {0}, if (sub$percQ2.5[2]<0) {sub$percQ2.5[2]} else {0}, if (sub$percQ97.5[2]<0) {sub$percQ97.5[2]} else {0}, if (sub$percQ97.5[2]<0) {sub$percQ97.5[2]} else {0}),y=(yyy1+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[2]<0) {sub$percQ5[2]} else {0}, if (sub$percQ5[2]<0) {sub$percQ5[2]} else {0}, if (sub$percQ95[2]<0) {sub$percQ95[2]} else {0}, if (sub$percQ95[2]<0) {sub$percQ95[2]} else {0}),y=(yyy2+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[2]<0) {sub$percQ10[2]} else {0}, if (sub$percQ10[2]<0) {sub$percQ10[2]} else {0}, if (sub$percQ90[2]<0) {sub$percQ90[2]} else {0}, if (sub$percQ90[2]<0) {sub$percQ90[2]} else {0}),y=(yyy3+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[2]>0) {sub$percQ2.5[2]} else {0}, if (sub$percQ2.5[2]>0) {sub$percQ2.5[2]} else {0}, if (sub$percQ97.5[2]>0) {sub$percQ97.5[2]} else {0}, if (sub$percQ97.5[2]>0) {sub$percQ97.5[2]} else {0}),y=(yyy1+11), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[2]>0) {sub$percQ5[2]} else {0}, if (sub$percQ5[2]>0) {sub$percQ5[2]} else {0}, if (sub$percQ95[2]>0) {sub$percQ95[2]} else {0}, if (sub$percQ95[2]>0) {sub$percQ95[2]} else {0}),y=(yyy2+11), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[2]>0) {sub$percQ10[2]} else {0}, if (sub$percQ10[2]>0) {sub$percQ10[2]} else {0}, if (sub$percQ90[2]>0) {sub$percQ90[2]} else {0}, if (sub$percQ90[2]>0) {sub$percQ90[2]} else {0}),y=(yyy3+11), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[3]<0) {sub$percQ2.5[3]} else {0}, if (sub$percQ2.5[3]<0) {sub$percQ2.5[3]} else {0}, if (sub$percQ97.5[3]<0) {sub$percQ97.5[3]} else {0}, if (sub$percQ97.5[3]<0) {sub$percQ97.5[3]} else {0}),y=(yyy1+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[3]<0) {sub$percQ5[3]} else {0}, if (sub$percQ5[3]<0) {sub$percQ5[3]} else {0}, if (sub$percQ95[3]<0) {sub$percQ95[3]} else {0}, if (sub$percQ95[3]<0) {sub$percQ95[3]} else {0}),y=(yyy2+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[3]<0) {sub$percQ10[3]} else {0}, if (sub$percQ10[3]<0) {sub$percQ10[3]} else {0}, if (sub$percQ90[3]<0) {sub$percQ90[3]} else {0}, if (sub$percQ90[3]<0) {sub$percQ90[3]} else {0}),y=(yyy3+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[3]>0) {sub$percQ2.5[3]} else {0}, if (sub$percQ2.5[3]>0) {sub$percQ2.5[3]} else {0}, if (sub$percQ97.5[3]>0) {sub$percQ97.5[3]} else {0}, if (sub$percQ97.5[3]>0) {sub$percQ97.5[3]} else {0}),y=(yyy1+10), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[3]>0) {sub$percQ5[3]} else {0}, if (sub$percQ5[3]>0) {sub$percQ5[3]} else {0}, if (sub$percQ95[3]>0) {sub$percQ95[3]} else {0}, if (sub$percQ95[3]>0) {sub$percQ95[3]} else {0}),y=(yyy2+10), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[3]>0) {sub$percQ10[3]} else {0}, if (sub$percQ10[3]>0) {sub$percQ10[3]} else {0}, if (sub$percQ90[3]>0) {sub$percQ90[3]} else {0}, if (sub$percQ90[3]>0) {sub$percQ90[3]} else {0}),y=(yyy3+10), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[4]<0) {sub$percQ2.5[4]} else {0}, if (sub$percQ2.5[4]<0) {sub$percQ2.5[4]} else {0}, if (sub$percQ97.5[4]<0) {sub$percQ97.5[4]} else {0}, if (sub$percQ97.5[4]<0) {sub$percQ97.5[4]} else {0}),y=(yyy1+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[4]<0) {sub$percQ5[4]} else {0}, if (sub$percQ5[4]<0) {sub$percQ5[4]} else {0}, if (sub$percQ95[4]<0) {sub$percQ95[4]} else {0}, if (sub$percQ95[4]<0) {sub$percQ95[4]} else {0}),y=(yyy2+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[4]<0) {sub$percQ10[4]} else {0}, if (sub$percQ10[4]<0) {sub$percQ10[4]} else {0}, if (sub$percQ90[4]<0) {sub$percQ90[4]} else {0}, if (sub$percQ90[4]<0) {sub$percQ90[4]} else {0}),y=(yyy3+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[4]>0) {sub$percQ2.5[4]} else {0}, if (sub$percQ2.5[4]>0) {sub$percQ2.5[4]} else {0}, if (sub$percQ97.5[4]>0) {sub$percQ97.5[4]} else {0}, if (sub$percQ97.5[4]>0) {sub$percQ97.5[4]} else {0}),y=(yyy1+9), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[4]>0) {sub$percQ5[4]} else {0}, if (sub$percQ5[4]>0) {sub$percQ5[4]} else {0}, if (sub$percQ95[4]>0) {sub$percQ95[4]} else {0}, if (sub$percQ95[4]>0) {sub$percQ95[4]} else {0}),y=(yyy2+9), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[4]>0) {sub$percQ10[4]} else {0}, if (sub$percQ10[4]>0) {sub$percQ10[4]} else {0}, if (sub$percQ90[4]>0) {sub$percQ90[4]} else {0}, if (sub$percQ90[4]>0) {sub$percQ90[4]} else {0}),y=(yyy3+9), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[5]<0) {sub$percQ2.5[5]} else {0}, if (sub$percQ2.5[5]<0) {sub$percQ2.5[5]} else {0}, if (sub$percQ97.5[5]<0) {sub$percQ97.5[5]} else {0}, if (sub$percQ97.5[5]<0) {sub$percQ97.5[5]} else {0}),y=(yyy1+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[5]<0) {sub$percQ5[5]} else {0}, if (sub$percQ5[5]<0) {sub$percQ5[5]} else {0}, if (sub$percQ95[5]<0) {sub$percQ95[5]} else {0}, if (sub$percQ95[5]<0) {sub$percQ95[5]} else {0}),y=(yyy2+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[5]<0) {sub$percQ10[5]} else {0}, if (sub$percQ10[5]<0) {sub$percQ10[5]} else {0}, if (sub$percQ90[5]<0) {sub$percQ90[5]} else {0}, if (sub$percQ90[5]<0) {sub$percQ90[5]} else {0}),y=(yyy3+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[5]>0) {sub$percQ2.5[5]} else {0}, if (sub$percQ2.5[5]>0) {sub$percQ2.5[5]} else {0}, if (sub$percQ97.5[5]>0) {sub$percQ97.5[5]} else {0}, if (sub$percQ97.5[5]>0) {sub$percQ97.5[5]} else {0}),y=(yyy1+8), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[5]>0) {sub$percQ5[5]} else {0}, if (sub$percQ5[5]>0) {sub$percQ5[5]} else {0}, if (sub$percQ95[5]>0) {sub$percQ95[5]} else {0}, if (sub$percQ95[5]>0) {sub$percQ95[5]} else {0}),y=(yyy2+8), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[5]>0) {sub$percQ10[5]} else {0}, if (sub$percQ10[5]>0) {sub$percQ10[5]} else {0}, if (sub$percQ90[5]>0) {sub$percQ90[5]} else {0}, if (sub$percQ90[5]>0) {sub$percQ90[5]} else {0}),y=(yyy3+8), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[6]<0) {sub$percQ2.5[6]} else {0}, if (sub$percQ2.5[6]<0) {sub$percQ2.5[6]} else {0}, if (sub$percQ97.5[6]<0) {sub$percQ97.5[6]} else {0}, if (sub$percQ97.5[6]<0) {sub$percQ97.5[6]} else {0}),y=(yyy1+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[6]<0) {sub$percQ5[6]} else {0}, if (sub$percQ5[6]<0) {sub$percQ5[6]} else {0}, if (sub$percQ95[6]<0) {sub$percQ95[6]} else {0}, if (sub$percQ95[6]<0) {sub$percQ95[6]} else {0}),y=(yyy2+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[6]<0) {sub$percQ10[6]} else {0}, if (sub$percQ10[6]<0) {sub$percQ10[6]} else {0}, if (sub$percQ90[6]<0) {sub$percQ90[6]} else {0}, if (sub$percQ90[6]<0) {sub$percQ90[6]} else {0}),y=(yyy3+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[6]>0) {sub$percQ2.5[6]} else {0}, if (sub$percQ2.5[6]>0) {sub$percQ2.5[6]} else {0}, if (sub$percQ97.5[6]>0) {sub$percQ97.5[6]} else {0}, if (sub$percQ97.5[6]>0) {sub$percQ97.5[6]} else {0}),y=(yyy1+7), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[6]>0) {sub$percQ5[6]} else {0}, if (sub$percQ5[6]>0) {sub$percQ5[6]} else {0}, if (sub$percQ95[6]>0) {sub$percQ95[6]} else {0}, if (sub$percQ95[6]>0) {sub$percQ95[6]} else {0}),y=(yyy2+7), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[6]>0) {sub$percQ10[6]} else {0}, if (sub$percQ10[6]>0) {sub$percQ10[6]} else {0}, if (sub$percQ90[6]>0) {sub$percQ90[6]} else {0}, if (sub$percQ90[6]>0) {sub$percQ90[6]} else {0}),y=(yyy3+7), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[7]<0) {sub$percQ2.5[7]} else {0}, if (sub$percQ2.5[7]<0) {sub$percQ2.5[7]} else {0}, if (sub$percQ97.5[7]<0) {sub$percQ97.5[7]} else {0}, if (sub$percQ97.5[7]<0) {sub$percQ97.5[7]} else {0}),y=(yyy1+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[7]<0) {sub$percQ5[7]} else {0}, if (sub$percQ5[7]<0) {sub$percQ5[7]} else {0}, if (sub$percQ95[7]<0) {sub$percQ95[7]} else {0}, if (sub$percQ95[7]<0) {sub$percQ95[7]} else {0}),y=(yyy2+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[7]<0) {sub$percQ10[7]} else {0}, if (sub$percQ10[7]<0) {sub$percQ10[7]} else {0}, if (sub$percQ90[7]<0) {sub$percQ90[7]} else {0}, if (sub$percQ90[7]<0) {sub$percQ90[7]} else {0}),y=(yyy3+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[7]>0) {sub$percQ2.5[7]} else {0}, if (sub$percQ2.5[7]>0) {sub$percQ2.5[7]} else {0}, if (sub$percQ97.5[7]>0) {sub$percQ97.5[7]} else {0}, if (sub$percQ97.5[7]>0) {sub$percQ97.5[7]} else {0}),y=(yyy1+6), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[7]>0) {sub$percQ5[7]} else {0}, if (sub$percQ5[7]>0) {sub$percQ5[7]} else {0}, if (sub$percQ95[7]>0) {sub$percQ95[7]} else {0}, if (sub$percQ95[7]>0) {sub$percQ95[7]} else {0}),y=(yyy2+6), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[7]>0) {sub$percQ10[7]} else {0}, if (sub$percQ10[7]>0) {sub$percQ10[7]} else {0}, if (sub$percQ90[7]>0) {sub$percQ90[7]} else {0}, if (sub$percQ90[7]>0) {sub$percQ90[7]} else {0}),y=(yyy3+6), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[8]<0) {sub$percQ2.5[8]} else {0}, if (sub$percQ2.5[8]<0) {sub$percQ2.5[8]} else {0}, if (sub$percQ97.5[8]<0) {sub$percQ97.5[8]} else {0}, if (sub$percQ97.5[8]<0) {sub$percQ97.5[8]} else {0}),y=(yyy1+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[8]<0) {sub$percQ5[8]} else {0}, if (sub$percQ5[8]<0) {sub$percQ5[8]} else {0}, if (sub$percQ95[8]<0) {sub$percQ95[8]} else {0}, if (sub$percQ95[8]<0) {sub$percQ95[8]} else {0}),y=(yyy2+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[8]<0) {sub$percQ10[8]} else {0}, if (sub$percQ10[8]<0) {sub$percQ10[8]} else {0}, if (sub$percQ90[8]<0) {sub$percQ90[8]} else {0}, if (sub$percQ90[8]<0) {sub$percQ90[8]} else {0}),y=(yyy3+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[8]>0) {sub$percQ2.5[8]} else {0}, if (sub$percQ2.5[8]>0) {sub$percQ2.5[8]} else {0}, if (sub$percQ97.5[8]>0) {sub$percQ97.5[8]} else {0}, if (sub$percQ97.5[8]>0) {sub$percQ97.5[8]} else {0}),y=(yyy1+5), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[8]>0) {sub$percQ5[8]} else {0}, if (sub$percQ5[8]>0) {sub$percQ5[8]} else {0}, if (sub$percQ95[8]>0) {sub$percQ95[8]} else {0}, if (sub$percQ95[8]>0) {sub$percQ95[8]} else {0}),y=(yyy2+5), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[8]>0) {sub$percQ10[8]} else {0}, if (sub$percQ10[8]>0) {sub$percQ10[8]} else {0}, if (sub$percQ90[8]>0) {sub$percQ90[8]} else {0}, if (sub$percQ90[8]>0) {sub$percQ90[8]} else {0}),y=(yyy3+5), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[9]<0) {sub$percQ2.5[9]} else {0}, if (sub$percQ2.5[9]<0) {sub$percQ2.5[9]} else {0}, if (sub$percQ97.5[9]<0) {sub$percQ97.5[9]} else {0}, if (sub$percQ97.5[9]<0) {sub$percQ97.5[9]} else {0}),y=(yyy1+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[9]<0) {sub$percQ5[9]} else {0}, if (sub$percQ5[9]<0) {sub$percQ5[9]} else {0}, if (sub$percQ95[9]<0) {sub$percQ95[9]} else {0}, if (sub$percQ95[9]<0) {sub$percQ95[9]} else {0}),y=(yyy2+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[9]<0) {sub$percQ10[9]} else {0}, if (sub$percQ10[9]<0) {sub$percQ10[9]} else {0}, if (sub$percQ90[9]<0) {sub$percQ90[9]} else {0}, if (sub$percQ90[9]<0) {sub$percQ90[9]} else {0}),y=(yyy3+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[9]>0) {sub$percQ2.5[9]} else {0}, if (sub$percQ2.5[9]>0) {sub$percQ2.5[9]} else {0}, if (sub$percQ97.5[9]>0) {sub$percQ97.5[9]} else {0}, if (sub$percQ97.5[9]>0) {sub$percQ97.5[9]} else {0}),y=(yyy1+4), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[9]>0) {sub$percQ5[9]} else {0}, if (sub$percQ5[9]>0) {sub$percQ5[9]} else {0}, if (sub$percQ95[9]>0) {sub$percQ95[9]} else {0}, if (sub$percQ95[9]>0) {sub$percQ95[9]} else {0}),y=(yyy2+4), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[9]>0) {sub$percQ10[9]} else {0}, if (sub$percQ10[9]>0) {sub$percQ10[9]} else {0}, if (sub$percQ90[9]>0) {sub$percQ90[9]} else {0}, if (sub$percQ90[9]>0) {sub$percQ90[9]} else {0}),y=(yyy3+4), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[10]<0) {sub$percQ2.5[10]} else {0}, if (sub$percQ2.5[10]<0) {sub$percQ2.5[10]} else {0}, if (sub$percQ97.5[10]<0) {sub$percQ97.5[10]} else {0}, if (sub$percQ97.5[10]<0) {sub$percQ97.5[10]} else {0}),y=(yyy1+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[10]<0) {sub$percQ5[10]} else {0}, if (sub$percQ5[10]<0) {sub$percQ5[10]} else {0}, if (sub$percQ95[10]<0) {sub$percQ95[10]} else {0}, if (sub$percQ95[10]<0) {sub$percQ95[10]} else {0}),y=(yyy2+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[10]<0) {sub$percQ10[10]} else {0}, if (sub$percQ10[10]<0) {sub$percQ10[10]} else {0}, if (sub$percQ90[10]<0) {sub$percQ90[10]} else {0}, if (sub$percQ90[10]<0) {sub$percQ90[10]} else {0}),y=(yyy3+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[10]>0) {sub$percQ2.5[10]} else {0}, if (sub$percQ2.5[10]>0) {sub$percQ2.5[10]} else {0}, if (sub$percQ97.5[10]>0) {sub$percQ97.5[10]} else {0}, if (sub$percQ97.5[10]>0) {sub$percQ97.5[10]} else {0}),y=(yyy1+3), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[10]>0) {sub$percQ5[10]} else {0}, if (sub$percQ5[10]>0) {sub$percQ5[10]} else {0}, if (sub$percQ95[10]>0) {sub$percQ95[10]} else {0}, if (sub$percQ95[10]>0) {sub$percQ95[10]} else {0}),y=(yyy2+3), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[10]>0) {sub$percQ10[10]} else {0}, if (sub$percQ10[10]>0) {sub$percQ10[10]} else {0}, if (sub$percQ90[10]>0) {sub$percQ90[10]} else {0}, if (sub$percQ90[10]>0) {sub$percQ90[10]} else {0}),y=(yyy3+3), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[11]<0) {sub$percQ2.5[11]} else {0}, if (sub$percQ2.5[11]<0) {sub$percQ2.5[11]} else {0}, if (sub$percQ97.5[11]<0) {sub$percQ97.5[11]} else {0}, if (sub$percQ97.5[11]<0) {sub$percQ97.5[11]} else {0}),y=(yyy1+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[11]<0) {sub$percQ5[11]} else {0}, if (sub$percQ5[11]<0) {sub$percQ5[11]} else {0}, if (sub$percQ95[11]<0) {sub$percQ95[11]} else {0}, if (sub$percQ95[11]<0) {sub$percQ95[11]} else {0}),y=(yyy2+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[11]<0) {sub$percQ10[11]} else {0}, if (sub$percQ10[11]<0) {sub$percQ10[11]} else {0}, if (sub$percQ90[11]<0) {sub$percQ90[11]} else {0}, if (sub$percQ90[11]<0) {sub$percQ90[11]} else {0}),y=(yyy3+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[11]>0) {sub$percQ2.5[11]} else {0}, if (sub$percQ2.5[11]>0) {sub$percQ2.5[11]} else {0}, if (sub$percQ97.5[11]>0) {sub$percQ97.5[11]} else {0}, if (sub$percQ97.5[11]>0) {sub$percQ97.5[11]} else {0}),y=(yyy1+2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[11]>0) {sub$percQ5[11]} else {0}, if (sub$percQ5[11]>0) {sub$percQ5[11]} else {0}, if (sub$percQ95[11]>0) {sub$percQ95[11]} else {0}, if (sub$percQ95[11]>0) {sub$percQ95[11]} else {0}),y=(yyy2+2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[11]>0) {sub$percQ10[11]} else {0}, if (sub$percQ10[11]>0) {sub$percQ10[11]} else {0}, if (sub$percQ90[11]>0) {sub$percQ90[11]} else {0}, if (sub$percQ90[11]>0) {sub$percQ90[11]} else {0}),y=(yyy3+2), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[12]<0) {sub$percQ2.5[12]} else {0}, if (sub$percQ2.5[12]<0) {sub$percQ2.5[12]} else {0}, if (sub$percQ97.5[12]<0) {sub$percQ97.5[12]} else {0}, if (sub$percQ97.5[12]<0) {sub$percQ97.5[12]} else {0}),y=(yyy1+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[12]<0) {sub$percQ5[12]} else {0}, if (sub$percQ5[12]<0) {sub$percQ5[12]} else {0}, if (sub$percQ95[12]<0) {sub$percQ95[12]} else {0}, if (sub$percQ95[12]<0) {sub$percQ95[12]} else {0}),y=(yyy2+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[12]<0) {sub$percQ10[12]} else {0}, if (sub$percQ10[12]<0) {sub$percQ10[12]} else {0}, if (sub$percQ90[12]<0) {sub$percQ90[12]} else {0}, if (sub$percQ90[12]<0) {sub$percQ90[12]} else {0}),y=(yyy3+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[12]>0) {sub$percQ2.5[12]} else {0}, if (sub$percQ2.5[12]>0) {sub$percQ2.5[12]} else {0}, if (sub$percQ97.5[12]>0) {sub$percQ97.5[12]} else {0}, if (sub$percQ97.5[12]>0) {sub$percQ97.5[12]} else {0}),y=(yyy1+1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[12]>0) {sub$percQ5[12]} else {0}, if (sub$percQ5[12]>0) {sub$percQ5[12]} else {0}, if (sub$percQ95[12]>0) {sub$percQ95[12]} else {0}, if (sub$percQ95[12]>0) {sub$percQ95[12]} else {0}),y=(yyy2+1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[12]>0) {sub$percQ10[12]} else {0}, if (sub$percQ10[12]>0) {sub$percQ10[12]} else {0}, if (sub$percQ90[12]>0) {sub$percQ90[12]} else {0}, if (sub$percQ90[12]>0) {sub$percQ90[12]} else {0}),y=(yyy3+1), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[13]<0) {sub$percQ2.5[13]} else {0}, if (sub$percQ2.5[13]<0) {sub$percQ2.5[13]} else {0}, if (sub$percQ97.5[13]<0) {sub$percQ97.5[13]} else {0}, if (sub$percQ97.5[13]<0) {sub$percQ97.5[13]} else {0}),y=(yyy1+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[13]<0) {sub$percQ5[13]} else {0}, if (sub$percQ5[13]<0) {sub$percQ5[13]} else {0}, if (sub$percQ95[13]<0) {sub$percQ95[13]} else {0}, if (sub$percQ95[13]<0) {sub$percQ95[13]} else {0}),y=(yyy2+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[13]<0) {sub$percQ10[13]} else {0}, if (sub$percQ10[13]<0) {sub$percQ10[13]} else {0}, if (sub$percQ90[13]<0) {sub$percQ90[13]} else {0}, if (sub$percQ90[13]<0) {sub$percQ90[13]} else {0}),y=(yyy3+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[13]>0) {sub$percQ2.5[13]} else {0}, if (sub$percQ2.5[13]>0) {sub$percQ2.5[13]} else {0}, if (sub$percQ97.5[13]>0) {sub$percQ97.5[13]} else {0}, if (sub$percQ97.5[13]>0) {sub$percQ97.5[13]} else {0}),y=(yyy1+0), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[13]>0) {sub$percQ5[13]} else {0}, if (sub$percQ5[13]>0) {sub$percQ5[13]} else {0}, if (sub$percQ95[13]>0) {sub$percQ95[13]} else {0}, if (sub$percQ95[13]>0) {sub$percQ95[13]} else {0}),y=(yyy2+0), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[13]>0) {sub$percQ10[13]} else {0}, if (sub$percQ10[13]>0) {sub$percQ10[13]} else {0}, if (sub$percQ90[13]>0) {sub$percQ90[13]} else {0}, if (sub$percQ90[13]>0) {sub$percQ90[13]} else {0}),y=(yyy3+0), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[14]<0) {sub$percQ2.5[14]} else {0}, if (sub$percQ2.5[14]<0) {sub$percQ2.5[14]} else {0}, if (sub$percQ97.5[14]<0) {sub$percQ97.5[14]} else {0}, if (sub$percQ97.5[14]<0) {sub$percQ97.5[14]} else {0}),y=(yyy1-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[14]<0) {sub$percQ5[14]} else {0}, if (sub$percQ5[14]<0) {sub$percQ5[14]} else {0}, if (sub$percQ95[14]<0) {sub$percQ95[14]} else {0}, if (sub$percQ95[14]<0) {sub$percQ95[14]} else {0}),y=(yyy2-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[14]<0) {sub$percQ10[14]} else {0}, if (sub$percQ10[14]<0) {sub$percQ10[14]} else {0}, if (sub$percQ90[14]<0) {sub$percQ90[14]} else {0}, if (sub$percQ90[14]<0) {sub$percQ90[14]} else {0}),y=(yyy3-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[14]>0) {sub$percQ2.5[14]} else {0}, if (sub$percQ2.5[14]>0) {sub$percQ2.5[14]} else {0}, if (sub$percQ97.5[14]>0) {sub$percQ97.5[14]} else {0}, if (sub$percQ97.5[14]>0) {sub$percQ97.5[14]} else {0}),y=(yyy1-1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[14]>0) {sub$percQ5[14]} else {0}, if (sub$percQ5[14]>0) {sub$percQ5[14]} else {0}, if (sub$percQ95[14]>0) {sub$percQ95[14]} else {0}, if (sub$percQ95[14]>0) {sub$percQ95[14]} else {0}),y=(yyy2-1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[14]>0) {sub$percQ10[14]} else {0}, if (sub$percQ10[14]>0) {sub$percQ10[14]} else {0}, if (sub$percQ90[14]>0) {sub$percQ90[14]} else {0}, if (sub$percQ90[14]>0) {sub$percQ90[14]} else {0}),y=(yyy3-1), col = "dodgerblue", border = NA)

abline(h=7.5,lty=2, col="grey60")

############
sub <- subset(ESTs, Model == "Weighted_noRandom")
yy <- c(13:0)
est <-sub$percEstimate
min(sub$percQ0.5)
max(sub$percQ99.5)
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-10.2,7.3), ylim=c(0,13.5),cex=2, bty="n")
#polygon(x=c(-100,-100,0,0),
#        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
segments(x0=0,y0=0,x1=0,y1=13.3,lty=2, lwd=2,col="grey60")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
legend("top", legend=("c, Weighted Meta-analysis noRand"), bty="n", cex=1.3)
mm <- cbind(yy,est)
mm_neg <- as.data.frame(mm[ which(est < 0),])
mm_pos <- as.data.frame(mm[ which(est > 0),])
#points(yy ~ est, pch="l",cex=2.2,col="white")
points(mm_neg$yy ~ mm_neg$est, pch="l",cex=2,col="firebrick2")
points(mm_pos$yy ~ mm_pos$est, pch="l",cex=2,col="dodgerblue")
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)

polygon(x=c(if (sub$percQ2.5[1]<0) {sub$percQ2.5[1]} else {0}, if (sub$percQ2.5[1]<0) {sub$percQ2.5[1]} else {0}, if (sub$percQ97.5[1]<0) {sub$percQ97.5[1]} else {0}, if (sub$percQ97.5[1]<0) {sub$percQ97.5[1]} else {0}),y=(yyy1+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[1]<0) {sub$percQ5[1]} else {0}, if (sub$percQ5[1]<0) {sub$percQ5[1]} else {0}, if (sub$percQ95[1]<0) {sub$percQ95[1]} else {0}, if (sub$percQ95[1]<0) {sub$percQ95[1]} else {0}),y=(yyy2+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[1]<0) {sub$percQ10[1]} else {0}, if (sub$percQ10[1]<0) {sub$percQ10[1]} else {0}, if (sub$percQ90[1]<0) {sub$percQ90[1]} else {0}, if (sub$percQ90[1]<0) {sub$percQ90[1]} else {0}),y=(yyy3+12), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[1]>0) {sub$percQ2.5[1]} else {0}, if (sub$percQ2.5[1]>0) {sub$percQ2.5[1]} else {0}, if (sub$percQ97.5[1]>0) {sub$percQ97.5[1]} else {0}, if (sub$percQ97.5[1]>0) {sub$percQ97.5[1]} else {0}),y=(yyy1+12), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[1]>0) {sub$percQ5[1]} else {0}, if (sub$percQ5[1]>0) {sub$percQ5[1]} else {0}, if (sub$percQ95[1]>0) {sub$percQ95[1]} else {0}, if (sub$percQ95[1]>0) {sub$percQ95[1]} else {0}),y=(yyy2+12), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[1]>0) {sub$percQ10[1]} else {0}, if (sub$percQ10[1]>0) {sub$percQ10[1]} else {0}, if (sub$percQ90[1]>0) {sub$percQ90[1]} else {0}, if (sub$percQ90[1]>0) {sub$percQ90[1]} else {0}),y=(yyy3+12), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[2]<0) {sub$percQ2.5[2]} else {0}, if (sub$percQ2.5[2]<0) {sub$percQ2.5[2]} else {0}, if (sub$percQ97.5[2]<0) {sub$percQ97.5[2]} else {0}, if (sub$percQ97.5[2]<0) {sub$percQ97.5[2]} else {0}),y=(yyy1+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[2]<0) {sub$percQ5[2]} else {0}, if (sub$percQ5[2]<0) {sub$percQ5[2]} else {0}, if (sub$percQ95[2]<0) {sub$percQ95[2]} else {0}, if (sub$percQ95[2]<0) {sub$percQ95[2]} else {0}),y=(yyy2+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[2]<0) {sub$percQ10[2]} else {0}, if (sub$percQ10[2]<0) {sub$percQ10[2]} else {0}, if (sub$percQ90[2]<0) {sub$percQ90[2]} else {0}, if (sub$percQ90[2]<0) {sub$percQ90[2]} else {0}),y=(yyy3+11), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[2]>0) {sub$percQ2.5[2]} else {0}, if (sub$percQ2.5[2]>0) {sub$percQ2.5[2]} else {0}, if (sub$percQ97.5[2]>0) {sub$percQ97.5[2]} else {0}, if (sub$percQ97.5[2]>0) {sub$percQ97.5[2]} else {0}),y=(yyy1+11), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[2]>0) {sub$percQ5[2]} else {0}, if (sub$percQ5[2]>0) {sub$percQ5[2]} else {0}, if (sub$percQ95[2]>0) {sub$percQ95[2]} else {0}, if (sub$percQ95[2]>0) {sub$percQ95[2]} else {0}),y=(yyy2+11), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[2]>0) {sub$percQ10[2]} else {0}, if (sub$percQ10[2]>0) {sub$percQ10[2]} else {0}, if (sub$percQ90[2]>0) {sub$percQ90[2]} else {0}, if (sub$percQ90[2]>0) {sub$percQ90[2]} else {0}),y=(yyy3+11), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[3]<0) {sub$percQ2.5[3]} else {0}, if (sub$percQ2.5[3]<0) {sub$percQ2.5[3]} else {0}, if (sub$percQ97.5[3]<0) {sub$percQ97.5[3]} else {0}, if (sub$percQ97.5[3]<0) {sub$percQ97.5[3]} else {0}),y=(yyy1+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[3]<0) {sub$percQ5[3]} else {0}, if (sub$percQ5[3]<0) {sub$percQ5[3]} else {0}, if (sub$percQ95[3]<0) {sub$percQ95[3]} else {0}, if (sub$percQ95[3]<0) {sub$percQ95[3]} else {0}),y=(yyy2+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[3]<0) {sub$percQ10[3]} else {0}, if (sub$percQ10[3]<0) {sub$percQ10[3]} else {0}, if (sub$percQ90[3]<0) {sub$percQ90[3]} else {0}, if (sub$percQ90[3]<0) {sub$percQ90[3]} else {0}),y=(yyy3+10), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[3]>0) {sub$percQ2.5[3]} else {0}, if (sub$percQ2.5[3]>0) {sub$percQ2.5[3]} else {0}, if (sub$percQ97.5[3]>0) {sub$percQ97.5[3]} else {0}, if (sub$percQ97.5[3]>0) {sub$percQ97.5[3]} else {0}),y=(yyy1+10), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[3]>0) {sub$percQ5[3]} else {0}, if (sub$percQ5[3]>0) {sub$percQ5[3]} else {0}, if (sub$percQ95[3]>0) {sub$percQ95[3]} else {0}, if (sub$percQ95[3]>0) {sub$percQ95[3]} else {0}),y=(yyy2+10), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[3]>0) {sub$percQ10[3]} else {0}, if (sub$percQ10[3]>0) {sub$percQ10[3]} else {0}, if (sub$percQ90[3]>0) {sub$percQ90[3]} else {0}, if (sub$percQ90[3]>0) {sub$percQ90[3]} else {0}),y=(yyy3+10), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[4]<0) {sub$percQ2.5[4]} else {0}, if (sub$percQ2.5[4]<0) {sub$percQ2.5[4]} else {0}, if (sub$percQ97.5[4]<0) {sub$percQ97.5[4]} else {0}, if (sub$percQ97.5[4]<0) {sub$percQ97.5[4]} else {0}),y=(yyy1+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[4]<0) {sub$percQ5[4]} else {0}, if (sub$percQ5[4]<0) {sub$percQ5[4]} else {0}, if (sub$percQ95[4]<0) {sub$percQ95[4]} else {0}, if (sub$percQ95[4]<0) {sub$percQ95[4]} else {0}),y=(yyy2+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[4]<0) {sub$percQ10[4]} else {0}, if (sub$percQ10[4]<0) {sub$percQ10[4]} else {0}, if (sub$percQ90[4]<0) {sub$percQ90[4]} else {0}, if (sub$percQ90[4]<0) {sub$percQ90[4]} else {0}),y=(yyy3+9), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[4]>0) {sub$percQ2.5[4]} else {0}, if (sub$percQ2.5[4]>0) {sub$percQ2.5[4]} else {0}, if (sub$percQ97.5[4]>0) {sub$percQ97.5[4]} else {0}, if (sub$percQ97.5[4]>0) {sub$percQ97.5[4]} else {0}),y=(yyy1+9), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[4]>0) {sub$percQ5[4]} else {0}, if (sub$percQ5[4]>0) {sub$percQ5[4]} else {0}, if (sub$percQ95[4]>0) {sub$percQ95[4]} else {0}, if (sub$percQ95[4]>0) {sub$percQ95[4]} else {0}),y=(yyy2+9), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[4]>0) {sub$percQ10[4]} else {0}, if (sub$percQ10[4]>0) {sub$percQ10[4]} else {0}, if (sub$percQ90[4]>0) {sub$percQ90[4]} else {0}, if (sub$percQ90[4]>0) {sub$percQ90[4]} else {0}),y=(yyy3+9), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[5]<0) {sub$percQ2.5[5]} else {0}, if (sub$percQ2.5[5]<0) {sub$percQ2.5[5]} else {0}, if (sub$percQ97.5[5]<0) {sub$percQ97.5[5]} else {0}, if (sub$percQ97.5[5]<0) {sub$percQ97.5[5]} else {0}),y=(yyy1+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[5]<0) {sub$percQ5[5]} else {0}, if (sub$percQ5[5]<0) {sub$percQ5[5]} else {0}, if (sub$percQ95[5]<0) {sub$percQ95[5]} else {0}, if (sub$percQ95[5]<0) {sub$percQ95[5]} else {0}),y=(yyy2+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[5]<0) {sub$percQ10[5]} else {0}, if (sub$percQ10[5]<0) {sub$percQ10[5]} else {0}, if (sub$percQ90[5]<0) {sub$percQ90[5]} else {0}, if (sub$percQ90[5]<0) {sub$percQ90[5]} else {0}),y=(yyy3+8), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[5]>0) {sub$percQ2.5[5]} else {0}, if (sub$percQ2.5[5]>0) {sub$percQ2.5[5]} else {0}, if (sub$percQ97.5[5]>0) {sub$percQ97.5[5]} else {0}, if (sub$percQ97.5[5]>0) {sub$percQ97.5[5]} else {0}),y=(yyy1+8), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[5]>0) {sub$percQ5[5]} else {0}, if (sub$percQ5[5]>0) {sub$percQ5[5]} else {0}, if (sub$percQ95[5]>0) {sub$percQ95[5]} else {0}, if (sub$percQ95[5]>0) {sub$percQ95[5]} else {0}),y=(yyy2+8), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[5]>0) {sub$percQ10[5]} else {0}, if (sub$percQ10[5]>0) {sub$percQ10[5]} else {0}, if (sub$percQ90[5]>0) {sub$percQ90[5]} else {0}, if (sub$percQ90[5]>0) {sub$percQ90[5]} else {0}),y=(yyy3+8), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[6]<0) {sub$percQ2.5[6]} else {0}, if (sub$percQ2.5[6]<0) {sub$percQ2.5[6]} else {0}, if (sub$percQ97.5[6]<0) {sub$percQ97.5[6]} else {0}, if (sub$percQ97.5[6]<0) {sub$percQ97.5[6]} else {0}),y=(yyy1+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[6]<0) {sub$percQ5[6]} else {0}, if (sub$percQ5[6]<0) {sub$percQ5[6]} else {0}, if (sub$percQ95[6]<0) {sub$percQ95[6]} else {0}, if (sub$percQ95[6]<0) {sub$percQ95[6]} else {0}),y=(yyy2+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[6]<0) {sub$percQ10[6]} else {0}, if (sub$percQ10[6]<0) {sub$percQ10[6]} else {0}, if (sub$percQ90[6]<0) {sub$percQ90[6]} else {0}, if (sub$percQ90[6]<0) {sub$percQ90[6]} else {0}),y=(yyy3+7), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[6]>0) {sub$percQ2.5[6]} else {0}, if (sub$percQ2.5[6]>0) {sub$percQ2.5[6]} else {0}, if (sub$percQ97.5[6]>0) {sub$percQ97.5[6]} else {0}, if (sub$percQ97.5[6]>0) {sub$percQ97.5[6]} else {0}),y=(yyy1+7), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[6]>0) {sub$percQ5[6]} else {0}, if (sub$percQ5[6]>0) {sub$percQ5[6]} else {0}, if (sub$percQ95[6]>0) {sub$percQ95[6]} else {0}, if (sub$percQ95[6]>0) {sub$percQ95[6]} else {0}),y=(yyy2+7), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[6]>0) {sub$percQ10[6]} else {0}, if (sub$percQ10[6]>0) {sub$percQ10[6]} else {0}, if (sub$percQ90[6]>0) {sub$percQ90[6]} else {0}, if (sub$percQ90[6]>0) {sub$percQ90[6]} else {0}),y=(yyy3+7), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[7]<0) {sub$percQ2.5[7]} else {0}, if (sub$percQ2.5[7]<0) {sub$percQ2.5[7]} else {0}, if (sub$percQ97.5[7]<0) {sub$percQ97.5[7]} else {0}, if (sub$percQ97.5[7]<0) {sub$percQ97.5[7]} else {0}),y=(yyy1+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[7]<0) {sub$percQ5[7]} else {0}, if (sub$percQ5[7]<0) {sub$percQ5[7]} else {0}, if (sub$percQ95[7]<0) {sub$percQ95[7]} else {0}, if (sub$percQ95[7]<0) {sub$percQ95[7]} else {0}),y=(yyy2+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[7]<0) {sub$percQ10[7]} else {0}, if (sub$percQ10[7]<0) {sub$percQ10[7]} else {0}, if (sub$percQ90[7]<0) {sub$percQ90[7]} else {0}, if (sub$percQ90[7]<0) {sub$percQ90[7]} else {0}),y=(yyy3+6), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[7]>0) {sub$percQ2.5[7]} else {0}, if (sub$percQ2.5[7]>0) {sub$percQ2.5[7]} else {0}, if (sub$percQ97.5[7]>0) {sub$percQ97.5[7]} else {0}, if (sub$percQ97.5[7]>0) {sub$percQ97.5[7]} else {0}),y=(yyy1+6), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[7]>0) {sub$percQ5[7]} else {0}, if (sub$percQ5[7]>0) {sub$percQ5[7]} else {0}, if (sub$percQ95[7]>0) {sub$percQ95[7]} else {0}, if (sub$percQ95[7]>0) {sub$percQ95[7]} else {0}),y=(yyy2+6), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[7]>0) {sub$percQ10[7]} else {0}, if (sub$percQ10[7]>0) {sub$percQ10[7]} else {0}, if (sub$percQ90[7]>0) {sub$percQ90[7]} else {0}, if (sub$percQ90[7]>0) {sub$percQ90[7]} else {0}),y=(yyy3+6), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[8]<0) {sub$percQ2.5[8]} else {0}, if (sub$percQ2.5[8]<0) {sub$percQ2.5[8]} else {0}, if (sub$percQ97.5[8]<0) {sub$percQ97.5[8]} else {0}, if (sub$percQ97.5[8]<0) {sub$percQ97.5[8]} else {0}),y=(yyy1+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[8]<0) {sub$percQ5[8]} else {0}, if (sub$percQ5[8]<0) {sub$percQ5[8]} else {0}, if (sub$percQ95[8]<0) {sub$percQ95[8]} else {0}, if (sub$percQ95[8]<0) {sub$percQ95[8]} else {0}),y=(yyy2+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[8]<0) {sub$percQ10[8]} else {0}, if (sub$percQ10[8]<0) {sub$percQ10[8]} else {0}, if (sub$percQ90[8]<0) {sub$percQ90[8]} else {0}, if (sub$percQ90[8]<0) {sub$percQ90[8]} else {0}),y=(yyy3+5), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[8]>0) {sub$percQ2.5[8]} else {0}, if (sub$percQ2.5[8]>0) {sub$percQ2.5[8]} else {0}, if (sub$percQ97.5[8]>0) {sub$percQ97.5[8]} else {0}, if (sub$percQ97.5[8]>0) {sub$percQ97.5[8]} else {0}),y=(yyy1+5), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[8]>0) {sub$percQ5[8]} else {0}, if (sub$percQ5[8]>0) {sub$percQ5[8]} else {0}, if (sub$percQ95[8]>0) {sub$percQ95[8]} else {0}, if (sub$percQ95[8]>0) {sub$percQ95[8]} else {0}),y=(yyy2+5), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[8]>0) {sub$percQ10[8]} else {0}, if (sub$percQ10[8]>0) {sub$percQ10[8]} else {0}, if (sub$percQ90[8]>0) {sub$percQ90[8]} else {0}, if (sub$percQ90[8]>0) {sub$percQ90[8]} else {0}),y=(yyy3+5), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[9]<0) {sub$percQ2.5[9]} else {0}, if (sub$percQ2.5[9]<0) {sub$percQ2.5[9]} else {0}, if (sub$percQ97.5[9]<0) {sub$percQ97.5[9]} else {0}, if (sub$percQ97.5[9]<0) {sub$percQ97.5[9]} else {0}),y=(yyy1+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[9]<0) {sub$percQ5[9]} else {0}, if (sub$percQ5[9]<0) {sub$percQ5[9]} else {0}, if (sub$percQ95[9]<0) {sub$percQ95[9]} else {0}, if (sub$percQ95[9]<0) {sub$percQ95[9]} else {0}),y=(yyy2+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[9]<0) {sub$percQ10[9]} else {0}, if (sub$percQ10[9]<0) {sub$percQ10[9]} else {0}, if (sub$percQ90[9]<0) {sub$percQ90[9]} else {0}, if (sub$percQ90[9]<0) {sub$percQ90[9]} else {0}),y=(yyy3+4), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[9]>0) {sub$percQ2.5[9]} else {0}, if (sub$percQ2.5[9]>0) {sub$percQ2.5[9]} else {0}, if (sub$percQ97.5[9]>0) {sub$percQ97.5[9]} else {0}, if (sub$percQ97.5[9]>0) {sub$percQ97.5[9]} else {0}),y=(yyy1+4), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[9]>0) {sub$percQ5[9]} else {0}, if (sub$percQ5[9]>0) {sub$percQ5[9]} else {0}, if (sub$percQ95[9]>0) {sub$percQ95[9]} else {0}, if (sub$percQ95[9]>0) {sub$percQ95[9]} else {0}),y=(yyy2+4), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[9]>0) {sub$percQ10[9]} else {0}, if (sub$percQ10[9]>0) {sub$percQ10[9]} else {0}, if (sub$percQ90[9]>0) {sub$percQ90[9]} else {0}, if (sub$percQ90[9]>0) {sub$percQ90[9]} else {0}),y=(yyy3+4), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[10]<0) {sub$percQ2.5[10]} else {0}, if (sub$percQ2.5[10]<0) {sub$percQ2.5[10]} else {0}, if (sub$percQ97.5[10]<0) {sub$percQ97.5[10]} else {0}, if (sub$percQ97.5[10]<0) {sub$percQ97.5[10]} else {0}),y=(yyy1+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[10]<0) {sub$percQ5[10]} else {0}, if (sub$percQ5[10]<0) {sub$percQ5[10]} else {0}, if (sub$percQ95[10]<0) {sub$percQ95[10]} else {0}, if (sub$percQ95[10]<0) {sub$percQ95[10]} else {0}),y=(yyy2+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[10]<0) {sub$percQ10[10]} else {0}, if (sub$percQ10[10]<0) {sub$percQ10[10]} else {0}, if (sub$percQ90[10]<0) {sub$percQ90[10]} else {0}, if (sub$percQ90[10]<0) {sub$percQ90[10]} else {0}),y=(yyy3+3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[10]>0) {sub$percQ2.5[10]} else {0}, if (sub$percQ2.5[10]>0) {sub$percQ2.5[10]} else {0}, if (sub$percQ97.5[10]>0) {sub$percQ97.5[10]} else {0}, if (sub$percQ97.5[10]>0) {sub$percQ97.5[10]} else {0}),y=(yyy1+3), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[10]>0) {sub$percQ5[10]} else {0}, if (sub$percQ5[10]>0) {sub$percQ5[10]} else {0}, if (sub$percQ95[10]>0) {sub$percQ95[10]} else {0}, if (sub$percQ95[10]>0) {sub$percQ95[10]} else {0}),y=(yyy2+3), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[10]>0) {sub$percQ10[10]} else {0}, if (sub$percQ10[10]>0) {sub$percQ10[10]} else {0}, if (sub$percQ90[10]>0) {sub$percQ90[10]} else {0}, if (sub$percQ90[10]>0) {sub$percQ90[10]} else {0}),y=(yyy3+3), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[11]<0) {sub$percQ2.5[11]} else {0}, if (sub$percQ2.5[11]<0) {sub$percQ2.5[11]} else {0}, if (sub$percQ97.5[11]<0) {sub$percQ97.5[11]} else {0}, if (sub$percQ97.5[11]<0) {sub$percQ97.5[11]} else {0}),y=(yyy1+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[11]<0) {sub$percQ5[11]} else {0}, if (sub$percQ5[11]<0) {sub$percQ5[11]} else {0}, if (sub$percQ95[11]<0) {sub$percQ95[11]} else {0}, if (sub$percQ95[11]<0) {sub$percQ95[11]} else {0}),y=(yyy2+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[11]<0) {sub$percQ10[11]} else {0}, if (sub$percQ10[11]<0) {sub$percQ10[11]} else {0}, if (sub$percQ90[11]<0) {sub$percQ90[11]} else {0}, if (sub$percQ90[11]<0) {sub$percQ90[11]} else {0}),y=(yyy3+2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[11]>0) {sub$percQ2.5[11]} else {0}, if (sub$percQ2.5[11]>0) {sub$percQ2.5[11]} else {0}, if (sub$percQ97.5[11]>0) {sub$percQ97.5[11]} else {0}, if (sub$percQ97.5[11]>0) {sub$percQ97.5[11]} else {0}),y=(yyy1+2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[11]>0) {sub$percQ5[11]} else {0}, if (sub$percQ5[11]>0) {sub$percQ5[11]} else {0}, if (sub$percQ95[11]>0) {sub$percQ95[11]} else {0}, if (sub$percQ95[11]>0) {sub$percQ95[11]} else {0}),y=(yyy2+2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[11]>0) {sub$percQ10[11]} else {0}, if (sub$percQ10[11]>0) {sub$percQ10[11]} else {0}, if (sub$percQ90[11]>0) {sub$percQ90[11]} else {0}, if (sub$percQ90[11]>0) {sub$percQ90[11]} else {0}),y=(yyy3+2), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[12]<0) {sub$percQ2.5[12]} else {0}, if (sub$percQ2.5[12]<0) {sub$percQ2.5[12]} else {0}, if (sub$percQ97.5[12]<0) {sub$percQ97.5[12]} else {0}, if (sub$percQ97.5[12]<0) {sub$percQ97.5[12]} else {0}),y=(yyy1+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[12]<0) {sub$percQ5[12]} else {0}, if (sub$percQ5[12]<0) {sub$percQ5[12]} else {0}, if (sub$percQ95[12]<0) {sub$percQ95[12]} else {0}, if (sub$percQ95[12]<0) {sub$percQ95[12]} else {0}),y=(yyy2+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[12]<0) {sub$percQ10[12]} else {0}, if (sub$percQ10[12]<0) {sub$percQ10[12]} else {0}, if (sub$percQ90[12]<0) {sub$percQ90[12]} else {0}, if (sub$percQ90[12]<0) {sub$percQ90[12]} else {0}),y=(yyy3+1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[12]>0) {sub$percQ2.5[12]} else {0}, if (sub$percQ2.5[12]>0) {sub$percQ2.5[12]} else {0}, if (sub$percQ97.5[12]>0) {sub$percQ97.5[12]} else {0}, if (sub$percQ97.5[12]>0) {sub$percQ97.5[12]} else {0}),y=(yyy1+1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[12]>0) {sub$percQ5[12]} else {0}, if (sub$percQ5[12]>0) {sub$percQ5[12]} else {0}, if (sub$percQ95[12]>0) {sub$percQ95[12]} else {0}, if (sub$percQ95[12]>0) {sub$percQ95[12]} else {0}),y=(yyy2+1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[12]>0) {sub$percQ10[12]} else {0}, if (sub$percQ10[12]>0) {sub$percQ10[12]} else {0}, if (sub$percQ90[12]>0) {sub$percQ90[12]} else {0}, if (sub$percQ90[12]>0) {sub$percQ90[12]} else {0}),y=(yyy3+1), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[13]<0) {sub$percQ2.5[13]} else {0}, if (sub$percQ2.5[13]<0) {sub$percQ2.5[13]} else {0}, if (sub$percQ97.5[13]<0) {sub$percQ97.5[13]} else {0}, if (sub$percQ97.5[13]<0) {sub$percQ97.5[13]} else {0}),y=(yyy1+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[13]<0) {sub$percQ5[13]} else {0}, if (sub$percQ5[13]<0) {sub$percQ5[13]} else {0}, if (sub$percQ95[13]<0) {sub$percQ95[13]} else {0}, if (sub$percQ95[13]<0) {sub$percQ95[13]} else {0}),y=(yyy2+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[13]<0) {sub$percQ10[13]} else {0}, if (sub$percQ10[13]<0) {sub$percQ10[13]} else {0}, if (sub$percQ90[13]<0) {sub$percQ90[13]} else {0}, if (sub$percQ90[13]<0) {sub$percQ90[13]} else {0}),y=(yyy3+0), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[13]>0) {sub$percQ2.5[13]} else {0}, if (sub$percQ2.5[13]>0) {sub$percQ2.5[13]} else {0}, if (sub$percQ97.5[13]>0) {sub$percQ97.5[13]} else {0}, if (sub$percQ97.5[13]>0) {sub$percQ97.5[13]} else {0}),y=(yyy1+0), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[13]>0) {sub$percQ5[13]} else {0}, if (sub$percQ5[13]>0) {sub$percQ5[13]} else {0}, if (sub$percQ95[13]>0) {sub$percQ95[13]} else {0}, if (sub$percQ95[13]>0) {sub$percQ95[13]} else {0}),y=(yyy2+0), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[13]>0) {sub$percQ10[13]} else {0}, if (sub$percQ10[13]>0) {sub$percQ10[13]} else {0}, if (sub$percQ90[13]>0) {sub$percQ90[13]} else {0}, if (sub$percQ90[13]>0) {sub$percQ90[13]} else {0}),y=(yyy3+0), col = "dodgerblue", border = NA)

polygon(x=c(if (sub$percQ2.5[14]<0) {sub$percQ2.5[14]} else {0}, if (sub$percQ2.5[14]<0) {sub$percQ2.5[14]} else {0}, if (sub$percQ97.5[14]<0) {sub$percQ97.5[14]} else {0}, if (sub$percQ97.5[14]<0) {sub$percQ97.5[14]} else {0}),y=(yyy1-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ5[14]<0) {sub$percQ5[14]} else {0}, if (sub$percQ5[14]<0) {sub$percQ5[14]} else {0}, if (sub$percQ95[14]<0) {sub$percQ95[14]} else {0}, if (sub$percQ95[14]<0) {sub$percQ95[14]} else {0}),y=(yyy2-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ10[14]<0) {sub$percQ10[14]} else {0}, if (sub$percQ10[14]<0) {sub$percQ10[14]} else {0}, if (sub$percQ90[14]<0) {sub$percQ90[14]} else {0}, if (sub$percQ90[14]<0) {sub$percQ90[14]} else {0}),y=(yyy3-1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$percQ2.5[14]>0) {sub$percQ2.5[14]} else {0}, if (sub$percQ2.5[14]>0) {sub$percQ2.5[14]} else {0}, if (sub$percQ97.5[14]>0) {sub$percQ97.5[14]} else {0}, if (sub$percQ97.5[14]>0) {sub$percQ97.5[14]} else {0}),y=(yyy1-1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ5[14]>0) {sub$percQ5[14]} else {0}, if (sub$percQ5[14]>0) {sub$percQ5[14]} else {0}, if (sub$percQ95[14]>0) {sub$percQ95[14]} else {0}, if (sub$percQ95[14]>0) {sub$percQ95[14]} else {0}),y=(yyy2-1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$percQ10[14]>0) {sub$percQ10[14]} else {0}, if (sub$percQ10[14]>0) {sub$percQ10[14]} else {0}, if (sub$percQ90[14]>0) {sub$percQ90[14]} else {0}, if (sub$percQ90[14]>0) {sub$percQ90[14]} else {0}),y=(yyy3-1), col = "dodgerblue", border = NA)

abline(h=7.5,lty=2, col="grey60")

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