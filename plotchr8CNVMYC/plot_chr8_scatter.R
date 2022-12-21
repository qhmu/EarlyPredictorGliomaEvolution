setwd('~/Documents/Research/Longitudinal/')
myc = c(128745680,128757197)
library(data.table)
#PS118
dp = fread('PS118.T04BIR.chr8.depth')
dp = dp[dp$V3>10 & dp$V4>10 & dp$V5>10,]

dp$myc = ifelse(dp$V2>=myc[1] & dp$V2<=myc[2],'y','n')
bnr = 46220882
inr = 54051287
rnr = 54492541
sizeFactor1 = inr/bnr
sizeFactor2 = rnr/bnr


par(cex=0.8, mai=c(0.5,0.35,0.05,0.1),fig=c(0,1,0.,0.63),mgp=c(1.8,0.7,0))
dp1 = dp[sample(1:nrow(dp),20000),]
plot(dp1$V2/1e6, log2(dp1$V5/dp1$V3/sizeFactor2), pch = 20, cex = 0.1, col=rgb(0,0,0,0.5), 
     xlab = 'chr8', ylab = '', ylim=c(-2,2), yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp1$V2[dp1$myc=='y']/1e6, log2(dp1$V5[dp1$myc=='y']/dp1$V3[dp1$myc=='y']/sizeFactor2), pch = 20, cex = 0.2, col=rgb(1,0,0,0.5), ylim=c(-2,2)); 
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')

par(cex=0.8, mai=c(0.05,0.35,0.05,0.1),fig=c(0,1,0.66,1),mgp=c(1.8,0.7,0), new = T)
dp1 = dp[sample(1:nrow(dp),20000),]
plot(dp1$V2/1e6, log2(dp1$V4/dp1$V3/sizeFactor1), pch = 20, cex = 0.1, col=rgb(0,0,0,0.5), 
      ylab = '', ylim=c(-2,2), xaxt = 'n', yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp1$V2[dp1$myc=='y']/1e6, log2(dp1$V4[dp1$myc=='y']/dp1$V3[dp1$myc=='y']/sizeFactor1), pch = 20, cex = 0.2, col=rgb(1,0,0,0.5), ylim=c(-2,2)); 
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')

##PSX01
bnr = 71125113
inr = 144021252
rnr = 189362975
dp = fread('PSX01.1149P742.chr8.depth')
dp = dp[dp$V3>10 & dp$V4>10 & dp$V5>10,]

dp$myc = ifelse(dp$V2>=myc[1] & dp$V2<=myc[2],'y','n')

sizeFactor1 = inr/bnr
sizeFactor2 = rnr/bnr

par(cex=0.8, mai=c(0.5,0.35,0.05,0.1),fig=c(0,1,0.,0.63),mgp=c(1.8,0.7,0))
dp1 = dp[sample(1:nrow(dp),20000),]
plot(dp1$V2/1e6, log2(dp1$V5/dp1$V3/sizeFactor2), pch = 20, cex = 0.1, col=rgb(0,0,0,0.5), 
     xlab = 'chr8', ylab = '', ylim=c(-2,2), yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp1$V2[dp1$myc=='y']/1e6, log2(dp1$V5[dp1$myc=='y']/dp1$V3[dp1$myc=='y']/sizeFactor2), pch = 20, cex = 0.5, col=rgb(1,0,0,0.5), ylim=c(-2,2)); 
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')

par(cex=0.8, mai=c(0.05,0.35,0.05,0.1),fig=c(0,1,0.66,1),mgp=c(1.8,0.7,0), new = T)
dp1 = dp[sample(1:nrow(dp),20000),]
plot(dp1$V2/1e6, log2(dp1$V4/dp1$V3/sizeFactor1), pch = 20, cex = 0.1, col=rgb(0,0,0,0.5), 
     ylab = '', ylim=c(-2,2), xaxt = 'n', yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp1$V2[dp1$myc=='y']/1e6, log2(dp1$V4[dp1$myc=='y']/dp1$V3[dp1$myc=='y']/sizeFactor1), pch = 20, cex = 0.5, col=rgb(1,0,0,0.5), ylim=c(-2,2)); 
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')

#PSX02
bnr = 44587661
inr = 186748306
rnr = 148599754
dp = fread('PSX02.p387p890.chr8.depth')
dp = dp[dp$V3>10 & dp$V4>10 & dp$V5>10,]

dp$myc = ifelse(dp$V2>=myc[1] & dp$V2<=myc[2],'y','n')

sizeFactor1 = inr/bnr
sizeFactor2 = rnr/bnr
dp1 = dp[sample(1:nrow(dp),20000),]

par(cex=0.8, mai=c(0.5,0.35,0.05,0.1),fig=c(0,1,0.,0.63),mgp=c(1.8,0.7,0))

plot(dp1$V2/1e6, log2(dp1$V5/dp1$V3/sizeFactor2), pch = 20, cex = 0.1, col=rgb(0,0,0,0.5), 
     xlab = 'chr8', ylab = '', ylim=c(-2,2), yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp1$V2[dp1$myc=='y']/1e6, log2(dp1$V5[dp1$myc=='y']/dp1$V3[dp1$myc=='y']/sizeFactor2), pch = 20, cex = 0.5, col=rgb(1,0,0,0.5), ylim=c(-2,2)); 
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')

par(cex=0.8, mai=c(0.05,0.35,0.05,0.1),fig=c(0,1,0.66,1),mgp=c(1.8,0.7,0), new = T)
dp1 = dp[sample(1:nrow(dp),20000),]
plot(dp1$V2/1e6, log2(dp1$V4/dp1$V3/sizeFactor1), pch = 20, cex = 0.1, col=rgb(0,0,0,0.5), 
     ylab = '', ylim=c(-2,2), xaxt = 'n', yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp1$V2[dp1$myc=='y']/1e6, log2(dp1$V4[dp1$myc=='y']/dp1$V3[dp1$myc=='y']/sizeFactor1), pch = 20, cex = 0.5, col=rgb(1,0,0,0.5), ylim=c(-2,2)); 
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')


#PS146
bnr = 648471031
inr = 1055598788
rnr = 1032498332
dp = fread('PS146_T32BIRG.chr8.min10.depth')
#dp = dp[dp$V3>10 & dp$V4>10 & dp$V5>10,]

dp0 = dp[dp$V2>=myc[1] & dp$V2<=myc[2],]
dp0 = dp0[sample(1:nrow(dp0),4),]

sizeFactor1 = inr/bnr
sizeFactor2 = rnr/bnr

par(cex=0.8, mai=c(0.5,0.35,0.05,0.1),fig=c(0,1,0.,0.63),mgp=c(1.8,0.7,0))
dp1 = dp[sample(1:nrow(dp),20000),]
dp1$myc = ifelse(dp1$V2>=myc[1] & dp1$V2<=myc[2],'y','n')
plot(dp1$V2/1e6, log2(dp1$V5/dp1$V3/sizeFactor2), pch = 20, cex = 0.1, col=rgb(0,0,0,1/3), 
     xlab = 'chr8', ylab = '', ylim=c(-2,2), yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp0$V2/1e6, log2(dp0$V5/dp0$V3/sizeFactor2), pch = 20, cex = 0.5, col=rgb(1,0,0,0.5), ylim=c(-2,2)); 
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')

par(cex=0.8, mai=c(0.05,0.35,0.05,0.1),fig=c(0,1,0.66,1),mgp=c(1.8,0.7,0), new = T)
dp1 = dp[sample(1:nrow(dp),20000),]
plot(dp1$V2/1e6, log2(dp1$V4/dp1$V3/sizeFactor1), pch = 20, cex = 0.1, col=rgb(0,0,0,1/3), 
     ylab = '', ylim=c(-2,2), xaxt = 'n', yaxt = 'n'); 
axis(side = 2,at = c(-2,0,2), labels = c(-2,0,2))
points(dp0$V2/1e6, log2(dp0$V4/dp0$V3/sizeFactor1), pch = 20, cex = 0.5, col=rgb(1,0,0,0.5), ylim=c(-2,2));  
abline(h=0,lty=3,lwd=1.5,col='#a1d76a')

