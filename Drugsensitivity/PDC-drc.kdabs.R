library(drc)

library(reshape2)


cls = c("#b2df8a","#33a02c","#cab2d6","#6a3d9a","#fb9a99","#e31a1c")

inputfile<-"~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_5/PDC.drc.dat"
data<-read.table(inputfile,sep="\t",header=T,quote=NULL, stringsAsFactors = F)
data = data[data$Cell %in% c('T2-4','T2-4TR'),]
names(data) = names(dt)
dt = rbind(dt, data)
modelC3 <- drm(value~Conc,Cell, fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")), data=dt)
plot(modelC3, type = 'bar', ylim=c(0,1.1),
     col = cls[c(4,3,2,1,6,5)],
     lty=c(3,1,3,1,3,1),
     lwd = 2, legend = F, pch = c(20,17,20,17,20,17),cex=0.75,
     ylab="Response ", xlab="Concentration (µM)", broken = T)

#legend('topright',legend = c('T2-4TR','T2-4','MYCKD','control'), pch=c(20,17,20,17), col = cls,lty=c(1,3,1,3),bty='n', cex = 0.65)
