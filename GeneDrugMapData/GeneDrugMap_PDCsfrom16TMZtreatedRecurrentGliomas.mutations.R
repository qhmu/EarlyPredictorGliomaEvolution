rm(list=ls())
dt = read.delim('~/Documents/Research/Longitudinal/GeneDrugMap_PDCsfrom16TMZtreatedRecurrentGliomas.mutations.txt')
dt = dt[order(dt$TMZ.AUC, decreasing = T),]
rownames(dt) = dt$ID
dt1 = as.data.frame(t(dt[,5:24]))
dt2 = dt[,1:4]

library(pheatmap)
pheatmap(dt1[c('IDH1','TP53','ATRX','PTEN','EGFR','NF1','PIK3CA','PIK3R1','MLH1',
               'EGFRamp',"PDGFRAamp",'PTENdel','CDKN2Adel','CDK4amp','MDM2amp','MDM4amp',"RB1del","MYCgain",
               'FGFR3.TACC3'),], annotation_col = subset(dt2, select = 'TMZ.AUC'), border_color = 'white',
         cluster_rows = F, cluster_cols = F,legend = F, color = c('#cccccc','red'), na_col = '#888888', 
         fontsize_row = 6, fontsize_col = 6)

#AUCs
dauc = read.delim('~/Documents/Research/Longitudinal/GeneDrugMap_PDCs.AUCs.txt', check.names = F)
dauc = dauc[dauc$`Sample ID` %in% dt$ID,]
names(dt)[4] = 'Temozolomide'
dauc = merge( dt[,c(1,4)],dauc, by.y = 'Sample ID', by.x = 'ID')
rownames(dauc) = dauc$ID
dauc=dauc[,-3]
## z score
# for (i in 2:ncol(dauc)){
#   dauc[,i] = scale(dauc[,i])[,1]
# }
#\z score
#quantile norm
# library(preprocessCore)
# dauc.qt <-normalize.quantiles(as.matrix(dauc[,-1]))
# dauc1 = as.data.frame( dauc.qt)
# dauc1$ID = dauc$ID
# dauc1 = dauc1[,c(ncol(dauc1), 1:(ncol(dauc1)-1))]
# names(dauc1) = names(dauc)
# rownames(dauc1) = rownames(dauc);
# dauc = dauc1
#\quantile norm

ix = apply(dauc['P43.T',-1],2,mean)
dauc = dauc[,c(1,which(!is.na(ix))+1)]
ix = as.numeric(dauc['P43.T',-1])
ix = order(ix, names(dauc)[-1],decreasing = F)
dauc = dauc[,c(1,(ix+1))]

library(reshape2)
m = melt(dauc, id.vars = 'ID')
drgs = c("Temozolomide","Imatinib (Gleevec)","Olaparib (AZD2281)")
m$drg = ifelse(m$variable %in% drgs,'MYC','no')
# lvls = levels(m$variable)[c(1,3,2,4:56)]
# m$variable = factor(m$variable, levels = lvls)
#m = m[!m$variable %in% c('Dacomitinib (PF299804_PF-00299804)','Afatinib (BIBW2992)'),]
ggplot(m, aes(x = variable, y = value, color = drg))+
  geom_boxplot(outlier.shape = NA, show.legend = F)+
  geom_point(data = m[m$ID=='P43.T',], mapping = aes(x =variable, y = value ), color = 'red', pch ='*',size=5)+
  #coord_flip()+
  scale_color_manual(values = c('tomato','black'))+
  labs(x = '', y ='AUC')+
  theme_bw()+
  theme(axis.text.x = element_text(color = 'black', size= 6, angle = 90, hjust=1))+
  theme(axis.text.y = element_text(color = 'black', size= 7),
        panel.border = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(colour = "black",size=0.4))
drgs = #c('ID',"TMZ.AUC","Ibrutinib","Imatinib..Gleevec.","ABT.888..Veliparib.","Olaparib..AZD2281.")
dauc=dauc[,names(dauc) %in% drgs]
library(reshape2)
m = melt(dauc, id.vars = 'ID')
#m$variable=factor(m$variable, levels = c("TMZ.AUC","Imatinib..Gleevec.","Ibrutinib","ABT.888..Veliparib.","Olaparib..AZD2281.") )
library(ggbeeswarm)
ggplot(m, aes(x = variable, y = value))+
  #geom_boxplot(outlier.shape = NA, aes(fill = variable), show.legend = F, width=0.75)+
  geom_quasirandom(width = 0.25, color = '#999999')+
  geom_point(data = m[m$ID=='P43.T',], mapping = aes(x =variable, y = value ), color = 'red', pch =18,size=4)+
  theme_classic()+theme(axis.text = element_text(color ='black'),axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels =c('TMZ','Imatinib','Ibrutinib','Veliparib','Olaparib'))+
  labs(x = '', y ='AUC')


