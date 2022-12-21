setwd('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/')

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(fgsea)
hmpathways <- gmtPathways("h.all.v7.0.symbols.gmt")
plotfgsea <- function(pathway, stats,fgseaRes, gene.set.name, class.name, posClass, negClass){
  stopifnot(!missing(pathway))
  stopifnot(!missing(stats))
  stopifnot(!missing(fgseaRes))
  stopifnot(gene.set.name %in% fgseaRes$pathway )
  stats = sort(stats, decreasing = T)
  metric.range <- c(min(stats), max(stats))
  gsea.enrichment.score= round(fgseaRes$ES[which(fgseaRes$pathway==gene.set.name)],2)
  gsea.normalized.enrichment.score= round(fgseaRes$NES[which(fgseaRes$pathway==gene.set.name)],2)
  gsea.p.value = fgseaRes$pval[which(fgseaRes$pathway==gene.set.name)]
  gsea.p.value = formatC(gsea.p.value, format = "e", digits = 2)
  gsea.fdr = fgseaRes$padj[which(fgseaRes$pathway==gene.set.name)]
  gsea.fdr = formatC(gsea.fdr, format = "e", digits = 2)
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * abs(statsAdj)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  enrichment.score.range <- c(min(toPlot$y), max(toPlot$y))
  
  def.par <- par(no.readonly = TRUE)
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.3, 0.2, 2))
  #layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(toPlot$x, toPlot$y, type = "l", 
       col = "#A4C361",lwd = 2,lty=1,pch=19, cex=0.25,
       xaxt = "n",xaxs = "i", xlab = "", ylab = "Enrichment score",
       ylim = enrichment.score.range,
       main = list(gene.set.name, font = 1, cex = 1),
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(stats) * 0.01, plot.coordinates[3] * 0.98,
         paste( "FDR:", gsea.fdr,#"p-value:", gsea.p.value,
                "\nNES:",gsea.normalized.enrichment.score, "\n"), adj = c(0, 0))
  } else {
    text(length(stats) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste( "FDR:", gsea.fdr,#"p-value:", gsea.p.value,
                "\nNES:",gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
  }
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(stats)))
  abline(v = toPlot$x, lwd = 0.75)
  
  par(mar = c(0, 5, 0, 2))
  rank.colors <- stats - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(stats)))
  box()
  text(length(stats) / 2, 0.7,
       labels = ifelse(!missing(class.name), class.name, gsea.template))
  text(length(stats) * 0.01, 0.7, ifelse(!missing(posClass), posClass, 'Positive'), adj = c(0, NA))
  text(length(stats) * 0.99, 0.7, ifelse(!missing(negClass), negClass, 'Negative'), adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  
  #rank.metric <- rle(round(unname(stats), digits = 2))
  #plot(stats, type = "n", xaxs = "i",
  #     xlab = "Rank in ordered gene list", xlim = c(0, length(stats)),
  #     ylim = metric.range, yaxs = "i",
  #     ylab = 'Ranking metric',
  #ylab = if(gsea.metric == "None") {"Ranking metric"} else {gsea.metric},
  #     panel.first = abline(h = seq(metric.range[1] / 2,
  #                                  metric.range[2] - metric.range[1] / 4,
  #                                  metric.range[2] / 2), col = "gray95", lty = 2))
  
  #barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",xaxt='n',
  #        xlab = "", xlim = c(0, length(stats)),
  #        ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
  #        #ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric),
  #        ylab = "",space = 0, add = TRUE)
  #v = (sum(rank.metric$lengths[rank.metric$values>0]) + sum(rank.metric$lengths[rank.metric$values>=0]))/2
  #v = round(v,0)
  #abline(v=v,lty=2, col = 'lightgrey')
  #text(x = v,y=0, labels = paste0('Zero cross at ',v))
  #box()
  
  # Reset to default
  par(def.par)
}


gei = read.delim('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/revision/panglioma544.initial.RPKM20220830qtnorm.txt')
gei0 = gei

ger = read.delim('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/revision/panglioma544.recurrence.RPKM20220830qtnorm.txt')
ger0 = ger

batch2 = read.delim('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/glass2022_newRNAseq_153samples.txt')
zscale = T
#------------
if (zscale ==T){
  #initial 
  m1 = gei[,names(gei) %in% batch2$pairID]
  m1z = as.data.frame(t(scale(t(log2(m1+1)))))
  
  m2= gei[,!names(gei) %in% batch2$pairID]
  m2z = as.data.frame(t(scale(t(log2(m2+1)))))
  
  mm12z = cbind(m1z, m2z)
  mm12z = mm12z[apply(mm12z, 1, function(x) sum(is.na(x)))==0,]
  gei = mm12z 
  #recurrence
  m1r = ger[,names(ger) %in% batch2$pairID]
  m1rz = as.data.frame(t(scale(t(log2(m1r+1)))))
  
  m2r= ger[,!names(ger) %in% batch2$pairID]
  m2rz = as.data.frame(t(scale(t(log2(m2r+1)))))
  
  mm12rz = cbind(m1rz, m2rz)
  mm12rz = mm12rz[apply(mm12rz, 1, function(x) sum(is.na(x)))==0,]
  ger = mm12rz 
  
}else{
  gei = log2(gei+1)
  ger = log2(ger+1)
}

#--------------
cl = read.delim('../revision/panglioma544.clinical.information.txt')
cl$ID1 = paste0(cl$Patient.ID,"_I")
cl$ID2 = paste0(cl$Patient.ID,"_R")
cl$Subtype = paste(cl$IDH.mutation, cl$X1p.19q.co.deletion, sep = "-")


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#IDHwt
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cl_wt = cl[which(cl$Subtype=='IDHwt-noncodel' & cl$TMZ.Treatment.Before.Recurrence==1),]
prehm_wt = cl_wt$ID1[cl_wt$Hypermutation.of.Recurrent.Tumor=='Yes']
prenhm_wt = cl_wt$ID1[cl_wt$Hypermutation.of.Recurrent.Tumor=='No']

hm_wt = cl_wt$ID2[cl_wt$Hypermutation.of.Recurrent.Tumor=='Yes']
nhm_wt = cl_wt$ID2[cl_wt$Hypermutation.of.Recurrent.Tumor=='No']

plotgeneiwt <-function(gene = 'MGMT', offset = 2, mytest = 'wilcox.test'){ #examine a single gene
  x = as.numeric(gei[gene,names(gei)%in% prehm_wt]) #Pre-HM
  y = as.numeric(gei[gene,names(gei)%in% prenhm_wt]) #pre-NHM
  df = data.frame(var = c(rep('pre-HM', length(x)),rep('pre-NHM', length(y)) ), val = c(x,y))
  df$SID = c(names(gei)[names(gei)%in% prehm_wt], names(gei)[names(gei)%in% prenhm_wt])
  df$isnew = ifelse(df$SID %in% batch2$pairID, 'isnew','notnew')
  df$CaseID = substr(df$SID,1,5)
  df$Timepoint = substr(df$SID, 7,7)
  
  xr = as.numeric(ger[gene,names(ger)%in% hm_wt]) #HM
  yr = as.numeric(ger[gene,names(ger)%in% nhm_wt]) #NHM
  dfr = data.frame(var = c(rep('HM', length(xr)),rep('NHM', length(yr)) ), val = c(xr,yr))
  dfr$SID = c(names(ger)[names(ger)%in% hm_wt], names(ger)[names(ger)%in% nhm_wt])
  dfr$isnew = ifelse(dfr$SID %in% batch2$pairID, 'isnew','notnew')
  dfr$CaseID = substr(dfr$SID,1,5)
  dfr$Timepoint = substr(dfr$SID, 7,7)
  
  df = rbind(df, dfr)
  df$var = factor(df$var, levels = c('pre-HM','HM','pre-NHM','NHM'))
  df$hmstatus = ifelse(df$var %in% c('HM','pre-HM'),'H','N')
  p<- 
    ggplot(df, aes(x = var , y = val, color = hmstatus))+
    geom_boxplot(width = 0.4, show.legend = F, outlier.shape = NA)+
    geom_line(aes(group = CaseID), color = '#ddddddaa', size=1)+
    geom_jitter( show.legend = F, size=0.5,width = 0.1)+
    stat_compare_means(label = 'p.format', method = mytest, show.legend = F,
                       comparisons = list(c('pre-HM','HM'),c('NHM','pre-NHM'),c('NHM','HM'),c('pre-HM','pre-NHM')))+
    scale_color_manual(values = c('#e41a1c','#377eb8'))+
    theme_classic()+theme(axis.text = element_text(color = 'black'))+
    lims(y=c(NA, offset*max(df$val,na.rm = T)))+
    labs(x = '', y = paste0('Scaled ',gene, ' expression'), color = '',title = 'IDHwt')
  print(p)
}

plotgeneiwt(gene = 'MGMT')
plotgeneiwt(gene = 'MKI67')



#MGMT methylation and expression and hypermutation
mgmt = read.delim('~/Dropbox/communter//Rev1/MYCbySubtype/MYC_methylation_subtype.txt', 
                      stringsAsFactors = F,na.strings = c("NA","#N/A","NaN"), row.names = 1)
x = as.numeric(gei['MGMT',names(gei)%in% prehm_wt]) #Pre-HM
y = as.numeric(gei['MGMT',names(gei)%in% prenhm_wt]) #pre-NHM
df = data.frame(var = c(rep('pre-HM', length(x)),rep('pre-NHM', length(y)) ), val = c(x,y))
df$SID = c(names(gei)[names(gei) %in%prehm_wt ], names(gei)[names(gei) %in%prenhm_wt ])
df$MGMTmethy = mgmt$MGMT_methylation_I[match(substr(df$SID,1,5),rownames(mgmt))]
df$MGMTmethy=ifelse(is.na(df$MGMTmethy),'unknown',ifelse(df$MGMTmethy==1,'methylated','unmethylated'))
df$MGMTmethy = factor(df$MGMTmethy, levels = c('methylated','unmethylated','unknown'))
df$TMZ = cl_wt$TMZ.Treatment.Before.Recurrence[match(substr(df$SID,1,5), cl_wt$Patient.ID)]
df_wt_mgmt = df
ggplot(df[which(df$TMZ==1 ),], aes(x = MGMTmethy , y = val, color = var))+
  geom_boxplot(width = 0.5, show.legend = T, outlier.shape = NA)+
  geom_jitter( show.legend = T,position = position_dodge(width = 0.5))+
  scale_shape_manual(values = c(17,19,4))+
  stat_compare_means(label = 'p.format',show.legend = F,label.y.npc = 0.95)+
  scale_color_manual(values = c('#e41a1c','#377eb8','#cccccc'), guide = 'none')+
  theme_classic()+theme(axis.text = element_text(color = 'black'))+
  lims(y=c(NA, 1.8*max(df$val,na.rm = T)))+
  #facet_wrap(~MGMTmethy)+
  #scale_y_continuous(expand = c(0.1,1))+
  labs(x = '', y =  'Scaled MGMT expression', color = '')

#differential gene expression analysis
deg_wt_ini = data.frame(gene = rownames(gei), med_hm = 0,med_nhm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_wt_ini)){
  x = as.numeric(gei[i,names(gei)%in% prehm_wt]) #Pre-HM
  y = as.numeric(gei[i,names(gei)%in% prenhm_wt]) #pre-NHM
  deg_wt_ini$med_hm[i] = mean(x)
  deg_wt_ini$med_nhm[i] = mean(y)
  
  deg_wt_ini$log2fc[i] = deg_wt_ini$med_hm[i] - deg_wt_ini$med_nhm[i]
  
  #deg_wt_ini$pval[i] = t.test(x,y)$p.value
  #deg_wt_ini$tstat[i] = t.test(x,y)$statistic
  
  deg_wt_ini$pval[i] = wilcox.test(x,y)$p.value
  deg_wt_ini$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_wt_ini$type = ifelse(deg_wt_ini$pval>0.05, 'nosig',ifelse(deg_wt_ini$log2fc>0,'up','down'))

deg_wt_rec = data.frame(gene = rownames(ger), med_hm = 0,med_nhm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_wt_rec)){
  x = as.numeric(ger[i,names(ger)%in% hm_wt]) #HM
  y = as.numeric(ger[i,names(ger)%in% nhm_wt]) #NHM
  deg_wt_rec$med_hm[i] = mean(x)
  deg_wt_rec$med_nhm[i] = mean(y)
  
  deg_wt_rec$log2fc[i] = deg_wt_rec$med_hm[i] - deg_wt_rec$med_nhm[i]
  #deg_wt_rec$pval[i] = t.test(x,y)$p.value
  #deg_wt_rec$tstat[i] = t.test(x,y)$statistic
  
  deg_wt_rec$pval[i] = wilcox.test(x,y)$p.value
  deg_wt_rec$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_wt_rec$type = ifelse(deg_wt_rec$pval>0.05, 'nosig',ifelse(deg_wt_rec$log2fc>0,'up','down'))

deg_wt_hmpair = data.frame(gene = rownames(gei), med_prehm = 0,med_hm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_wt_hmpair)){
  x = as.numeric(gei[i,names(gei)%in% prehm_wt]) #pre-HM
  y = as.numeric(ger[i,names(ger)%in% hm_wt]) #HM
  deg_wt_hmpair$med_prehm[i] = mean(x)
  deg_wt_hmpair$med_hm[i] = mean(y)
  
  deg_wt_hmpair$log2fc[i] = deg_wt_hmpair$med_hm[i] - deg_wt_hmpair$med_prehm[i]
  #deg_wt_hmpair$pval[i] = t.test(x,y)$p.value
  #deg_wt_hmpair$tstat[i] = t.test(x,y)$statistic
  
  deg_wt_hmpair$pval[i] = wilcox.test(x,y)$p.value
  deg_wt_hmpair$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_wt_hmpair$type = ifelse(deg_wt_hmpair$pval>0.05, 'nosig',ifelse(deg_wt_hmpair$log2fc>0,'up','down'))

deg_wt_ini$timepoint = 'initial'; deg_wt_rec$timepoint = 'recurrence';deg_wt_hmpair$timepoint = 'pair'

deg_wt = deg_wt_ini[,c('gene','med_hm','med_nhm','log2fc','pval','tstat','type')]
names(deg_wt) = paste0(names(deg_wt),"_ini")
deg_wt$log2fc_rec = deg_wt_rec$log2fc
deg_wt$pval_rec = deg_wt_rec$pval
deg_wt$type_rec = deg_wt_rec$type
deg_wt$pval_pair = deg_wt_hmpair$pval
deg_wt$type_pair = deg_wt_hmpair$type

write.table(deg_wt, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.deg.ranksum.txt', row.names = F, quote = F, sep = "\t")

ggplot(deg_wt, aes(x = log2fc_ini, y = log2fc_rec))+
  #geom_point( alpha = 0.1)+
  geom_hex(binwidth = c(.03, 0.03),show.legend = T)+
  #stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200)+
  geom_hline(yintercept = 0,lty=2, col = '#999999')+
  geom_vline(xintercept = 0, lty = 2, col = '#999999')+
  scale_fill_viridis_c()+
  #scale_fill_continuous(low = "white", high = "dodgerblue4")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x = 'pre-HM versus pre-NHM', y = 'HM versus NHM', fill = '# genes')

ggplot(deg_wt, aes(x = sign(log2fc_ini)*-log10(pval_ini), y = sign(log2fc_rec)*-log10(pval_rec)))+
  #geom_point( alpha = 0.1)+
  geom_hex(binwidth = c(.25, 0.25),show.legend = T)+
  #stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200)+
  geom_hline(yintercept = 0,lty=2, col = '#999999')+
  geom_vline(xintercept = 0, lty = 2, col = '#999999')+
  scale_fill_viridis_c()+
  #scale_fill_continuous(low = "white", high = "dodgerblue4")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x = 'pre-HM versus pre-NHM', y = 'HM versus NHM', fill = '# genes')


cor.test(deg_wt$log2fc_ini, deg_wt$log2fc_rec,alternative = 'two.sided')$p.value
pheatmap(table(deg_wt$type_ini, deg_wt$type_rec),display_numbers = T,number_format = '%.0f', legend = F, cluster_rows = F, cluster_cols = F)

#fgsea
deg_wt_ini = deg_wt_ini[order(deg_wt_ini$tstat, decreasing = T),]

fcrnk = deg_wt_ini$tstat
names(fcrnk) = deg_wt_ini$gene  
fcrnk = fcrnk[!is.na(fcrnk)]
fgseaRes <- fgsea(pathways = hmpathways, stats = fcrnk, minSize=15,maxSize=500)

plotfgsea(pathway = hmpathways[["HALLMARK_MYC_TARGETS_V2"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MYC_TARGETS_V2',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')
plotfgsea(pathway = hmpathways[["HALLMARK_MITOTIC_SPINDLE"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MITOTIC_SPINDLE',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')

#
fgseaRes = fgseaRes[order(-fgseaRes$NES, decreasing = T),]
fgseaRes$pathway = factor(fgseaRes$pathway, levels = fgseaRes$pathway)
fgseaRes$pathway1 = gsub("HALLMARK_","", fgseaRes$pathway)
fgseaRes$pathway1  = factor(fgseaRes$pathway1, levels = fgseaRes$pathway1)
fgseaRes$pval[fgseaRes$pval< 1e-5] = 1e-5
ggplot(aes(x = pathway1, y = NES, fill = sign(NES)*-log10(pval)),data = fgseaRes[fgseaRes$pval<0.05,])+
  geom_bar(stat = 'identity', color = 'black') + 
  scale_fill_gradient2(low = 'blue',high='red')+
  coord_flip() + 
  theme_minimal()+labs(fill = 'Significance', x= '')+
  theme(axis.text.x = element_text(colour = "black"),#panel.grid = element_blank(),
        axis.text.y = element_text(colour = "black"))

gns2mark = deg_wt$gene[deg_wt$gene %in% c(hmpathways$HALLMARK_MYC_TARGETS_V1,hmpathways$HALLMARK_MYC_TARGETS_V2,
                                          hmpathways$HALLMARK_E2F_TARGETS,hmpathways$HALLMARK_G2M_CHECKPOINT,hmpathways$HALLMARK_MITOTIC_SPINDLE)];
gns2mark = c(gns2mark,'MGMT')
deg_wt$lab = ifelse(deg_wt$gene %in% gns2mark,deg_wt$gene, NA)

deg_wt_predictors_up = deg_wt[deg_wt$type_ini=='up' & deg_wt$type_rec=='up' & deg_wt$type_pair!='down',]
deg_wt_predictors_dn = deg_wt[deg_wt$type_ini=='down' & deg_wt$type_rec=='down' & deg_wt$type_pair!='up',]
deg_wt_predictors = rbind(deg_wt_predictors_up, deg_wt_predictors_dn)
deg_wt_predictors$pcomb = 1
for (i in 1:nrow(deg_wt_predictors)){
  deg_wt_predictors$pcomb[i] = pchisq( -2*sum(log(c(deg_wt_predictors$pval_ini[i],deg_wt_predictors$pval_rec[i]))), df=4, lower.tail=FALSE)
}
deg_wt_predictors$comb.padj = p.adjust(deg_wt_predictors$pcomb)
deg_wt_predictors$signedpcomb = sign(deg_wt_predictors$log2fc_ini)*-log10(deg_wt_predictors$pcomb)
deg_wt_predictors = deg_wt_predictors[order(deg_wt_predictors$signedpcomb, decreasing = F),]
deg_wt_predictors$idx = 1:nrow(deg_wt_predictors)
deg_wt_predictors$myctagrtes = ifelse(deg_wt_predictors$gene_ini %in% c(hmpathways$HALLMARK_MYC_TARGETS_V1,hmpathways$HALLMARK_MYC_TARGETS_V2),deg_wt_predictors$gene_ini,NA)
deg_wt_predictors$myctagrtes[deg_wt_predictors$gene_ini=='MKI67']='MKI67'
deg_wt_predictors$myctagrtes[deg_wt_predictors$gene_ini=='MGMT']='MGMT'

ggplot()+
  geom_bar(data = deg_wt_predictors, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_wt_predictors[!is.na(deg_wt_predictors$myctagrtes),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_wt_predictors[!is.na(deg_wt_predictors$myctagrtes),],
                  size=2.5,segment.size=0.25,max.iter =1000,max.overlaps = 30,
                  mapping = aes(x = idx, y = signedpcomb,label = lab))+
  theme_classic()+coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')


ggplot()+
  geom_bar(data = deg_wt_predictors, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_wt_predictors[!is.na(deg_wt_predictors$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_wt_predictors[!is.na(deg_wt_predictors$lab),],
                  size=3,segment.size=0.25,max.iter =1000,max.overlaps = 30,
             mapping = aes(x = idx, y = signedpcomb,label = lab))+
  theme_classic()+#coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')

deg_wt_predictors_pcrt = deg_wt_predictors[deg_wt_predictors$comb.padj<0.05,]
deg_wt_predictors_pcrt$idx = 1:nrow(deg_wt_predictors_pcrt)

ggplot()+
  geom_bar(data = deg_wt_predictors_pcrt, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_wt_predictors_pcrt,pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_wt_predictors_pcrt,
                  size=2.5,segment.size=0.1,max.iter =1000,max.overlaps = 10,nudge_x = -1,hjust=1,direction     = "y",
                  nudge_y = ifelse(deg_wt_predictors_pcrt$log2fc_rec>0,5,-15),
                  mapping = aes(x = idx, y = signedpcomb,label = gene_ini))+
  theme_classic()+coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')

deg_wt_predictors_pcrt$lab[deg_wt_predictors_pcrt$gene_ini=='CDT1']='CDT1'
ggplot()+
  geom_bar(data = deg_wt_predictors_pcrt, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_wt_predictors_pcrt[!is.na(deg_wt_predictors_pcrt$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_wt_predictors_pcrt[!is.na(deg_wt_predictors_pcrt$lab),],
                  size=2.5,segment.size=0.1,max.iter =1000,max.overlaps = 10,nudge_x = -1,hjust=1,direction     = "y",
                  nudge_y = ifelse(deg_wt_predictors_pcrt$log2fc_rec[!is.na(deg_wt_predictors_pcrt$lab)]>0,5,-5),
                  mapping = aes(x = idx, y = signedpcomb,label = lab))+
  theme_classic()+coord_flip()+lims(y = c(-9,9))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  labs(x = 'Genes', y = 'signed log10(P)', fill = '')


#deg_wt$type = ifelse(deg_wt$pval>0.05, 'nosig',ifelse(deg_wt$log2fc>0,'up','down'))

library(ggbeeswarm)
pos <- position_quasirandom()
deg_wt$pval_ini[deg_wt$pval_ini<1e-7]=1e-7
ggplot(deg_wt[!is.na(deg_wt$pval_ini) & deg_wt$log2fc_ini!=0,], aes(x = log2fc_ini, y = -log10(pval_ini), label = lab))+
  geom_point(aes( color = type_ini), size = 0.2, alpha = 0.7, show.legend = F)+
  geom_vline(xintercept = 0, lty=2)+
  geom_text_repel(size=2.,segment.size = 0.25, max.iter = 1e4,max.overlaps = 30,position = pos)+ #color = "white",  bg.color = "grey30", colorbg.r = 0.15
  scale_color_manual(values = c('blue','#dddddd','red'))+
  theme_classic()+#lims(x = c(-2.,2))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'log2 (pre-HM/pre-NHM)', y = '-log10(P)')

#the scatter plot is weak. do heatmap instead

# upgns = NULL
# for (pwy in fgseaRes$pathway[which(fgseaRes$padj<0.05 &fgseaRes$NES >0)]){
#   upgns = append(upgns, hmpathways[[as.character(pwy)]])
# }
# upgns = unique(upgns)
# upgns = upgns[upgns %in% deg_wt$gene[deg_wt$pval<0.05]]
# 
# dngns = NULL
# for (pwy in fgseaRes$pathway[which(fgseaRes$padj<0.05 &fgseaRes$NES <0)]){
#   dngns = append(dngns, hmpathways[[as.character(pwy)]])
# }
# dngns = unique(dngns)
# dngns = dngns[dngns %in% deg_wt$gene[deg_wt$pval<0.05]]

upgns = deg_wt$gene[which(deg_wt$type=='up' & deg_wt$log2fc>0.5)]
dngns = deg_wt$gene[which(deg_wt$type=='down' & deg_wt$log2fc< -0.5 )]

prenhm_wt2 = cl_wt$ID1[cl_wt$Hypermutation.of.Recurrent.Tumor=='No' & cl_wt$TMZ.Treatment.Before.Recurrence==1]
m1h = gei[upgns,names(gei) %in% prehm_wt]
m1n = gei[upgns,names(gei) %in% prenhm_wt2]

m2h = gei[dngns,names(gei) %in% prehm_wt]
m2n = gei[dngns,names(gei) %in% prenhm_wt2]

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
gns2mark =c('MKI67','PLK1','POLE','NUMA1','POLD1','PLK1',
            'DNMT1','MGMT','SERPINE1','DUSP1','CLU','TMEM50A','SDHD')
mat = rbind(cbind(m1h,m1n), cbind(m2h,m2n))
ppp <-Heatmap(mat, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
        name = "scaled_expr", cluster_columns = F,
        clustering_method_rows = 'ward.D2',
        row_split = 2,column_split = c(rep('pre-HM',24),rep('pre-NHM',130)),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,#border_gp = gpar(color = 'yellow'),
        #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
        show_column_names = F, #width = unit(8, "cm"),
        heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gns2mark), 
                                 labels = rownames(mat)[rownames(mat) %in% gns2mark], 
                                 labels_gp = gpar(fontsize = 6))) 
ppp

nhm_wt2 = cl_wt$ID2[cl_wt$Hypermutation.of.Recurrent.Tumor=='No' & cl_wt$TMZ.Treatment.Before.Recurrence==1]
m1rh = ger[upgns,names(ger) %in% hm_wt]; m1rn = ger[upgns,names(ger) %in% nhm_wt2]
m2rh = ger[dngns,names(ger) %in% hm_wt]; m2rn = ger[dngns,names(ger) %in% nhm_wt2]
mat2 = rbind(cbind(m1rh,m1rn), cbind(m2rh,m2rn))
pppr <-Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
              name = "scaled_expr", cluster_columns = F,
              clustering_method_rows = 'ward.D2',
              row_split = 2,column_split = c(rep('HM',ncol(m1rh)),rep('NHM',ncol(m1rn))),
              row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,#border_gp = gpar(color = 'yellow'),
              #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
              show_column_names = F, #width = unit(8, "cm"),
              heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat2) %in% gns2mark), 
                                 labels = rownames(mat)[rownames(mat2) %in% gns2mark], 
                                 labels_gp = gpar(fontsize = 6))) 
pppr

library(dendextend)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

ego_up_wt <- enrichGO(gene       = deg_wt$gene[deg_wt$type=='up' & deg_wt$type_rec=='up'],
                   universe      = deg_wt$gene,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego_up_wt),n=10)
dotplot(ego_up_wt, size = "GeneRatio", x='p.adjust',showCategory = 10,label_format = 40)+scale_x_log10()
cnetplot(ego_up_wt, categorySize="pvalue",circular = F, colorEdge = TRUE)

ego_dw_wt<- enrichGO(gene       = deg_wt$gene[deg_wt$type=='down' & deg_wt$type_rec=='down'],
                   universe      = deg_wt$gene,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego_dw_wt),n=10)
dotplot(ego_dw_wt, size = "GeneRatio", x='p.adjust',showCategory = 10,label_format = 40)+scale_x_log10()
cnetplot(ego_dw_wt, categorySize="pvalue",circular = F, colorEdge = TRUE)


d1 = row_dend(ppp)[[1]] #the down-regulated ones
cutree(d1,2) ->tc1 #cut into 2 clusters
genes1.1 = names(tc1[tc1==1])
ego1.1 <- enrichGO(gene       = genes1.1,
                universe      = rownames(gei),
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                keyType       = 'SYMBOL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1)
head(summary(ego1.1),n=10)
dotplot(ego1.1)
cnetplot(ego1.1, categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
#gn2label =c('ERCC2','DNMT1','POLD1','POLE','MCM5','E2F1','AKT1')
gn2label = c('LCA5','LCA5L','IFT57','TRAF3IP1','TTC30A')
genes1.2 = names(tc1[tc1==2])
ego1.2 <- enrichGO(gene       = genes1.2,
                universe      = rownames(gei),
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                keyType       = 'SYMBOL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1)
head(summary(ego1.2))
dotplot(ego1.2)
cnetplot(ego1.2, categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
#gn2label = c(gn2label, 'MKI67','PLK1','CENPE','UBE2C','KIF18B','AURKB')
gn2label = c(gn2label, 'PEX2','PEX12','ZFAND6','PEX7','MGMT')

d2 = row_dend(ppp)[[2]] #the up regulated ones
cutree(d2,1) ->tc2  #divide into one subclusters
genes2.1 = names(tc2[tc2==1])
ego2.1 <- enrichGO(gene          = genes2.1,
                universe      = rownames(gei),
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                keyType       = 'SYMBOL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1)
head(summary(ego2.1))
dotplot(ego2.1)
cnetplot(ego2.1, categorySize="pvalue",circular = F, colorEdge = TRUE)
#gn2label = c(gn2label,c('ANXA1','CD44','CD274','CCR2','IL10','IL7R','CXCL3'))
gn2label = c(gn2label,c('POLD1','POLE','FANCC','MCM5','PIF1','RAD54L','KIF18B','DNMT1','BRD1','BRD4',
                        'PTBP1','SAFB','KMT2B','MBD3','PRKD2','SNRPB','WIZ','SUGP2'))

# genes2.2 = names(tc2[tc2==2])
# ego2.2 <- enrichGO(gene          = genes2.2,
#                    universe      = rownames(gei),
#                    OrgDb         = org.Hs.eg.db,
#                    ont           = "ALL",
#                    keyType       = 'SYMBOL',
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.05,
#                    qvalueCutoff  = 0.1)
# head(summary(ego2.2)) 
# dotplot(ego2.2)
# cnetplot(ego2.2, categorySize="pvalue",circular = F, colorEdge = TRUE)
# gn2label = c(gn2label,c('CDKN1A','MID1','VCL'))

# genes2.3 = names(tc2[tc2==3])
# ego2.3 <- enrichGO(gene          = genes2.3,
#                    universe      = rownames(gei),
#                    OrgDb         = org.Hs.eg.db,
#                    ont           = "ALL",
#                    keyType       = 'SYMBOL',
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.05,
#                    qvalueCutoff  = 0.1)
# head(summary(ego2.3))
# dotplot(ego2.3)
# cnetplot(ego2.3, categorySize="pvalue",circular = F, colorEdge = TRUE)
# gn2label = c(gn2label,c('MT3','MGMT','CDKN1A','ALDH1A1'))

pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt_deg.heatmap0830.pdf',width =5.2,height = 3.2)
ppp <-Heatmap(mat, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
        name = "scaled_expr", cluster_columns = F,cluster_rows = rev(hclust(dist(mat),method = 'ward.D2')),
        #row_dend_reorder = F,
        #clustering_method_rows = 'ward.D',
        row_split = 2,column_split = c(rep('pre-HM',24),rep('pre-NHM',130)),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,#border_gp = gpar(color = 'yellow'),
        #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
        show_column_names = FALSE, #width = unit(8, "cm"),
        heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gn2label), 
                                 labels = rownames(mat)[rownames(mat) %in% gn2label], 
                                 labels_gp = gpar(fontsize = 5))) 
# ppp1 <-Heatmap(mat, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
#                name = "scaled_expr", cluster_columns = F,cluster_rows = rev(row_dend(ppp)),
#                #clustering_method_rows = 'ward.D',
#                #row_split = 2,
#                column_split = c(rep('pre-HM',24),rep('pre-NHM',130)),
#                row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,#border_gp = gpar(color = 'yellow'),
#                #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
#                show_column_names = FALSE, #width = unit(8, "cm"),
#                heatmap_legend_param = list(title = "Exp")) +
#   rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gn2label), 
#                                  labels = rownames(mat)[rownames(mat) %in% gn2label], 
#                                  labels_gp = gpar(fontsize = 5))) 
ppp
dev.off()

pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt_deg.recurrence.heatmap0830.pdf',width =5.2,height = 3.2)
pppr <-Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
               name = "scaled_expr", cluster_columns = F,#cluster_rows = rev(hclust(dist(mat2),method = 'ward.D2')),
               #row_dend_reorder = F,
               #clustering_method_rows = 'ward.D',
               row_split = 2,column_split = c(rep('HM',28),rep('NHM',137)),
               row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,#border_gp = gpar(color = 'yellow'),
               #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
               show_column_names = FALSE, #width = unit(8, "cm"),
               heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gn2label), 
                                 labels = rownames(mat)[rownames(mat) %in% gn2label], 
                                 labels_gp = gpar(fontsize = 5)))  
pppr
dev.off()

deg_wt$lab = ifelse(deg_wt$gene %in% gn2label & deg_wt$pval<0.05, deg_wt$gene, NA)
ggplot(deg_wt[!is.na(deg_wt$pval) & deg_wt$log2fc!=0,], aes(x = log2fc, y = -log10(pval), label = lab))+
  geom_point(aes( color = type), size = 0.2, alpha = 0.7, show.legend = F)+
  geom_vline(xintercept = 0, lty=2)+
  geom_text_repel(size=2.,segment.size = 0.25, max.iter = 1e4,max.overlaps = 30,position = pos)+ #color = "white",  bg.color = "grey30", colorbg.r = 0.15
  scale_color_manual(values = c('blue','#dddddd','red'))+
  theme_classic()+#lims(x = c(-2.,2))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'Expression difference', y = '-log10(P)')

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#IDHnon
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cl_non = cl[which(cl$Subtype=='IDHmut-noncodel'& cl$TMZ.Treatment.Before.Recurrence==1),] #
prehm_non = cl_non$ID1[cl_non$Hypermutation.of.Recurrent.Tumor=='Yes']
prenhm_non = cl_non$ID1[cl_non$Hypermutation.of.Recurrent.Tumor=='No']

hm_non = cl_non$ID2[cl_non$Hypermutation.of.Recurrent.Tumor=='Yes']
nhm_non = cl_non$ID2[cl_non$Hypermutation.of.Recurrent.Tumor=='No']

plotgeneinon <-function(gene = 'MGMT',offset = 2,mytest = 't.test'){ #examine a single gene
  x = as.numeric(gei[gene,names(gei)%in% prehm_non]) #Pre-HM
  y = as.numeric(gei[gene,names(gei)%in% prenhm_non]) #pre-NHM
  df = data.frame(var = c(rep('pre-HM', length(x)),rep('pre-NHM', length(y)) ), val = c(x,y))
  df$SID = c(names(gei)[names(gei)%in% prehm_non], names(gei)[names(gei)%in% prenhm_non])
  df$isnew = ifelse(df$SID %in% batch2$pairID, 'isnew','notnew')
  df$CaseID = substr(df$SID,1,5)
  df$Timepoint = substr(df$SID, 7,7)
  
  xr = as.numeric(ger[gene,names(ger)%in% hm_non]) #HM
  yr = as.numeric(ger[gene,names(ger)%in% nhm_non]) #NHM
  dfr = data.frame(var = c(rep('HM', length(xr)),rep('NHM', length(yr)) ), val = c(xr,yr))
  dfr$SID = c(names(ger)[names(ger)%in% hm_non], names(ger)[names(ger)%in% nhm_non])
  dfr$isnew = ifelse(dfr$SID %in% batch2$pairID, 'isnew','notnew')
  dfr$CaseID = substr(dfr$SID,1,5)
  dfr$Timepoint = substr(dfr$SID, 7,7)
  
  df = rbind(df, dfr)
  df$var = factor(df$var, levels = c('pre-HM','HM','pre-NHM','NHM'))
  df$hmstatus = ifelse(df$var %in% c('HM','pre-HM'),'H','N')
  p<- 
    ggplot(df, aes(x = var , y = val, color = hmstatus))+
    geom_boxplot(width = 0.4, show.legend = F, outlier.shape = NA)+
    geom_line(aes(group = CaseID), color = '#ddddddaa', size=1)+
    geom_jitter( show.legend = F, size=0.5,width = 0.1)+
    stat_compare_means(label = 'p.format', method = mytest, show.legend = F,
                       comparisons = list(c('pre-HM','HM'),c('NHM','pre-NHM'),c('NHM','HM'),c('pre-HM','pre-NHM')))+
    scale_color_manual(values = c('#e41a1c','#377eb8'))+
    theme_classic()+theme(axis.text = element_text(color = 'black'))+
    lims(y=c(NA, offset*max(df$val,na.rm = T)))+
    labs(x = '', y = paste0('Scaled ',gene, ' expression'), color = '',title = 'IDHmut-non-codel')
  print(p)
}
plotgeneinon(gene = 'MKI67')

plotgeneinon(gene = 'MGMT')
plotgeneinon(gene = 'TIPIN')
plotgeneinon(gene = 'RCBTB1',offset = 2)

#differential gene expression analysis
deg_non_ini = data.frame(gene = rownames(gei), med_hm = 0,med_nhm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_non_ini)){
  x = as.numeric(gei[i,names(gei)%in% prehm_non]) #Pre-HM
  y = as.numeric(gei[i,names(gei)%in% prenhm_non]) #pre-NHM
  deg_non_ini$med_hm[i] = mean(x)
  deg_non_ini$med_nhm[i] = mean(y)
  
  deg_non_ini$log2fc[i] = deg_non_ini$med_hm[i] - deg_non_ini$med_nhm[i]
  
  #deg_non_ini$pval[i] = t.test(x,y)$p.value
  #deg_non_ini$tstat[i] = t.test(x,y)$statistic
  
  deg_non_ini$pval[i] = wilcox.test(x,y)$p.value
  deg_non_ini$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_non_ini$type = ifelse(deg_non_ini$pval>0.05, 'nosig',ifelse(deg_non_ini$log2fc>0,'up','down'))

deg_non_rec = data.frame(gene = rownames(ger), med_hm = 0,med_nhm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_non_rec)){
  x = as.numeric(ger[i,names(ger)%in% hm_non]) #HM
  y = as.numeric(ger[i,names(ger)%in% nhm_non]) #NHM
  deg_non_rec$med_hm[i] = mean(x)
  deg_non_rec$med_nhm[i] = mean(y)
  
  deg_non_rec$log2fc[i] = deg_non_rec$med_hm[i] - deg_non_rec$med_nhm[i]
  #deg_non_rec$pval[i] = t.test(x,y)$p.value
  #deg_non_rec$tstat[i] = t.test(x,y)$statistic
  
  deg_non_rec$pval[i] = wilcox.test(x,y)$p.value
  deg_non_rec$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_non_rec$type = ifelse(deg_non_rec$pval>0.05, 'nosig',ifelse(deg_non_rec$log2fc>0,'up','down'))

deg_non_hmpair = data.frame(gene = rownames(gei), med_prehm = 0,med_hm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_non_hmpair)){
  x = as.numeric(gei[i,names(gei)%in% prehm_non]) #pre-HM
  y = as.numeric(ger[i,names(ger)%in% hm_non]) #HM
  deg_non_hmpair$med_prehm[i] = mean(x)
  deg_non_hmpair$med_hm[i] = mean(y)
  
  deg_non_hmpair$log2fc[i] = deg_non_hmpair$med_hm[i] - deg_non_hmpair$med_prehm[i]
  #deg_non_hmpair$pval[i] = t.test(x,y)$p.value
  #deg_non_hmpair$tstat[i] = t.test(x,y)$statistic
  
  deg_non_hmpair$pval[i] = wilcox.test(x,y)$p.value
  deg_non_hmpair$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_non_hmpair$type = ifelse(deg_non_hmpair$pval>0.05, 'nosig',ifelse(deg_non_hmpair$log2fc>0,'up','down'))

deg_non_ini$timepoint = 'initial'; deg_non_rec$timepoint = 'recurrence';deg_non_hmpair$timepoint = 'pair'

deg_non = deg_non_ini[,c('gene','med_hm','med_nhm','log2fc','pval','tstat','type')]
names(deg_non) = paste0(names(deg_non),"_ini")
deg_non$log2fc_rec = deg_non_rec$log2fc
deg_non$pval_rec = deg_non_rec$pval
deg_non$type_rec = deg_non_rec$type
deg_non$pval_pair = deg_non_hmpair$pval
deg_non$type_pair = deg_non_hmpair$type

write.table(deg_non, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel.deg.ranksum.txt', row.names = F, quote = F, sep = "\t")

ggplot(deg_non, aes(x = log2fc_ini, y = log2fc_rec))+
  #geom_point( alpha = 0.1)+
  geom_hex(binwidth = c(.03, 0.03),show.legend = T)+
  #stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200)+
  geom_hline(yintercept = 0,lty=2, col = '#999999')+
  geom_vline(xintercept = 0, lty = 2, col = '#999999')+
  scale_fill_viridis_c()+
  #scale_fill_continuous(low = "white", high = "dodgerblue4")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x = 'pre-HM versus pre-NHM', y = 'HM versus NHM', fill = '# genes')

ggplot(deg_non, aes(x = sign(log2fc_ini)*-log10(pval_ini), y = sign(log2fc_rec)*-log10(pval_rec)))+
  #geom_point( alpha = 0.1)+
  geom_hex(binwidth = c(.25, 0.25),show.legend = T)+
  #stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200)+
  geom_hline(yintercept = 0,lty=2, col = '#999999')+
  geom_vline(xintercept = 0, lty = 2, col = '#999999')+
  scale_fill_viridis_c()+
  #scale_fill_continuous(low = "white", high = "dodgerblue4")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x = 'pre-HM versus pre-NHM', y = 'HM versus NHM', fill = '# genes')


cor.test(deg_non$log2fc_ini, deg_non$log2fc_rec,alternative = 'two.sided')$p.value
pheatmap(table(deg_non$type_ini, deg_non$type_rec),display_numbers = T,number_format = '%.0f', legend = F, cluster_rows = F, cluster_cols = F)

#fgsea
deg_non_ini = deg_non_ini[order(deg_non_ini$tstat, decreasing = T),]

fcrnk = deg_non_ini$tstat
names(fcrnk) = deg_non_ini$gene  
fcrnk = fcrnk[!is.na(fcrnk)]
fgseaRes <- fgsea(pathways = hmpathways, stats = fcrnk, minSize=15,maxSize=500)

plotfgsea(pathway = hmpathways[["HALLMARK_MYC_TARGETS_V2"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MYC_TARGETS_V2',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')
plotfgsea(pathway = hmpathways[["HALLMARK_MITOTIC_SPINDLE"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MITOTIC_SPINDLE',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')

#
fgseaRes = fgseaRes[order(-fgseaRes$NES, decreasing = T),]
fgseaRes$pathway = factor(fgseaRes$pathway, levels = fgseaRes$pathway)
fgseaRes$pathway1 = gsub("HALLMARK_","", fgseaRes$pathway)
fgseaRes$pathway1  = factor(fgseaRes$pathway1, levels = fgseaRes$pathway1)
fgseaRes$pval[fgseaRes$pval< 1e-5] = 1e-5
ggplot(aes(x = pathway1, y = NES, fill = sign(NES)*-log10(pval)),data = fgseaRes[fgseaRes$pval<0.05,])+
  geom_bar(stat = 'identity', color = 'black') + 
  scale_fill_gradient2(low = 'blue',high='red')+
  coord_flip() + 
  theme_minimal()+labs(fill = 'Significance', x= '')+
  theme(axis.text.x = element_text(colour = "black"),#panel.grid = element_blank(),
        axis.text.y = element_text(colour = "black"))

gns2mark = deg_non$gene[deg_non$gene %in% c(hmpathways$HALLMARK_MYC_TARGETS_V1,hmpathways$HALLMARK_MYC_TARGETS_V2,
                                            hmpathways$HALLMARK_E2F_TARGETS,hmpathways$HALLMARK_G2M_CHECKPOINT,hmpathways$HALLMARK_MITOTIC_SPINDLE)];
gns2mark = c(gns2mark,'MGMT')
deg_non$lab = ifelse(deg_non$gene %in% gns2mark,deg_non$gene, NA)

deg_non_predictors_up = deg_non[deg_non$type_ini=='up' & deg_non$type_rec=='up' & deg_non$type_pair!='down',]
deg_non_predictors_dn = deg_non[deg_non$type_ini=='down' & deg_non$type_rec=='down' & deg_non$type_pair!='up',]
deg_non_predictors = rbind(deg_non_predictors_up, deg_non_predictors_dn)
deg_non_predictors$pcomb = 1
for (i in 1:nrow(deg_non_predictors)){
  deg_non_predictors$pcomb[i] = pchisq( -2*sum(log(c(deg_non_predictors$pval_ini[i],deg_non_predictors$pval_rec[i]))), df=4, lower.tail=FALSE)
}

deg_non_predictors$signedpcomb = sign(deg_non_predictors$log2fc_ini)*-log10(deg_non_predictors$pcomb)
deg_non_predictors = deg_non_predictors[order(deg_non_predictors$signedpcomb, decreasing = F),]
deg_non_predictors$idx = 1:nrow(deg_non_predictors)
deg_non_predictors$myctagrtes = ifelse(deg_non_predictors$gene_ini %in% c(hmpathways$HALLMARK_MYC_TARGETS_V1,hmpathways$HALLMARK_MYC_TARGETS_V2),deg_non_predictors$gene_ini,NA)
ggplot()+
  geom_bar(data = deg_non_predictors, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_non_predictors[!is.na(deg_non_predictors$myctagrtes),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_non_predictors[!is.na(deg_non_predictors$myctagrtes),],
                  size=3,segment.size=0.25,max.iter =1000,max.overlaps = 100,
                  mapping = aes(x = idx, y = signedpcomb,label = lab))+
  theme_classic()+#coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')


ggplot()+
  geom_bar(data = deg_non_predictors, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_non_predictors[!is.na(deg_non_predictors$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_non_predictors[!is.na(deg_non_predictors$lab),],
                  size=3,segment.size=0.25,max.iter =1000,max.overlaps = 30,
                  mapping = aes(x = idx, y = signedpcomb,label = lab))+
  theme_classic()+#coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')

deg_non_predictors$comb.padj = p.adjust(deg_non_predictors$pcomb)
deg_non_predictors_pcrt = deg_non_predictors[deg_non_predictors$comb.padj<0.05,]
deg_non_predictors_pcrt$idx = 1:nrow(deg_non_predictors_pcrt)
ggplot()+
  geom_bar(data = deg_non_predictors_pcrt, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_non_predictors_pcrt,#[!is.na(deg_non_predictors_pcrt$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_non_predictors_pcrt,#[!is.na(deg_non_predictors_pcrt$lab),],
                  size=2.5,segment.size=0.1,max.iter =1000,max.overlaps = 10,#nudge_x = -1,
                  hjust=1,direction= "y",
                  nudge_y = ifelse(deg_non_predictors_pcrt$log2fc_rec>0,5,-5),
                  mapping = aes(x = idx, y = signedpcomb,label = gene_ini))+
  theme_classic()+coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')

deg_non_predictors_pcrt$lab[deg_non_predictors_pcrt$gene_ini=='RCBTB1'] = 'RCBTB1'
deg_non_predictors_pcrt$lab[deg_non_predictors_pcrt$gene_ini=='CLDN6'] = 'CLDN6'

ggplot()+
  geom_bar(data = deg_non_predictors_pcrt, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_non_predictors_pcrt[!is.na(deg_non_predictors_pcrt$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_non_predictors_pcrt[!is.na(deg_non_predictors_pcrt$lab),],
                  size=2.5,segment.size=0.1,max.iter =1000,max.overlaps = 10,#nudge_x = -1,
                  hjust=1,direction= "y",
                  nudge_y = ifelse(deg_non_predictors_pcrt$log2fc_rec[!is.na(deg_non_predictors_pcrt$lab)]>0,5,-5),
                  mapping = aes(x = idx, y = signedpcomb,label = gene_ini))+
  theme_classic()+coord_flip()+lims(y=c(-7,5))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  labs(x = 'Genes', y = 'signed log10(P)', fill = '')


ego_up_non <- enrichGO(gene       = deg_non$gene[deg_non$type=='up' & deg_non$type_rec=='up'],
                   universe      = deg_non$gene,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego_up_non),n=10)
dotplot(ego_up_non, size = "GeneRatio", x='p.adjust',showCategory = 10,label_format = 40)+scale_x_log10()
cnetplot(ego_up_non, categorySize="pvalue",circular = F, colorEdge = TRUE)

ego_dw_non<- enrichGO(gene       = deg_non$gene[deg_non$type=='down' & deg_non$type_rec=='down'],
                  universe      = deg_non$gene,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  keyType       = 'SYMBOL',
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1)
head(summary(ego_dw_non),n=10)
dotplot(ego_dw_non, size = "GeneRatio", x='p.adjust',showCategory = 10,label_format = 40)+scale_x_log10()
cnetplot(ego_dw_non, categorySize="pvalue",circular = F, colorEdge = TRUE)


#fgsea
deg_non_ini = deg_non_ini[order(deg_non_ini$tstat, decreasing = T),]

fcrnk = deg_non_ini$tstat
names(fcrnk) = deg_non_ini$gene  
fgseaRes <- fgsea(pathways = hmpathways, stats = fcrnk, minSize=15,maxSize=500)

plotfgsea(pathway = hmpathways[["HALLMARK_MYC_TARGETS_V2"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MYC_TARGETS_V2',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')
plotfgsea(pathway = hmpathways[["HALLMARK_E2F_TARGETS"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_E2F_TARGETS',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')

#
fgseaRes = fgseaRes[order(-fgseaRes$NES, decreasing = T),]
fgseaRes$pathway = factor(fgseaRes$pathway, levels = fgseaRes$pathway)
fgseaRes$pathway1 = gsub("HALLMARK_","", fgseaRes$pathway)
fgseaRes$pathway1  = factor(fgseaRes$pathway1, levels = fgseaRes$pathway1)
ggplot(aes(x = pathway1, y = NES, fill = NES),data = fgseaRes[fgseaRes$padj<0.05,])+
  geom_bar(stat = 'identity', color = 'black') + 
  scale_fill_gradient2(low = 'blue',high = 'red')+
  coord_flip() +
  theme_classic()+labs(x='')+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black", size = 6))


gns2mark = deg_non$gene[deg_non$type=='up'  &(deg_non$gene %in% c(hmpathways$HALLMARK_MYC_TARGETS_V1,
                                                                                       hmpathways$HALLMARK_MYC_TARGETS_V2,
                                                                                       hmpathways$HALLMARK_E2F_TARGETS,
                                                                                       hmpathways$HALLMARK_G2M_CHECKPOINT,
                                                                                       hmpathways$HALLMARK_MITOTIC_SPINDLE))];
gns2mark = c(gns2mark,'MGMT')
deg_non$lab = ifelse(deg_non$gene %in% gns2mark,deg_non$gene, NA)
library(ggbeeswarm)
pos <- position_quasirandom()
ggplot(deg_non[!is.na(deg_non$pval) & deg_non$log2fc!=0,], aes(x = log2fc, y = -log10(pval), label = lab))+
  geom_point(aes( color = type), size = 0.2, alpha = 0.7, show.legend = F)+
  geom_vline(xintercept = 0, lty=2)+
  geom_text_repel(size=2.,segment.size = 0.25, max.iter = 1e4,max.overlaps = 40,position = pos)+ #color = "white",  bg.color = "grey30", colorbg.r = 0.15
  scale_color_manual(values = c('blue','#dddddd','red'))+
  theme_classic()+#lims(x = c(-2.,2))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'log2 (pre-HM/pre-NHM)', y = '-log10(P)')

#MGMT methylation and expression and hypermutation
mgmt = read.delim('~/Dropbox/communter//Rev1/MYCbySubtype/MYC_methylation_subtype.txt', 
                  stringsAsFactors = F,na.strings = c("NA","#N/A","NaN"), row.names = 1)
x = as.numeric(gei['MGMT',names(gei)%in% prehm_non]) #Pre-HM
y = as.numeric(gei['MGMT',names(gei)%in% prenhm_non]) #pre-NHM
df = data.frame(var = c(rep('pre-HM', length(x)),rep('pre-NHM', length(y)) ), val = c(x,y))
df$SID = c(names(gei)[names(gei) %in%prehm_non ], names(gei)[names(gei) %in%prenhm_non ])
df$MGMTmethy = mgmt$MGMT_methylation_I[match(substr(df$SID,1,5),rownames(mgmt))]
df$MGMTmethy=ifelse(is.na(df$MGMTmethy),'unknown',ifelse(df$MGMTmethy==1,'methylated','unmethylated'))
df$TMZ = cl_non$TMZ.Treatment.Before.Recurrence[match(substr(df$SID,1,5), cl_non$Patient.ID)]
df$MGMTmethy = factor(df$MGMTmethy, levels = c('methylated','unmethylated','unknown'))
df_non_mgmt = df
ggplot(df[which(df$TMZ==1 ),], aes(x = MGMTmethy , y = val, color = var))+
  geom_boxplot(width = 0.5, show.legend = T, outlier.shape = NA)+
  geom_jitter( show.legend = T,position = position_dodge(width = 0.5))+
  scale_shape_manual(values = c(17,19,4))+
  stat_compare_means(label = 'p.format',show.legend = F,label.y.npc = 0.95)+
  scale_color_manual(values = c('#e41a1c','#377eb8','#cccccc'), guide = 'none')+
  theme_classic()+theme(axis.text = element_text(color = 'black'))+
  lims(y=c(NA, 1.8*max(df$val,na.rm = T)))+
  #facet_wrap(~MGMTmethy)+
  #scale_y_continuous(expand = c(0.1,1))+
  labs(x = '', y =  'Scaled MGMT expression', color = '')


# upgns = NULL
# for (pwy in fgseaRes$pathway[which(fgseaRes$padj<0.05 &fgseaRes$NES >0)]){
#   upgns = append(upgns, hmpathways[[as.character(pwy)]])
# }
# upgns = unique(upgns)
# upgns = upgns[upgns %in% deg_non$gene[deg_non$pval<0.05]]
# 
# dngns = NULL
# for (pwy in fgseaRes$pathway[which(fgseaRes$padj<0.05 &fgseaRes$NES <0)]){
#   dngns = append(dngns, hmpathways[[as.character(pwy)]])
# }
# dngns = unique(dngns)
# dngns = dngns[dngns %in% deg_non$gene[deg_non$pval<0.05]]

upgns = deg_non$gene[which(deg_non$type=='up' )]
dngns = deg_non$gene[which(deg_non$type=='down' )]

prenhm_non2 = cl_non$ID1[cl_non$Hypermutation.of.Recurrent.Tumor=='No' & cl_non$TMZ.Treatment.Before.Recurrence==1]
m1h = gei[upgns,names(gei) %in% prehm_non]
m1n = gei[upgns,names(gei) %in% prenhm_non2]

m2h = gei[dngns,names(gei) %in% prehm_non]
m2n = gei[dngns,names(gei) %in% prenhm_non2]

mat = rbind(cbind(m1h,m1n), cbind(m2h,m2n))
#mat = t(scale(t(log2(mat+1))))
gns2mark = c(upgns[upgns %in% hmpathways$HALLMARK_MYC_TARGETS_V1], dngns)
#pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHnon_deg.heatmap.pdf',width =4,height = 2.5)
ppp<- Heatmap(mat, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
        name = "scaled_expr", cluster_columns = F,
        clustering_method_rows = 'ward.D',
        row_split = 2,column_split = c(rep('pre-HM',14),rep('pre-NHM',35)),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
        #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
        show_column_names = FALSE, #width = unit(8, "cm"),
        heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gns2mark), 
                                 labels = rownames(mat)[rownames(mat) %in% gns2mark], 
                                 labels_gp = gpar(fontsize = 6))) 
ppp

nhm_non2 = cl_non$ID2[cl_non$Hypermutation.of.Recurrent.Tumor=='No' & cl_non$TMZ.Treatment.Before.Recurrence==1]
m1rh = ger[upgns,names(ger) %in% hm_non]; m1rn = ger[upgns,names(ger) %in% nhm_non2]
m2rh = ger[dngns,names(ger) %in% hm_non]; m2rn = ger[dngns,names(ger) %in% nhm_non2]
mat2 = rbind(cbind(m1rh,m1rn), cbind(m2rh,m2rn))

pppr<- Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
              name = "scaled_expr", cluster_columns = F,
              clustering_method_rows = 'ward.D',
              row_split = 2,column_split = c(rep('HM',14),rep('NHM',40)),
              row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
              #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
              show_column_names = FALSE, #width = unit(8, "cm"),
              heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gns2mark), 
                                 labels = rownames(mat)[rownames(mat) %in% gns2mark], 
                                 labels_gp = gpar(fontsize = 6))) 
pppr

#dev.off()
d1 = row_dend(ppp)[[1]] #the first large cluster
# cutree(d1,1) ->tc1 #as a whole
# genes1.1 = names(tc1[tc1==3])
# ego1.1 <- enrichGO(gene       = genes1.1,
#                    universe      = rownames(gei),
#                    OrgDb         = org.Hs.eg.db,
#                    ont           = "ALL",
#                    keyType       = 'SYMBOL',
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.05,
#                    qvalueCutoff  = 0.1)
# head(summary(ego1.1),n=10)
# dotplot(ego1.1)
# cnetplot(ego1.1, categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
#weird, will turn to MetaScape
gn2label = c("CDKN2D","RAB8A","MNAT1","MSH4","PRKACA","RFC2","SEC13","UBE2I","RBX1","GMNN","HAUS8","CENPS")



d2 = row_dend(ppp)[[2]] #the 2nd large  cluster
cutree(d2,2) ->tc2  #separate into two subclusters
genes2.1 = names(tc2[tc2==1])
ego2.1 <- enrichGO(gene          = genes2.1,
                   universe      = rownames(gei),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
summary(ego2.1) ->xxx
dotplot(ego2.1)
cnetplot(ego2.1, categorySize="pvalue",circular = F, colorEdge = TRUE)
gn2label = c(gn2label,c('FBLN2','OGN','FBN1','DCN'))

genes2.2 = names(tc2[tc2==2])
ego2.2 <- enrichGO(gene          = genes2.2,
                   universe      = rownames(gei),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego2.2))
dotplot(ego2.2)
cnetplot(ego2.2, categorySize="pvalue",circular = F, colorEdge = TRUE)
 
gn2label = c(gn2label,c('PTPRS','PTPRD','KLHL24','SYT11'))
pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel_deg.heatmap0830.pdf',width =3.5,height = 2.5)
Heatmap(mat, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
        name = "scaled_expr", cluster_columns = F,
        clustering_method_rows = 'ward.D',
        row_split = 2,column_split = c(rep('pre-HM',14),rep('pre-NHM',35)),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,#border_gp = gpar(color = 'yellow'),
        #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
        show_column_names = FALSE, #width = unit(8, "cm"),
        heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gn2label), 
                                 labels = rownames(mat)[rownames(mat) %in% gn2label], 
                                 labels_gp = gpar(fontsize = 5))) 
dev.off()


pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel_deg.recurrence.heatmap0902.pdf',width =3.5,height = 2.5)
pppr<- Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
               name = "scaled_expr", cluster_columns = F,
               clustering_method_rows = 'ward.D',
               row_split = 2,column_split = c(rep('HM',14),rep('NHM',40)),
               row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
               #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
               show_column_names = FALSE, #width = unit(8, "cm"),
               heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gn2label), 
                                 labels = rownames(mat)[rownames(mat) %in% gn2label], 
                                 labels_gp = gpar(fontsize = 6))) 
pppr
dev.off()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#IDHcocel
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cl_cod = cl[which((cl$Subtype=='IDHmut-codel') & cl$TMZ.Treatment.Before.Recurrence==1),]
prehm_cod = cl_cod$ID1[cl_cod$Hypermutation.of.Recurrent.Tumor=='Yes']
prenhm_cod = cl_cod$ID1[cl_cod$Hypermutation.of.Recurrent.Tumor=='No']

hm_cod = cl_cod$ID2[cl_cod$Hypermutation.of.Recurrent.Tumor=='Yes']
nhm_cod = cl_cod$ID2[cl_cod$Hypermutation.of.Recurrent.Tumor=='No']

plotgeneicod <-function(gene = 'MGMT',offset = 2,mytest = 't.test'){ #examine a single gene
  x = as.numeric(gei[gene,names(gei)%in% prehm_cod]) #Pre-HM
  y = as.numeric(gei[gene,names(gei)%in% prenhm_cod]) #pre-NHM
  df = data.frame(var = c(rep('pre-HM', length(x)),rep('pre-NHM', length(y)) ), val = c(x,y))
  df$SID = c(names(gei)[names(gei)%in% prehm_cod], names(gei)[names(gei)%in% prenhm_cod])
  df$isnew = ifelse(df$SID %in% batch2$pairID, 'isnew','notnew')
  df$CaseID = substr(df$SID,1,5)
  df$Timepoint = substr(df$SID, 7,7)
  
  xr = as.numeric(ger[gene,names(ger)%in% hm_cod]) #HM
  yr = as.numeric(ger[gene,names(ger)%in% nhm_cod]) #NHM
  dfr = data.frame(var = c(rep('HM', length(xr)),rep('NHM', length(yr)) ), val = c(xr,yr))
  dfr$SID = c(names(ger)[names(ger)%in% hm_cod], names(ger)[names(ger)%in% nhm_cod])
  dfr$isnew = ifelse(dfr$SID %in% batch2$pairID, 'isnew','notnew')
  dfr$CaseID = substr(dfr$SID,1,5)
  dfr$Timepoint = substr(dfr$SID, 7,7)
  
  df = rbind(df, dfr)
  df$var = factor(df$var, levels = c('pre-HM','HM','pre-NHM','NHM'))
  df$hmstatus = ifelse(df$var %in% c('HM','pre-HM'),'H','N')
  p<- 
    ggplot(df, aes(x = var , y = val, color = hmstatus))+
    geom_boxplot(width = 0.4, show.legend = F, outlier.shape = NA)+
    geom_line(aes(group = CaseID), color = '#ddddddaa', size=1)+
    geom_jitter( show.legend = F, size=0.5,width = 0.1)+
    stat_compare_means(label = 'p.format', method = mytest, show.legend = F,
                       comparisons = list(c('pre-HM','HM'),c('NHM','pre-NHM'),c('NHM','HM'),c('pre-HM','pre-NHM')))+
    scale_color_manual(values = c('#e41a1c','#377eb8'))+
    theme_classic()+theme(axis.text = element_text(color = 'black'))+
    lims(y=c(NA, offset*max(df$val,na.rm = T)))+
    labs(x = '', y = paste0('Scaled ',gene, ' expression'), color = '',title = 'IDHmut-codel')
  print(p)
}
plotgeneicod(gene = 'MKI67',mytest = 'wilcox.test')
plotgeneicod(gene = 'PLK1',mytest = 'wilcox.test')
plotgeneicod(gene = 'MGMT',mytest = 'wilcox.test')
plotgeneicod(gene = 'PRC1',mytest = 't.test')


#differential gene expression analysis
deg_cod_ini = data.frame(gene = rownames(gei), med_hm = 0,med_nhm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_cod_ini)){
  x = as.numeric(gei[i,names(gei)%in% prehm_cod]) #Pre-HM
  y = as.numeric(gei[i,names(gei)%in% prenhm_cod]) #pre-NHM
  deg_cod_ini$med_hm[i] = mean(x)
  deg_cod_ini$med_nhm[i] = mean(y)
  
  deg_cod_ini$log2fc[i] = deg_cod_ini$med_hm[i] - deg_cod_ini$med_nhm[i]
  
  #deg_cod_ini$pval[i] = t.test(x,y)$p.value
  #deg_cod_ini$tstat[i] = t.test(x,y)$statistic
  
  deg_cod_ini$pval[i] = wilcox.test(x,y)$p.value
  deg_cod_ini$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_cod_ini$type = ifelse(deg_cod_ini$pval>0.05, 'nosig',ifelse(deg_cod_ini$log2fc>0,'up','down'))

deg_cod_rec = data.frame(gene = rownames(ger), med_hm = 0,med_nhm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_cod_rec)){
  x = as.numeric(ger[i,names(ger)%in% hm_cod]) #HM
  y = as.numeric(ger[i,names(ger)%in% nhm_cod]) #NHM
  deg_cod_rec$med_hm[i] = mean(x)
  deg_cod_rec$med_nhm[i] = mean(y)
  
  deg_cod_rec$log2fc[i] = deg_cod_rec$med_hm[i] - deg_cod_rec$med_nhm[i]
  #deg_cod_rec$pval[i] = t.test(x,y)$p.value
  #deg_cod_rec$tstat[i] = t.test(x,y)$statistic
  
  deg_cod_rec$pval[i] = wilcox.test(x,y)$p.value
  deg_cod_rec$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_cod_rec$type = ifelse(deg_cod_rec$pval>0.05, 'nosig',ifelse(deg_cod_rec$log2fc>0,'up','down'))

deg_cod_hmpair = data.frame(gene = rownames(gei), med_prehm = 0,med_hm=0, log2fc = 0,pval = 0,tstat = 0 )
for (i in 1:nrow(deg_cod_hmpair)){
  x = as.numeric(gei[i,names(gei)%in% prehm_cod]) #pre-HM
  y = as.numeric(ger[i,names(ger)%in% hm_cod]) #HM
  deg_cod_hmpair$med_prehm[i] = mean(x)
  deg_cod_hmpair$med_hm[i] = mean(y)
  
  deg_cod_hmpair$log2fc[i] = deg_cod_hmpair$med_hm[i] - deg_cod_hmpair$med_prehm[i]
  #deg_cod_hmpair$pval[i] = t.test(x,y)$p.value
  #deg_cod_hmpair$tstat[i] = t.test(x,y)$statistic
  
  deg_cod_hmpair$pval[i] = wilcox.test(x,y)$p.value
  deg_cod_hmpair$tstat[i] = wilcox.test(x,y)$statistic
  
}
deg_cod_hmpair$type = ifelse(deg_cod_hmpair$pval>0.05, 'nosig',ifelse(deg_cod_hmpair$log2fc>0,'up','down'))

deg_cod_ini$timepoint = 'initial'; deg_cod_rec$timepoint = 'recurrence';deg_cod_hmpair$timepoint = 'pair'

deg_cod = deg_cod_ini[,c('gene','med_hm','med_nhm','log2fc','pval','tstat','type')]
names(deg_cod) = paste0(names(deg_cod),"_ini")
deg_cod$log2fc_rec = deg_cod_rec$log2fc
deg_cod$pval_rec = deg_cod_rec$pval
deg_cod$type_rec = deg_cod_rec$type
deg_cod$pval_pair = deg_cod_hmpair$pval
deg_cod$type_pair = deg_cod_hmpair$type

write.table(deg_cod, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel.deg.ranksum.txt', row.names = F, quote = F, sep = "\t")

ggplot(deg_cod, aes(x = log2fc_ini, y = log2fc_rec))+
  #geom_point( alpha = 0.1)+
  geom_hex(binwidth = c(.15, 0.15),show.legend = T)+
  #stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200)+
  geom_hline(yintercept = 0,lty=2, col = '#999999')+
  geom_vline(xintercept = 0, lty = 2, col = '#999999')+
  scale_fill_viridis_c()+
  #scale_fill_continuous(low = "white", high = "dodgerblue4")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x = 'pre-HM versus pre-NHM', y = 'HM versus NHM', fill = '# genes')

ggplot(deg_cod, aes(x = sign(log2fc_ini)*-log10(pval_ini), y = sign(log2fc_rec)*-log10(pval_rec)))+
  #geom_point( alpha = 0.1)+
  geom_hex(binwidth = c(.05, 0.05),show.legend = T)+
  #stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200)+
  geom_hline(yintercept = 0,lty=2, col = '#999999')+
  geom_vline(xintercept = 0, lty = 2, col = '#999999')+
  scale_fill_viridis_c()+
  #scale_fill_continuous(low = "white", high = "dodgerblue4")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x = 'pre-HM versus pre-NHM', y = 'HM versus NHM', fill = '# genes')


cor.test(deg_cod$log2fc_ini, deg_cod$log2fc_rec,alternative = 'two.sided')$p.value
pheatmap(table(deg_cod$type_ini, deg_cod$type_rec),display_numbers = T,number_format = '%.0f', legend = F, cluster_rows = F, cluster_cols = F)

#fgsea
deg_cod_ini = deg_cod_ini[order(deg_cod_ini$tstat, decreasing = T),]

fcrnk = deg_cod_ini$tstat
names(fcrnk) = deg_cod_ini$gene  
fcrnk = fcrnk[!is.na(fcrnk)]
fgseaRes <- fgsea(pathways = hmpathways, stats = fcrnk, minSize=15,maxSize=500)

plotfgsea(pathway = hmpathways[["HALLMARK_MYC_TARGETS_V2"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MYC_TARGETS_V2',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')
plotfgsea(pathway = hmpathways[["HALLMARK_MITOTIC_SPINDLE"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MITOTIC_SPINDLE',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')

#
fgseaRes = fgseaRes[order(-fgseaRes$NES, decreasing = T),]
fgseaRes$pathway = factor(fgseaRes$pathway, levels = fgseaRes$pathway)
fgseaRes$pathway1 = gsub("HALLMARK_","", fgseaRes$pathway)
fgseaRes$pathway1  = factor(fgseaRes$pathway1, levels = fgseaRes$pathway1)
fgseaRes$pval[fgseaRes$pval< 1e-5] = 1e-5
ggplot(aes(x = pathway1, y = NES, fill = sign(NES)*-log10(pval)),data = fgseaRes[fgseaRes$pval<0.05,])+
  geom_bar(stat = 'identity', color = 'black') + 
  scale_fill_gradient2(low = 'blue',high='red')+
  coord_flip() + 
  theme_minimal()+labs(fill = 'Significance', x= '')+
  theme(axis.text.x = element_text(colour = "black"),#panel.grid = element_blank(),
        axis.text.y = element_text(colour = "black"))

gns2mark = deg_cod$gene[deg_cod$gene %in% c(hmpathways$HALLMARK_MYC_TARGETS_V1,hmpathways$HALLMARK_MYC_TARGETS_V2,
                                            hmpathways$HALLMARK_E2F_TARGETS,hmpathways$HALLMARK_G2M_CHECKPOINT,hmpathways$HALLMARK_MITOTIC_SPINDLE)];
gns2mark = c(gns2mark,'MGMT')
deg_cod$lab = ifelse(deg_cod$gene %in% gns2mark & deg_cod$log2fc_ini>0,deg_cod$gene, NA)

deg_cod_predictors_up = deg_cod[deg_cod$type_ini=='up' & deg_cod$type_rec=='up' & deg_cod$type_pair!='down',]
deg_cod_predictors_dn = deg_cod[deg_cod$type_ini=='down' & deg_cod$type_rec=='down' & deg_cod$type_pair!='up',]
deg_cod_predictors = rbind(deg_cod_predictors_up, deg_cod_predictors_dn)
deg_cod_predictors$pcomb = 1
for (i in 1:nrow(deg_cod_predictors)){
  deg_cod_predictors$pcomb[i] = pchisq( -2*sum(log(c(deg_cod_predictors$pval_ini[i],deg_cod_predictors$pval_rec[i]))), df=4, lower.tail=FALSE)
}

deg_cod_predictors$signedpcomb = sign(deg_cod_predictors$log2fc_ini)*-log10(deg_cod_predictors$pcomb)
deg_cod_predictors = deg_cod_predictors[order(deg_cod_predictors$signedpcomb, decreasing = F),]
deg_cod_predictors$idx = 1:nrow(deg_cod_predictors)
deg_cod_predictors$myctagrtes = ifelse(deg_cod_predictors$gene_ini %in% c(hmpathways$HALLMARK_MYC_TARGETS_V1,hmpathways$HALLMARK_MYC_TARGETS_V2),deg_cod_predictors$gene_ini,NA)
ggplot()+
  geom_bar(data = deg_cod_predictors, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_cod_predictors[!is.na(deg_cod_predictors$myctagrtes),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_cod_predictors[!is.na(deg_cod_predictors$myctagrtes),],
                  size=3,segment.size=0.25,max.iter =1000,max.overlaps = 100,
                  mapping = aes(x = idx, y = signedpcomb,label = lab))+
  theme_classic()+#coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')


ggplot()+
  geom_bar(data = deg_cod_predictors, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_cod_predictors[!is.na(deg_cod_predictors$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_cod_predictors[!is.na(deg_cod_predictors$lab),],
                  size=3,segment.size=0.25,max.iter =1000,max.overlaps = 30,
                  mapping = aes(x = idx, y = signedpcomb,label = lab))+
  theme_classic()+#coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')


deg_cod_predictors$comb.padj = p.adjust(deg_cod_predictors$pcomb)
deg_cod_predictors_pcrt = deg_cod_predictors[deg_cod_predictors$comb.padj<0.05,]
deg_cod_predictors_pcrt$idx = 1:nrow(deg_cod_predictors_pcrt)
ggplot()+
  geom_bar(data = deg_cod_predictors_pcrt, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_cod_predictors_pcrt,#[!is.na(deg_cod_predictors_pcrt$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_cod_predictors_pcrt,#[!is.na(deg_cod_predictors_pcrt$lab),],
                  size=2.5,segment.size=0.1,max.iter =1000,max.overlaps = 10,#nudge_x = -1,
                  hjust=1,direction= "y",
                  nudge_y = ifelse(deg_cod_predictors_pcrt$log2fc_rec>0,5,-15),
                  mapping = aes(x = idx, y = signedpcomb,label = gene_ini))+
  theme_classic()+coord_flip()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(x = 'Rank of genes', y = 'signed log10(P)', fill = '')

deg_cod_predictors_pcrt$lab[deg_cod_predictors_pcrt$gene_ini == 'UPF3A'] ='UPF3A'
deg_cod_predictors_pcrt$lab[deg_cod_predictors_pcrt$gene_ini == 'CYP3A4'] ='CYP3A4'
deg_cod_predictors_pcrt$lab[deg_cod_predictors_pcrt$gene_ini == 'ANXA9'] ='ANXA9'
ggplot()+
  geom_bar(data = deg_cod_predictors_pcrt, stat = 'identity',width = 0.75,show.legend = F,
           mapping=aes(x = idx, y = signedpcomb,fill = signedpcomb>0))+
  scale_fill_manual(values = c('blue','red'))+
  geom_point(data = deg_cod_predictors_pcrt[!is.na(deg_cod_predictors_pcrt$lab),],pch=18,
             mapping = aes(x = idx, y = signedpcomb))+
  geom_text_repel(data = deg_cod_predictors_pcrt[!is.na(deg_cod_predictors_pcrt$lab),],
                  size=2.5,segment.size=0.1,max.iter =1000,max.overlaps = 10,#nudge_x = -1,
                  hjust=1,direction= "y",
                  nudge_y = ifelse(deg_cod_predictors_pcrt$log2fc_rec[!is.na(deg_cod_predictors_pcrt$lab)]>0,5,-15),
                  mapping = aes(x = idx, y = signedpcomb,label = gene_ini))+
  theme_classic()+coord_flip()+lims(y=c(-16,10))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  labs(x = 'genes', y = 'signed log10(P)', fill = '')


ego_up_cod <- enrichGO(gene       = deg_cod$gene[deg_cod$type=='up' & deg_cod$type_rec=='up'],
                       universe      = deg_cod$gene,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "ALL",
                       keyType       = 'SYMBOL',
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1)
head(summary(ego_up_cod),n=10)
dotplot(ego_up_cod, size = "GeneRatio", x='p.adjust',showCategory = 10,label_format = 40)+scale_x_log10()
cnetplot(ego_up_cod, categorySize="pvalue",circular = F, colorEdge = TRUE)

ego_dw_cod<- enrichGO(gene       = deg_cod$gene[deg_cod$type=='down' & deg_cod$type_rec=='down'],
                      universe      = deg_cod$gene,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      keyType       = 'SYMBOL',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.1)
head(summary(ego_dw_cod),n=10)
dotplot(ego_dw_cod, size = "GeneRatio", x='p.adjust',showCategory = 10,label_format = 40)+scale_x_log10()
cnetplot(ego_dw_cod, categorySize="pvalue",circular = F, colorEdge = TRUE)


#fgsea
deg_cod_ini = deg_cod_ini[order(deg_cod_ini$tstat, decreasing = T),]

fcrnk = deg_cod_ini$tstat
names(fcrnk) = deg_cod_ini$gene  
fgseaRes <- fgsea(pathways = hmpathways, stats = fcrnk, minSize=15,maxSize=500)

plotfgsea(pathway = hmpathways[["HALLMARK_MYC_TARGETS_V2"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_MYC_TARGETS_V2',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')
plotfgsea(pathway = hmpathways[["HALLMARK_E2F_TARGETS"]],stats = fcrnk,fgseaRes = fgseaRes,
          gene.set.name ='HALLMARK_E2F_TARGETS',class.name = '',posClass = 'pre-HM',negClass = 'pre-NHM')


#
fgseaRes = fgseaRes[order(sign(fgseaRes$NES)*log10(fgseaRes$padj), decreasing = T),]
fgseaRes$pathway = factor(fgseaRes$pathway, levels = fgseaRes$pathway)
fgseaRes$pathway1 = gsub("HALLMARK_","", fgseaRes$pathway)
fgseaRes$pathway1  = factor(fgseaRes$pathway1, levels = fgseaRes$pathway1)
ggplot(aes(x = pathway1, y = NES, fill = NES),data = fgseaRes[fgseaRes$padj<0.05,])+
  geom_bar(stat = 'identity', color = 'black') + 
  scale_fill_gradient2(low = 'blue',high = 'red')+
  coord_flip() +
  theme_classic()+labs(x='')+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))



upgns = deg_cod$gene[which(deg_cod$type=='up' & deg_cod$log2fc> 0.5)]
dngns = deg_cod$gene[which(deg_cod$type=='down'& deg_cod$log2fc< -0.5 )]

prenhm_cod2 = cl_cod$ID1[cl_cod$Hypermutation.of.Recurrent.Tumor=='No' & cl_cod$TMZ.Treatment.Before.Recurrence==1]
m1h = gei[upgns,names(gei) %in% prehm_cod]
m1n = gei[upgns,names(gei) %in% prenhm_cod2]

m2h = gei[dngns,names(gei) %in% prehm_cod]
m2n = gei[dngns,names(gei) %in% prenhm_cod2]

mat = rbind(cbind(m1h,m1n), cbind(m2h,m2n))
#mat = t(scale(t(log2(mat+1))))
gns2mark = c(upgns[upgns %in% hmpathways$HALLMARK_MYC_TARGETS_V1])
#pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHnon_deg.heatmap.pdf',width =4,height = 2.5)
ppp<- Heatmap(mat, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
              name = "scaled_expr", cluster_columns = F,
              clustering_method_rows = 'ward.D',
              row_split = 2,column_split = c(rep('pre-HM',2),rep('pre-NHM',14)),
              row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
              #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
              show_column_names = FALSE, #width = unit(8, "cm"),
              heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gns2mark), 
                                 labels = rownames(mat)[rownames(mat) %in% gns2mark], 
                                 labels_gp = gpar(fontsize = 6))) 
ppp

nhm_cod2 = cl_cod$ID2[cl_cod$Hypermutation.of.Recurrent.Tumor=='No' & cl_cod$TMZ.Treatment.Before.Recurrence==1]
m1rh = ger[upgns,names(ger) %in% hm_cod]; m1rn = ger[upgns,names(ger) %in% nhm_cod2]
m2rh = ger[dngns,names(ger) %in% hm_cod]; m2rn = ger[dngns,names(ger) %in% nhm_cod2]
mat2 = rbind(cbind(m1rh,m1rn), cbind(m2rh,m2rn))

pppr<- Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
               name = "scaled_expr", cluster_columns = F,
               clustering_method_rows = 'ward.D',
               row_split = 2,column_split = c(rep('HM',2),rep('NHM',19)),
               row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
               #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
               show_column_names = FALSE, #width = unit(8, "cm"),
               heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gns2mark), 
                                 labels = rownames(mat)[rownames(mat) %in% gns2mark], 
                                 labels_gp = gpar(fontsize = 6))) 
pppr


d1 = row_dend(ppp)[[1]] #the up-regulated ones
cutree(d1,2) ->tc1 #cut into 2 clusters
genes1.1 = names(tc1[tc1==1])
ego1.1 <- enrichGO(gene       = genes1.1,
                   universe      = rownames(gei),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego1.1),n=10)
dotplot(ego1.1)
cnetplot(ego1.1, categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
#gn2label =c('ERCC2','DNMT1','POLD1','POLE','MCM5','E2F1','AKT1')
gn2label = c('MTHFD1','CDT1','CDKN2C','POLA1','TK1')
genes1.2 = names(tc1[tc1==2])
ego1.2 <- enrichGO(gene       = genes1.2,
                   universe      = rownames(gei),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego1.2))
dotplot(ego1.2)
cnetplot(ego1.2, categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
#gn2label = c(gn2label, 'MKI67','PLK1','CENPE','UBE2C','KIF18B','AURKB')
gn2label = c(gn2label, 'MCM5','MKI67','NCAPH2','HNRNPL','STX10')

d2 = row_dend(ppp)[[2]] #the up regulated ones
cutree(d2,2) ->tc2  #divide into one subclusters
genes2.1 = names(tc2[tc2==1])
ego2.1 <- enrichGO(gene          = genes2.1,
                   universe      = rownames(gei),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego2.1))
dotplot(ego2.1)
cnetplot(ego2.1, categorySize="pvalue",circular = F, colorEdge = TRUE)
gn2label = c(gn2label, 'COLCA2','TYRP1')

genes2.2 = names(tc2[tc2==2])
ego2.2 <- enrichGO(gene          = genes2.2,
                   universe      = rownames(gei),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   keyType       = 'SYMBOL',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
head(summary(ego2.2))
dotplot(ego2.2)
cnetplot(ego2.2, categorySize="pvalue",circular = F, colorEdge = TRUE)

gn2label = c(gn2label,c('LACC1','PTPN3','PRRT1','GRIN3A','STX1A'))
pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel_deg.heatmap0830.pdf',width =3,height = 2.5)
Heatmap(mat, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
        name = "scaled_expr", cluster_columns = F,use_raster=F,
        clustering_method_rows = 'ward.D',
        row_split = 2,column_split = c(rep('pre-HM',2),rep('pre-NHM',14)),
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,#border_gp = gpar(color = 'yellow'),
        #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
        show_column_names = FALSE, #width = unit(8, "cm"),
        heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gn2label), 
                                 labels = rownames(mat)[rownames(mat) %in% gn2label], 
                                 labels_gp = gpar(fontsize = 5))) 
dev.off()

pdf('~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel_deg.recurrence.heatmap0902.pdf',width =3,height = 2.5)
pppr<- Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("#3288bd", "white", "#d53e4f")), 
               name = "scaled_expr", cluster_columns = F,use_raster=F,
               clustering_method_rows = 'ward.D',
               row_split = 2,column_split = c(rep('HM',2),rep('NHM',19)),
               row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
               #row_split = c(rep('A',122),rep('B',40)), column_split = c(rep('C',12),rep('D',93)),
               show_column_names = FALSE, #width = unit(8, "cm"),
               heatmap_legend_param = list(title = "Exp")) +
  rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gn2label), 
                                 labels = rownames(mat)[rownames(mat) %in% gn2label], 
                                 labels_gp = gpar(fontsize = 6))) 
pppr
dev.off()

write.table(deg_wt, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.deg.txt', row.names = F, quote = F, sep = "\t")
write.table(deg_non, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmut-noncodel.deg.txt', row.names = F, quote = F, sep = "\t")
write.table(deg_cod, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmut-codel.deg.txt', row.names = F, quote = F, sep = "\t")

write.table(deg_wt_predictors, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.RNApredictors.txt', row.names = F, quote = F, sep = "\t")
write.table(deg_non_predictors, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmut-noncodel.RNApredictors.txt', row.names = F, quote = F, sep = "\t")
write.table(deg_cod_predictors, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmut-codel.RNApredictors.txt', row.names = F, quote = F, sep = "\t")


###
#===========write out the data for GSEA analysis============
###

writeGCT <-function(expdata, gctfile ){
  num_rows = nrow(expdata); num_col = ncol(expdata) #expdata: each row is a gene, and each column is a sample. gene name as row name
  
  line1<-'#1.2'
  line2<-paste(num_rows,num_col,sep='\t')
  line3<-paste(names(expdata),collapse='\t')
  line3<-paste('NAME','Description',line3, sep='\t')
  write.table(x = line1, file = gctfile,quote=F,col.name=F,row.name=F) #first line for the gct file
  write.table(x = line2, file = gctfile,quote=F,col.names=F,row.names=F,append=T)  # 2nd line for the gct file
  write.table(x = line3, file = gctfile,quote=F,col.names=F,row.names=F,append=T)  # 3rd line for the gct file
  
  dat = cbind(rownames(expdata),rep('.',num_rows),expdata)
  write.table(x = dat, file = gctfile, sep='\t',quote=F,row.names = F, col.name=F, append=T)
}

writeCLS <-function(samplelevels, clsfile ){ #temporary, not quite generalizable
  num_sample = length(samplelevels) #level correspond to each sample. e.g: c('HM','HM','NHM','NHM','NHM')
  num_level = length(unique(samplelevels)) #unique levels
  
  line1<-paste(num_sample,num_level,1,sep = "\t")
  line2<-paste0('#\t',paste(unique(samplelevels),collapse = '\t'))
  line3<-paste(samplelevels,collapse='\t')
  write.table(x = line1, file = clsfile,quote=F,col.name=F,row.name=F) #first line for the gct file
  write.table(x = line2, file = clsfile,quote=F,col.names=F,row.names=F,append=T)  # 2nd line for the gct file
  write.table(x = line3, file = clsfile,quote=F,col.names=F,row.names=F,append=T)  # 3rd line for the gct file
}

#IDHwt
gei_wt = cbind(gei[,names(gei) %in% prehm_wt], gei[,names(gei) %in% prenhm_wt])
writeGCT(expdata = gei_wt, gctfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.geneexpression.preHMpreNHM.gct')
writeCLS(samplelevels = ifelse(names(gei_wt) %in% prehm_wt,'preHM','preNHM'),clsfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.geneexpression.preHMpreNHM.cls')


ger_wt = cbind(ger[,names(ger) %in% hm_wt], ger[,names(ger) %in% nhm_wt])
writeGCT(expdata = ger_wt, gctfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.geneexpression.HMNHM.gct')
writeCLS(samplelevels = ifelse(names(ger_wt) %in% hm_wt,'HM','NHM'),clsfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.geneexpression.HMNHM.cls')

#IDHmut-noncodel
gei_non = cbind(gei[,names(gei) %in% prehm_non], gei[,names(gei) %in% prenhm_non])
writeGCT(expdata = gei_non, gctfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel.geneexpression.preHMpreNHM.gct')
writeCLS(samplelevels = ifelse(names(gei_non) %in% prehm_non,'preHM','preNHM'),clsfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel.geneexpression.preHMpreNHM.cls')


ger_non = cbind(ger[,names(ger) %in% hm_non], ger[,names(ger) %in% nhm_non])
writeGCT(expdata = ger_non, gctfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel.geneexpression.HMNHM.gct')
writeCLS(samplelevels = ifelse(names(ger_non) %in% hm_non,'HM','NHM'),clsfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel.geneexpression.HMNHM.cls')

#IDHmut-codel
gei_cod = cbind(gei[,names(gei) %in% prehm_cod], gei[,names(gei) %in% prenhm_cod])
writeGCT(expdata = gei_cod, gctfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel.geneexpression.preHMpreNHM.gct')
writeCLS(samplelevels = ifelse(names(gei_cod) %in% prehm_cod,'preHM','preNHM'),clsfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel.geneexpression.preHMpreNHM.cls')


ger_cod = cbind(ger[,names(ger) %in% hm_cod], ger[,names(ger) %in% nhm_cod])
writeGCT(expdata = ger_cod, gctfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel.geneexpression.HMNHM.gct')
writeCLS(samplelevels = ifelse(names(ger_cod) %in% hm_cod,'HM','NHM'),clsfile = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel.geneexpression.HMNHM.cls')

#:::command lines for the GSEA analysis::::
#gsea-cli.sh GSEA -res /Users/kevinmu/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.geneexpression.preHMpreNHM.gct -cls /Users/kevinmu/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHwt.geneexpression.preHMpreNHM.cls#preHM_versus_preNHM -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.5.1.symbols.gmt -collapse Collapse -mode Max_probe -norm None -nperm 1000 -permute gene_set -rnd_seed timestamp -rnd_type no_balance -scoring_scheme weighted -rpt_label IDHWT.preHMvspreNHM -metric Diff_of_Classes -sort real -order descending -chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_Gene_Symbol_with_Remapping_MSigDB.v7.5.1.chip -create_gcts false -create_svgs true -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /Users/kevinmu/gsea_home/output/sep06
#gsea-cli.sh GSEA -res /Users/kevinmu/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel.geneexpression.preHMpreNHM.gct -cls /Users/kevinmu/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutnoncodel.geneexpression.preHMpreNHM.cls#preHM_versus_preNHM -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.5.1.symbols.gmt -collapse Collapse -mode Max_probe -norm None -nperm 1000 -permute gene_set -rnd_seed timestamp -rnd_type no_balance -scoring_scheme weighted -rpt_label IDHMUTNONCODEL.preHMvspreNHM -metric Diff_of_Classes -sort real -order descending -chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_Gene_Symbol_with_Remapping_MSigDB.v7.5.1.chip -create_gcts false -create_svgs true -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /Users/kevinmu/gsea_home/output/sep06
#gsea-cli.sh GSEA -res /Users/kevinmu/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel.geneexpression.preHMpreNHM.gct -cls /Users/kevinmu/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature Genetics/code/IDHmutcodel.geneexpression.preHMpreNHM.cls#preHM_versus_preNHM -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.5.1.symbols.gmt -collapse Collapse -mode Max_probe -norm None -nperm 1000 -permute gene_set -rnd_seed timestamp -rnd_type no_balance -scoring_scheme weighted -rpt_label IDHMUTCOD.preHMvspreNHM -metric Diff_of_Classes -sort real -order descending -chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_Gene_Symbol_with_Remapping_MSigDB.v7.5.1.chip -create_gcts false -create_svgs true -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /Users/kevinmu/gsea_home/output/sep06

idhwt_gsea_prehm = read.delim('~/gsea_home/output/sep06/IDHWT.preHMvspreNHM.Gsea.touse/gsea_report_for_preHM_1662465686320.tsv')
idhwt_gsea_prenhm = read.delim('~/gsea_home/output/sep06/IDHWT.preHMvspreNHM.Gsea.touse/gsea_report_for_preNHM_1662465686320.tsv')

idhwt_gsea = rbind(idhwt_gsea_prehm,idhwt_gsea_prenhm)
idhwt_gsea$FDR.q.val[idhwt_gsea$FDR.q.val==0] = 0.001
idhwt_gsea = idhwt_gsea[order(idhwt_gsea$NES),]
idhwt_gsea$NAME = gsub("HALLMARK_","",idhwt_gsea$NAME)
idhwt_gsea$NAME = factor(idhwt_gsea$NAME, levels = idhwt_gsea$NAME)
ggplot(idhwt_gsea[idhwt_gsea$FDR.q.val<0.05,], aes(y = NAME, x = NES))+
  geom_bar(stat = 'identity', aes(fill = sign(NES)*-log10(FDR.q.val)))+
  theme_classic()+scale_fill_gradient2(low = 'blue',high = 'red')+
  labs(y = '', fill = 'FDR q')+theme(axis.text = element_text(color = 'black'))

idhnon_gsea_prehm = read.delim('~/gsea_home/output/sep06/IDHMUTNONCODEL.preHMvspreNHM.Gsea.touse/gsea_report_for_preHM_1662465871344.tsv')
idhnon_gsea_prenhm = read.delim('~/gsea_home/output/sep06/IDHMUTNONCODEL.preHMvspreNHM.Gsea.touse/gsea_report_for_preNHM_1662465871344.tsv')

idhnon_gsea = rbind(idhnon_gsea_prehm,idhnon_gsea_prenhm)
idhnon_gsea$FDR.q.val[idhnon_gsea$FDR.q.val==0] = 0.001
idhnon_gsea = idhnon_gsea[order(idhnon_gsea$NES),]
idhnon_gsea$NAME = gsub("HALLMARK_","",idhnon_gsea$NAME)
idhnon_gsea$NAME = factor(idhnon_gsea$NAME, levels = idhnon_gsea$NAME)
ggplot(idhnon_gsea[idhnon_gsea$FDR.q.val<0.05,], aes(y = NAME, x = NES))+
  geom_bar(stat = 'identity', aes(fill = sign(NES)*-log10(FDR.q.val)))+
  theme_classic()+scale_fill_gradient2(low = 'blue',high = 'red')+
  labs(y = '', fill = 'FDR q')+theme(axis.text = element_text(color = 'black'))

idhcod_gsea_prehm = read.delim('~/gsea_home/output/sep06/IDHMUTCOD.preHMvspreNHM.Gsea.touse/gsea_report_for_preHM_1662465097011.tsv')
idhcod_gsea_prenhm = read.delim('~/gsea_home/output/sep06/IDHMUTCOD.preHMvspreNHM.Gsea.touse/gsea_report_for_preNHM_1662465097011.tsv')

idhcod_gsea = rbind(idhcod_gsea_prehm,idhcod_gsea_prenhm)
idhcod_gsea$FDR.q.val[idhcod_gsea$FDR.q.val==0] = 0.001
idhcod_gsea = idhcod_gsea[order(idhcod_gsea$NES),]
idhcod_gsea$NAME = gsub("HALLMARK_","",idhcod_gsea$NAME)
idhcod_gsea$NAME = factor(idhcod_gsea$NAME, levels = idhcod_gsea$NAME)
ggplot(idhcod_gsea[idhcod_gsea$FDR.q.val<0.05,], aes(y = NAME, x = NES))+
  geom_bar(stat = 'identity', aes(fill = sign(NES)*-log10(FDR.q.val)))+
  theme_classic()+scale_fill_gradient2(low = 'blue',high = 'red')+
  labs(y = '', fill = 'FDR q')+theme(axis.text = element_text(color = 'black'))
