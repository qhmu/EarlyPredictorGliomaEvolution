#setwd('~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/')
library(ggplot2)
library(ggarchery)
options(stringsAsFactors = F)
myf <-function(n, direction = 'right'){ #helper function to define the coordinates for the oval-shape plot
  r = n
  if(direction=='right'){
    df1 = data.frame( angle = seq(length.out=n, from = pi/12, to = 11*pi/12),  group = 'A')
    df1$y = -r*cos(df1$angle)
    df1$x = r*sin(df1$angle)
    #df1$x = df1$x + 2*(max(df1$x)-df1$x)
    df1$index = 1:nrow(df1)
    return(df1)
  }else{
    df2 = data.frame( angle = seq(length.out=n, from = 13*pi/12, to = 23*pi/12),  group = 'B')
    df2$y = r*cos(df2$angle)
    df2$x = r*sin(df2$angle)
    #df2$x = df2$x - 2*(df2$x - min(df2$x))
    df2$index = 1:nrow(df2)
    return(df2)
  }
}

ini = read.delim('~/Dropbox/2020-PanGliomaEvolution/clinicalmutationCNAs/Pairedgliomas.genomics_initial_20220204.txt', na.strings = c('NA','#N/A'))
rec = read.delim('~/Dropbox/2020-PanGliomaEvolution/clinicalmutationCNAs/Pairedgliomas.genomics_recurrence_20220204.txt', na.strings = c('NA','#N/A'))
names(ini)[which(names(ini)=='chr17p_NLOH')]='CNLOH_17p'
names(rec)[which(names(rec)=='chr17p_NLOH')]='CNLOH_17p'

stopifnot(identical(ini$Patient_ID, rec$Patient_ID))

c710 = read.delim('~/Downloads/panglioma.chr710.new.tsv')
ini$g7l10 = c710$new_chr710_ini[match(ini$Patient_ID, c710$Patient_ID)]
rec$g7l10 = c710$new_chr710_rec[match(rec$Patient_ID, c710$Patient_ID)]

gns = c("IDH1.2","TP53","TP53del","ATRX",#"TERTp",
        "CIC","FUBP1", 
        "MDM2amp","MDM4amp","CDKN2Adel","CDK4amp","RB1","RB1del",
        "EGFR","EGFRamp","PDGFRA","PDGFRAamp","METamp",
        #"TMZ_R",
        #"METex14","EGFRvIII","F3T3","ZM",
        "NF1","NF1del",
        "PTEN","PTENdel","PIK3CA","PIK3CG","PIK3R1", "MYC_gain","g7l10","CNLOH_17p","codel")#

ini$histoGrade =ifelse(ini$X1p19qcodel=="codel", ifelse(ini$WHO2016_grade_initial=='III',1,0), 
                             ifelse(ini$WHO2016_grade_initial=='IV',1,0))

rec$histoGrade =ifelse(rec$X1p19qcodel=="codel", ifelse(rec$WHO2016_grade_recurrence=='III',1,0), 
                             ifelse(rec$WHO2016_grade_recurrence=='IV',1,0))

ini$WHO2021_grade_initial = substr(ini$WHO2021Classification_initial,start = nchar(ini$WHO2021Classification_initial), stop = nchar(ini$WHO2021Classification_initial))
ini$WHO2021_grade_initial = as.integer(ini$WHO2021_grade_initial)
ini$highestgrade2021 =ifelse(ini$X1p19qcodel=="codel", ifelse(ini$WHO2021_grade_initial==3,1,0), 
                             ifelse(ini$WHO2021_grade_initial==4,1,0))

rec$WHO2021_grade_recurrence= substr(rec$WHO2021Classification_recurrence,start = nchar(ini$WHO2021Classification_recurrence), stop = nchar(ini$WHO2021Classification_recurrence))
rec$WHO2021_grade_recurrence = as.integer(rec$WHO2021_grade_recurrence)
rec$highestgrade2021 =ifelse(rec$X1p19qcodel=="codel", ifelse(rec$WHO2021_grade_recurrence==3,1,0), 
                             ifelse(rec$WHO2021_grade_recurrence==4,1,0))


ini$HM = ifelse(ini$HM_I=='Yes',1,0)
rec$HM = ifelse(rec$HM_R=='Yes',1,0)

##IDH wildtype
gns_wt = gns[!gns%in% c('IDH1.2',"codel","CIC",'FUBP1','METamp','CNLOH_17p','PIK3CG')] #these features are irrelavant to the IDHwt
gns_wt = c(gns_wt,'histoGrade','HM')

df = expand.grid(early = gns_wt, late = gns_wt)
df$early = as.character(df$early)
df$late = as.character(df$late)

df$pvalue = NA; df$Coef = NA

df$early_mut = 0; df$early_wt = 0
df$late_mut= 0; df$late_wt =0
df$early_af = 0; df$late_af = 0

ini_wt = ini[which(ini$IDH =="IDHwt"),gns_wt]
rec_wt = rec[which(rec$IDH =="IDHwt"),gns_wt]

for (i in 1:nrow(df)){
  x = ini_wt[,df$early[i]]
  y = rec_wt[,df$late[i]]
  tb = table(x,y)
  if (dim(tb)[1] ==2 & dim(tb)[2]==2){
    df$pvalue[i] = fisher.test(tb)$p.value
    df$Coef[i] = log(fisher.test(tb)$estimate)
  }
  #md = glm(y ~ x, data = data.frame(x,y), family = binomial)
  #coef = summary(md)$coef
  # if (nrow(coef)==2){
  #   df$Coef[i] = coef[2,1]
  #   df$pvalue[i] = coef[2,4]
  # }else{
  #   df$Coef[i] = NA
  #   df$pvalue[i] = NA
  # }
  
  df$early_mut[i] = sum(x==1, na.rm = T); df$early_wt[i] =sum(x==0,na.rm = T)
  df$late_mut[i] = sum(y==1, na.rm = T); df$late_wt[i] =sum(y==0,na.rm = T)
  
  df$early_af[i] = sum(x==1,na.rm = T)/sum(!is.na(x))
  df$late_af[i] = sum(y==1,na.rm = T)/sum(!is.na(y))
}

idx = 1:length(gns_wt); names(idx)=gns_wt
df$early_y = idx[df$early]
df$late_y = idx[df$late]
library(ggplot2)
library(ggrepel)

df_pos = df#[df$late %in% c('histologicalgrade','highestgrade2021','hmstatus'),]
df_pos$lab1 = df_pos$early; df_pos$lab1[duplicated(df_pos$lab1)]=NA
df_pos$lab2 = df_pos$late; df_pos$lab2[duplicated(df_pos$lab2)]=NA
#df_pos = df_pos[df_pos$early!=df_pos$late,]

#change the shape
dfacc1 = myf(n = length(idx), direction = 'left'); dfacc2 = myf(n = length(idx), direction = 'right')
dfacc1$name = names(idx); dfacc2$name = names(idx)

df_pos$x1 = dfacc1$x[match(df_pos$early, dfacc1$name)]
df_pos$early_y = dfacc1$y[match(df_pos$early, dfacc1$name)]
df_pos$x2 = dfacc2$x[match(df_pos$late, dfacc2$name)]
df_pos$late_y = dfacc2$y[match(df_pos$late, dfacc2$name)]
#
ggplot(df_pos)+
  geom_point(aes(x = x1, y = early_y, size = early_af, fill = early), shape =21,show.legend = F)+
  geom_point(aes(x = x2, y = late_y, size = late_af, fill = late), shape=21,show.legend = F)+
  geom_segment(mapping = aes(x = x1, y = early_y, xend = x2, yend = late_y, alpha = -log10(pvalue),color = Coef>0), 
               data = df_pos[which(df_pos$pvalue<0.05),], show.legend = F)+
  scale_alpha(range = c(0.25,1))+#scale_color_manual(values = c('blue','red'))+
  scale_x_continuous(expand = expansion(mult = 0.4))+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  scale_size_area(limits = c(0,1))+
  geom_text_repel(data = df_pos, mapping = aes(x = x1, y = early_y, label = lab1), 
                  force_pull = 0, nudge_x = -15,direction = "y",
                  hjust = 0,angle = 0,size=2, segment.size = 0.25, max.iter = 1e4)+
  geom_text_repel(data = df_pos, mapping = aes(x = x2, y = late_y, label = lab2), 
                  force_pull = 0, nudge_x = 15,direction = "y",
                  hjust = 0,angle = 0,size=2, segment.size = 0.25, max.iter = 1e4)+
  theme_minimal()+theme(panel.grid = element_blank(), axis.text = element_blank(),
                        plot.title = element_text(hjust = 0.5))+labs(x = '', y = '', title = 'IDHwt')


df1_wt = df[df$early==df$late,]
df1_wt$alt_af = df1_wt$late_af - df1_wt$early_af
df1_wt$lab = ifelse(df1_wt$early %in% c('histoGrade','HM'), df1_wt$early, NA)
#df1_wt = df1_wt[df1_wt$early_af>0.05 | df1_wt$late_af>0.05,]
df1_wt$prop.p = NA
for (i in 1:nrow(df1_wt)){
  df1_wt$prop.p[i] = prop.test(x = c(df1_wt$early_mut[i], df1_wt$late_mut[i]),
                            n = c((df1_wt$early_mut[i] + df1_wt$early_wt[i]),(df1_wt$late_mut[i] + df1_wt$late_wt[i])))$p.value
}
df1_wt$prop.p[is.na(df1_wt$prop.p)]=1 #fix the NA issue
df1_wt$sigornot = ifelse(df1_wt$prop.p>=0.05,  'nosig', ifelse(df1_wt$late_af>df1_wt$early_af,'pos','neg'))
df1_wt$logp = -log10(df1_wt$prop.p); df1_wt$logp[df1_wt$logp>3]=3
ggplot(df1_wt, aes(x = 100*early_af, y = 100*late_af, label = early))+
  geom_point(aes(color = logp, size = -log10(prop.p)))+
  scale_colour_viridis_c(guide = "none")+
  geom_abline( intercept = 0,slope =1,lty=2)+
  geom_text_repel(size= 3, segment.size=0.25,nudge_y = 5,min.segment.length = 0, max.overlaps = 100)+
  theme_classic()+theme(axis.text = element_text(colour = 'black'))+
  lims(x = c(0,105), y =c(0,105))+
  labs(x = 'Frequency at initial (%)', y = 'Frequency at recurrence (%)', size = '-log10(P)')

df1_wt$n11 = df1_wt$n10 = df1_wt$n01 = df1_wt$n00 = 0
for (i in 1:nrow(df1_wt)){
  df1_wt$n00[i] = length(which(ini_wt[,df1_wt$early[i]]==0 &rec_wt[,df1_wt$early[i]]==0))
  df1_wt$n01[i] = length(which(ini_wt[,df1_wt$early[i]]==0 &rec_wt[,df1_wt$early[i]]==1))
  df1_wt$n10[i] = length(which(ini_wt[,df1_wt$early[i]]==1 &rec_wt[,df1_wt$early[i]]==0))
  df1_wt$n11[i] = length(which(ini_wt[,df1_wt$early[i]]==1 &rec_wt[,df1_wt$early[i]]==1))
}


for (i in 1:nrow(df1_wt)){ #use the df1_wt as an example
  gn_i = df1_wt$late[i]
  fts = df_pos$early[which(df_pos$late==gn_i & df_pos$pvalue<0.05)]
  data4mod = ini_wt
  data4mod$target = rec_wt[,gn_i]
  if (length(fts)>1){
    df = data4mod[,c(fts,'target')]
    #drop the features that are highly correlated with the feature itself, as this will be meaningless
    if (gn_i %in% fts){
      fts_selected = c(gn_i,'target')
      fts_cad = fts[!fts %in% c(gn_i, 'target')]
      for (ft in fts_cad){
        if (fisher.test(table(df[,c(ft,gn_i)]))$p.value>0.1){
          fts_selected = append(fts_selected, ft)
        }
      }
    }else{
      fts_selected=c(fts,'target')
    }
    df = df[,fts_selected]
    
    df = df[apply(as.data.frame(df[,-ncol(df)]),1,function(x) sum(!is.na(x))>0),]#remove the samples for which all values are missing
    df[is.na(df)]=0#fill the missing value with 0
    #
    mvfit_idhwt <- glm(target ~ ., data = df, family = binomial)#
    as.data.frame(summary(mvfit_idhwt)$coef) ->res_from_mod
    res_from_mod = res_from_mod[-1,] #discard the intercept
    res_from_mod$Feature = rownames(res_from_mod)
    names(res_from_mod)=c('Estimate','Stderr','Zvalue','Pvalue','Feature')
    res_from_mod$Target = gn_i
    
  }else if (length(fts)==1){
   df =  data4mod[,c(fts,'target')]
   tmp = fisher.test(table(df))
   res_from_mod = data.frame(Estimate = tmp$estimate, Stderr = NA, Zvalue = NA, Pvalue = tmp$p.value, Feature = fts, Target = gn_i)
  }
  if (i ==1){
    res_from_mod_wt = res_from_mod
  }else{
    res_from_mod_wt = rbind(res_from_mod_wt, res_from_mod)
  }
}
res_from_mod_wt0 = res_from_mod_wt
res_from_mod_wt0$univar_pval = 1
for (i in 1:nrow(res_from_mod_wt0)){
  res_from_mod_wt0$univar_pval[i] = df_pos$pvalue[which(df_pos$early==res_from_mod_wt0$Feature[i] & df_pos$late==res_from_mod_wt0$Target[i])]
}

df_pos_wt = df_pos
df_pos_wt$multivar_pval = NA
for (i in 1:nrow(df_pos_wt)){
  if (length(res_from_mod_wt0$Pvalue[which(res_from_mod_wt0$Feature==df_pos_wt$early[i] & res_from_mod_wt0$Target==df_pos_wt$late[i])])>0){
    df_pos_wt$multivar_pval[i] =res_from_mod_wt0$Pvalue[which(res_from_mod_wt0$Feature==df_pos_wt$early[i] & res_from_mod_wt0$Target==df_pos_wt$late[i])]
  }
}

res_from_mod_wt = res_from_mod_wt[res_from_mod_wt$Pvalue<0.05,]
res_from_mod_wt$direction  = ifelse(res_from_mod_wt$Estimate>0,'positive','negative')

res_from_mod_wt$early_af = df1_wt$early_af[match(res_from_mod_wt$Feature, df1_wt$early)]
res_from_mod_wt$late_af = df1_wt$late_af[match(res_from_mod_wt$Target, df1_wt$late)]

#write.table(df_pos[which(df_pos$pvalue<0.05),], file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/res_from_univar_mod_wt.txt',
#            row.names = F, quote = F, sep = "\t")

#write.table(res_from_mod_wt, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/res_from_multivar_mod_wt.txt',
#            row.names = F, quote = F, sep = "\t")
#write.table(df1_wt, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/node_feature_wt.txt',
#            row.names = F, quote = F, sep = "\t")

fts = df1_wt$early[order(-df1_wt$prop.p, -df1_wt$early_af)]
idx = 1:length(fts); names(idx) = fts
res_from_mod_wt$x1 = idx[res_from_mod_wt$Feature]; res_from_mod_wt$x2 = idx[res_from_mod_wt$Target]

#change the shape
dfacc1 = myf(n = length(idx), direction = 'left'); dfacc2 = myf(n = length(idx), direction = 'right')
dfacc1$name = names(idx); dfacc2$name = names(idx)

df1_wt$x1 = dfacc1$x[match(df1_wt$early, dfacc1$name)]
df1_wt$x2 = dfacc2$x[match(df1_wt$early, dfacc2$name)]
df1_wt$y = dfacc1$y[match(df1_wt$early, dfacc1$name)]

res_from_mod_wt$x1 = dfacc1$x[match(res_from_mod_wt$Feature, dfacc1$name)];
res_from_mod_wt$y1 = dfacc1$y[match(res_from_mod_wt$Feature, dfacc1$name)];
res_from_mod_wt$x2 = dfacc2$x[match(res_from_mod_wt$Target, dfacc2$name)];
res_from_mod_wt$y2 = dfacc2$y[match(res_from_mod_wt$Target, dfacc2$name)];
res_from_mod_wt$type = ifelse(res_from_mod_wt$Feature==res_from_mod_wt$Target,'self','cross')
#
library(ggarchery)
df1_wt$logprop = ifelse(df1_wt$prop.p<0.001, 3,-log10(df1_wt$prop.p))
df1_wt$early[df1_wt$early=='g7l10']='+7/-10'
df1_wt$late[df1_wt$late=='g7l10']='+7/-10'


uni_wt = df_pos[df_pos$pvalue<0.05,]
uni_wt = uni_wt[uni_wt$early %in% df1_wt$early & uni_wt$late %in% df1_wt$late,]
uni_wt$x1 = df1_wt$x1[match(uni_wt$early, df1_wt$early)]
uni_wt$x2 = df1_wt$x2[match(uni_wt$late, df1_wt$late)]
uni_wt$y1 = df1_wt$y[match(uni_wt$early, df1_wt$early)]
uni_wt$y2 = df1_wt$y[match(uni_wt$late, df1_wt$late)]

multi_wt = res_from_mod_wt[,c('Feature','Target','Pvalue','Estimate','x1','x2','y1','y2')]
uni_wt = uni_wt[,c('early','late','pvalue','Coef','x1','x2','y1','y2')]
names(uni_wt) = names(multi_wt)
multi_wt$uid = paste(multi_wt$Feature, multi_wt$Target, sep = "--")
uni_wt$uid = paste(uni_wt$Feature, uni_wt$Target, sep = "--")

IDHwt_uni_multi = rbind(uni_wt,multi_wt)
IDHwt_uni_multi$type = #ifelse(IDHwt_uni_multi$Feature==IDHwt_uni_multi$Target,'self',
                              ifelse(IDHwt_uni_multi$uid %in% multi_wt$uid,'multi', 'uni')#)
IDHwt_uni_multi$type[IDHwt_uni_multi$Feature==IDHwt_uni_multi$Target & IDHwt_uni_multi$type=='multi']='self'
ggplot()+
  # geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.25,
  #                   size=0.5,data = IDHwt_uni_multi[IDHwt_uni_multi$type=='uni',], show.legend = T,arrow_positions = 0.975,color= '#999999',arrow_fills='#999999',
  #                   arrows = arrow(angle = ifelse(IDHwt_uni_multi$Estimate[IDHwt_uni_multi$type=='uni']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
  #                   position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.5,
                    size=0.5,data = IDHwt_uni_multi[IDHwt_uni_multi$type=='self',], show.legend = T,arrow_positions = 0.975,
                    color= '#e5c494',arrow_fills='#e5c494',
                    arrows = arrow(angle = ifelse(IDHwt_uni_multi$Estimate[IDHwt_uni_multi$type=='self']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
                    position = position_attractsegment(start_shave = 0.01, end_shave = 0.01))+
  geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.75,
                    size=0.5,data = IDHwt_uni_multi[IDHwt_uni_multi$type=='multi',], show.legend = T,arrow_positions = 0.975,
                    color= 'indianred3',arrow_fills='indianred3',
                    arrows = arrow(angle = ifelse(IDHwt_uni_multi$Estimate[IDHwt_uni_multi$type=='multi']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
                    position = position_attractsegment(start_shave = 0.01, end_shave = 0.01))+
  # scale_alpha_manual(values = c(0.7,0.35,0.14))+
  # geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),
  #                           size=0.5,data = res_from_mod_wt, show.legend = T,arrow_positions = 0.975,color= 'indianred3',arrow_fills='indianred3',
  #                     arrows = arrow(angle = ifelse(res_from_mod_wt$Estimate>0,30,90), type = 'closed',length = unit(0.05, "inches")),
  #                     position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  #scale_alpha_manual(values = c(0.7,0.35))+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  scale_fill_viridis_c()+
  geom_point(data = df1_wt[df1_wt$early_af>0.01,], mapping = aes(x = x1, y = y, size = early_af, fill = logprop), shape =21,show.legend = F)+
  geom_point(data = df1_wt[df1_wt$late_af>0.01,], mapping = aes(x = x2, y = y, size = late_af, fill = logprop), shape=21,show.legend = T)+
  scale_size(limits = c(0,1))+
  scale_x_continuous(expand = expansion(mult = 0.3))+
  geom_text_repel(data = df1_wt[df1_wt$early_af>0.01,], mapping = aes(x = x1, y = y, label = early), force_pull = 0, 
                  nudge_x = -10,direction = "y", hjust = 0,angle = 0,size=2.5, segment.size = 0.25, max.iter = 1e4)+
  geom_text_repel(data = df1_wt[df1_wt$late_af>0.01,], mapping = aes(x = x2, y = y, label = late), force_pull = 0, 
                  nudge_x = 3,direction = "y", hjust = 0,angle = 0,size=2.5, segment.size = 0.25, max.iter = 1e4)+
  theme_minimal()+theme(panel.grid = element_blank(), axis.text = element_blank(),plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', title = 'IDHwt', color = 'Positive\npredictor')

#
##IDH mutant non codel
#
gns_non = gns[!gns%in% c('IDH1.2',"codel","CIC",'FUBP1','CDK4amp','MDM2amp','MDM4amp','g7l10','EGFRamp','PIK3CG','PIK3R1','EGFR','PDGFRA','NF1del','PTENdel')] # irrelavant to the IDHmut-noncodel
gns_non = c(gns_non,'histoGrade','HM')

df = expand.grid(early = gns_non, late = gns_non)
df$early = as.character(df$early)
df$late = as.character(df$late)

df$pvalue = NA; df$Coef = NA

df$early_mut = 0; df$early_non = 0
df$late_mut= 0; df$late_non =0
df$early_af = 0; df$late_af = 0

ini_non = ini[which(ini$IDH =="IDHmut" & ini$X1p19qcodel=="noncodel"),gns_non]
rec_non = rec[which(rec$IDH =="IDHmut" & rec$X1p19qcodel=="noncodel"),gns_non]

for (i in 1:nrow(df)){
  x = ini_non[,df$early[i]]
  y = rec_non[,df$late[i]]
  tb = table(x,y)
  if (dim(tb)[1] ==2 & dim(tb)[2]==2){
    df$pvalue[i] = fisher.test(tb)$p.value
    df$Coef[i] = log(fisher.test(tb)$estimate)
  }
  #md = glm(y ~ x, data = data.frame(x,y), family = binomial)
  #coef = summary(md)$coef
  # if (nrow(coef)==2){
  #   df$Coef[i] = coef[2,1]
  #   df$pvalue[i] = coef[2,4]
  # }else{
  #   df$Coef[i] = NA
  #   df$pvalue[i] = NA
  # }
  
  df$early_mut[i] = sum(x==1, na.rm = T); df$early_non[i] =sum(x==0,na.rm = T)
  df$late_mut[i] = sum(y==1, na.rm = T); df$late_non[i] =sum(y==0,na.rm = T)
  
  df$early_af[i] = sum(x==1,na.rm = T)/sum(!is.na(x))
  df$late_af[i] = sum(y==1,na.rm = T)/sum(!is.na(y))
}

idx = 1:length(gns_non); names(idx)=gns_non
df$early_y = idx[df$early]
df$late_y = idx[df$late]
library(ggplot2)
library(ggrepel)

df_pos = df#[df$late %in% c('histologicalgrade','highestgrade2021','hmstatus'),]
df_pos$lab1 = df_pos$early; df_pos$lab1[duplicated(df_pos$lab1)]=NA
df_pos$lab2 = df_pos$late; df_pos$lab2[duplicated(df_pos$lab2)]=NA
#df_pos = df_pos[df_pos$early!=df_pos$late,]

#change the shape
dfacc1 = myf(n = length(idx), direction = 'left'); dfacc2 = myf(n = length(idx), direction = 'right')
dfacc1$name = names(idx); dfacc2$name = names(idx)

df_pos$x1 = dfacc1$x[match(df_pos$early, dfacc1$name)]
df_pos$early_y = dfacc1$y[match(df_pos$early, dfacc1$name)]
df_pos$x2 = dfacc2$x[match(df_pos$late, dfacc2$name)]
df_pos$late_y = dfacc2$y[match(df_pos$late, dfacc2$name)]
#
ggplot(df_pos)+
  geom_point(aes(x = x1, y = early_y, size = early_af, fill = early), shape =21,show.legend = F)+
  geom_point(aes(x = x2, y = late_y, size = late_af, fill = late), shape=21,show.legend = F)+
  geom_segment(mapping = aes(x = x1, y = early_y, xend = x2, yend = late_y, alpha = -log10(pvalue),color = Coef>0), 
               data = df_pos[which(df_pos$pvalue<0.05),], show.legend = F)+
  scale_alpha(range = c(0.25,1))+#scale_color_manual(values = c('blue','red'))+
  scale_x_continuous(expand = expansion(mult = 0.4))+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  geom_text_repel(data = df_pos, mapping = aes(x = x1, y = early_y, label = lab1), 
                  force_pull = 0, nudge_x = -15,direction = "y",
                  hjust = 0,angle = 0,size=2, segment.size = 0.25, max.iter = 1e4)+
  geom_text_repel(data = df_pos, mapping = aes(x = x2, y = late_y, label = lab2), 
                  force_pull = 0, nudge_x = 15,direction = "y",
                  hjust = 0,angle = 0,size=2, segment.size = 0.25, max.iter = 1e4)+
  theme_minimal()+theme(panel.grid = element_blank(), axis.text = element_blank(),
                        plot.title = element_text(hjust = 0.5))+labs(x = '', y = '', title = 'IDHmut-noncodel')

df1_non = df[df$early==df$late,]
df1_non$alt_af = df1_non$late_af - df1_non$early_af
df1_non$lab = ifelse(df1_non$early %in% c('histoGrade','HM'), df1_non$early, NA)
#df1_non = df1_non[df1_non$early_af>0.05 | df1_non$late_af>0.05,]
df1_non$prop.p = NA
for (i in 1:nrow(df1_non)){
  df1_non$prop.p[i] = prop.test(x = c(df1_non$early_mut[i], df1_non$late_mut[i]),
                                n = c((df1_non$early_mut[i] + df1_non$early_non[i]),(df1_non$late_mut[i] + df1_non$late_non[i])))$p.value
}
df1_non$prop.p[is.na(df1_non$prop.p)]=1 #fix the NA issue
df1_non$sigornot = ifelse(df1_non$prop.p>=0.05,  'nosig', ifelse(df1_non$late_af>df1_non$early_af,'pos','neg'))
df1_non$logp = -log10(df1_non$prop.p); df1_non$logp[df1_non$logp>3]=3
ggplot(df1_non, aes(x = 100*early_af, y = 100*late_af, label = early))+
  geom_abline( intercept = 0,slope =1,lty=2,col = '#dddddd')+
  geom_point(aes(fill = logp, size = -log10(prop.p)),pch=21)+
  scale_fill_viridis_c(guide = 'none')+
  #geom_text_repel(size= 3, segment.size=0.25,nudge_y = 5,min.segment.length = 0, max.overlaps = 100)+
  theme_classic()+theme(axis.text = element_text(colour = 'black'))+
  lims(x = c(0,105), y =c(0,105))+
  labs(x = 'Frequency at initial (%)', y = 'Frequency at recurrence (%)', size = '-log10(P)')

ggplot(df1_non, aes(x = 100*early_af, y = 100*late_af, label = early))+
  geom_abline( intercept = 0,slope =1,lty=2,col = '#dddddd')+
  geom_point(aes(fill = logp, size = -log10(prop.p)),pch=21, show.legend = F)+
  scale_fill_viridis_c(guide = 'none')+
  #geom_text_repel(size= 3, segment.size=0.25,nudge_y = 5,min.segment.length = 0, max.overlaps = 100)+
  theme_classic()+theme(axis.text = element_text(color = 'black'), axis.title = element_blank())+
  scale_x_continuous(breaks = c(0,10), limits = c(0,14))+
  scale_y_continuous(breaks = c(0,10), limits = c(0,14))


df1_non$n11 = df1_non$n10 = df1_non$n01 = df1_non$n00 = 0
for (i in 1:nrow(df1_non)){
  df1_non$n00[i] = length(which(ini_non[,df1_non$early[i]]==0 &rec_non[,df1_non$early[i]]==0))
  df1_non$n01[i] = length(which(ini_non[,df1_non$early[i]]==0 &rec_non[,df1_non$early[i]]==1))
  df1_non$n10[i] = length(which(ini_non[,df1_non$early[i]]==1 &rec_non[,df1_non$early[i]]==0))
  df1_non$n11[i] = length(which(ini_non[,df1_non$early[i]]==1 &rec_non[,df1_non$early[i]]==1))
}


for (i in 1:nrow(df1_non)){ #use the df1_non as an example
  gn_i = df1_non$late[i]
  fts = df_pos$early[which(df_pos$late==gn_i & df_pos$pvalue<0.05)]
  data4mod = ini_non
  data4mod$target = rec_non[,gn_i]
  if (length(fts)>1){
    df = data4mod[,c(fts,'target')]
    #drop the features that are highly correlated with the feature itself, as this will be meaningless
    if (gn_i %in% fts){
      fts_selected = c(gn_i,'target')
      fts_cad = fts[!fts %in% c(gn_i, 'target')]
      for (ft in fts_cad){
        if (fisher.test(table(df[,c(ft,gn_i)]))$p.value>0.1){
          fts_selected = append(fts_selected, ft)
        }
      }
    }else{
      fts_selected=c(fts,'target')
    }
    df = df[,fts_selected]
    #
    if (length(fts_selected)>2){
      df = df[apply(as.data.frame(df[,-ncol(df)]),1,function(x) sum(!is.na(x))>0),]#remove the samples for which all values are missing
      df[is.na(df)]=0#fill the missing value with 0
      mvfit_idhnon <- glm(target ~ ., data = df, family = binomial)#
      as.data.frame(summary(mvfit_idhnon)$coef) -> res_from_mod
      res_from_mod = res_from_mod[-1,] #discard the intercept
      res_from_mod$Feature = rownames(res_from_mod)
      names(res_from_mod)=c('Estimate','Stderr','Zvalue','Pvalue','Feature')
      res_from_mod$Target = gn_i
    }else{
      tmp = fisher.test(table(df))
      res_from_mod = data.frame(Estimate = tmp$estimate, Stderr = NA, Zvalue = NA, Pvalue = tmp$p.value, Feature = fts, Target = gn_i)
    }
    
  }else if (length(fts)==1){
    df =  data4mod[,c(fts,'target')]
    tmp = fisher.test(table(df))
    res_from_mod = data.frame(Estimate = tmp$estimate, Stderr = NA, Zvalue = NA, Pvalue = tmp$p.value, Feature = fts, Target = gn_i)
  }
  if (i ==1){
    res_from_mod_non = res_from_mod
  }else{
    res_from_mod_non = rbind(res_from_mod_non, res_from_mod)
  }
}

res_from_mod_non0 = res_from_mod_non
res_from_mod_non0$univar_pval = 1
for (i in 1:nrow(res_from_mod_non0)){
  res_from_mod_non0$univar_pval[i] = df_pos$pvalue[which(df_pos$early==res_from_mod_non0$Feature[i] & df_pos$late==res_from_mod_non0$Target[i])]
}

df_pos_non = df_pos
df_pos_non$multivar_pval = NA
for (i in 1:nrow(df_pos_non)){
  if (length(res_from_mod_non0$Pvalue[which(res_from_mod_non0$Feature==df_pos_non$early[i] & res_from_mod_non0$Target==df_pos_non$late[i])])>0){
    df_pos_non$multivar_pval[i] =res_from_mod_non0$Pvalue[which(res_from_mod_non0$Feature==df_pos_non$early[i] & res_from_mod_non0$Target==df_pos_non$late[i])]
  }
}

res_from_mod_non = res_from_mod_non[res_from_mod_non$Pvalue<0.05,]
res_from_mod_non$direction  = ifelse(res_from_mod_non$Estimate>0,'positive','negative')

res_from_mod_non$early_af = df1_non$early_af[match(res_from_mod_non$Feature, df1_non$early)]
res_from_mod_non$late_af = df1_non$late_af[match(res_from_mod_non$Target, df1_non$late)]

#write.table(df_pos[which(df_pos$pvalue<0.05),], file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/res_from_univar_mod_non.txt',
#            row.names = F, quote = F, sep = "\t")

#write.table(res_from_mod_non, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/res_from_multivar_mod_non.txt',
#            row.names = F, quote = F, sep = "\t")
#write.table(df1_non, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/node_feature_non.txt',
#            row.names = F, quote = F, sep = "\t")

fts = df1_non$early[order(-df1_non$prop.p, -df1_non$early_af)]
idx = 1:length(fts); names(idx) = fts
res_from_mod_non$x1 = idx[res_from_mod_non$Feature]; res_from_mod_non$x2 = idx[res_from_mod_non$Target]

#change the shape
dfacc1 = myf(n = length(idx), direction = 'left'); dfacc2 = myf(n = length(idx), direction = 'right')
dfacc1$name = names(idx); dfacc2$name = names(idx)

df1_non$x1 = dfacc1$x[match(df1_non$early, dfacc1$name)]
df1_non$x2 = dfacc2$x[match(df1_non$early, dfacc2$name)]
df1_non$y = dfacc1$y[match(df1_non$early, dfacc1$name)]

res_from_mod_non$x1 = dfacc1$x[match(res_from_mod_non$Feature, dfacc1$name)];
res_from_mod_non$y1 = dfacc1$y[match(res_from_mod_non$Feature, dfacc1$name)];
res_from_mod_non$x2 = dfacc2$x[match(res_from_mod_non$Target, dfacc2$name)];
res_from_mod_non$y2 = dfacc2$y[match(res_from_mod_non$Target, dfacc2$name)];
res_from_mod_non$type = ifelse(res_from_mod_non$Feature==res_from_mod_non$Target,'self','cross')
#
library(ggarchery)
df1_non$logprop = ifelse(df1_non$prop.p<0.001, 3,-log10(df1_non$prop.p))
df1_non$early[df1_non$early=='CNLOH_17p']='17p_CNLOH'
df1_non$late[df1_non$late=='CNLOH_17p']='17p_CNLOH'


uni_non = df_pos[df_pos$pvalue<0.05,]
uni_non = uni_non[uni_non$early %in% df1_non$early & uni_non$late %in% df1_non$late,]
uni_non$x1 = df1_non$x1[match(uni_non$early, df1_non$early)]
uni_non$x2 = df1_non$x2[match(uni_non$late, df1_non$late)]
uni_non$y1 = df1_non$y[match(uni_non$early, df1_non$early)]
uni_non$y2 = df1_non$y[match(uni_non$late, df1_non$late)]

multi_non = res_from_mod_non[,c('Feature','Target','Pvalue','Estimate','x1','x2','y1','y2')]
uni_non = uni_non[,c('early','late','pvalue','Coef','x1','x2','y1','y2')]
names(uni_non) = names(multi_non)
multi_non$uid = paste(multi_non$Feature, multi_non$Target, sep = "--")
uni_non$uid = paste(uni_non$Feature, uni_non$Target, sep = "--")

IDHnon_uni_multi = rbind(uni_non,multi_non)
IDHnon_uni_multi$type = #ifelse(IDHnon_uni_multi$Feature==IDHnon_uni_multi$Target,'self',
  ifelse(IDHnon_uni_multi$uid %in% multi_non$uid,'multi', 'uni')#)
IDHnon_uni_multi$type[IDHnon_uni_multi$Feature==IDHnon_uni_multi$Target & IDHnon_uni_multi$type=='multi']='self'
ggplot()+
  # geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.25,
  #                   size=0.5,data = IDHnon_uni_multi[IDHnon_uni_multi$type=='uni',], show.legend = T,arrow_positions = 0.975,color= '#999999',arrow_fills='#999999',
  #                   arrows = arrow(angle = ifelse(IDHnon_uni_multi$Estimate[IDHnon_uni_multi$type=='uni']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
  #                   position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.5,
                    size=0.5,data = IDHnon_uni_multi[IDHnon_uni_multi$type=='self',], show.legend = T,arrow_positions = 0.975,color= '#e5c494',arrow_fills='#e5c494',
                    arrows = arrow(angle = ifelse(IDHnon_uni_multi$Estimate[IDHnon_uni_multi$type=='self']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
                    position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.75,
                    size=0.5,data = IDHnon_uni_multi[IDHnon_uni_multi$type=='multi',], show.legend = T,arrow_positions = 0.975,color= 'indianred3',arrow_fills='indianred3',
                    arrows = arrow(angle = ifelse(IDHnon_uni_multi$Estimate[IDHnon_uni_multi$type=='multi']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
                    position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  # scale_alpha_manual(values = c(0.7,0.35,0.14))+
  # geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),
  #                           size=0.5,data = res_from_mod_non, show.legend = T,arrow_positions = 0.975,color= 'indianred3',arrow_fills='indianred3',
  #                     arrows = arrow(angle = ifelse(res_from_mod_non$Estimate>0,30,90), type = 'closed',length = unit(0.05, "inches")),
  #                     position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  #scale_alpha_manual(values = c(0.7,0.35))+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  scale_fill_viridis_c()+
  geom_point(data = df1_non[df1_non$early_af>0.0,], mapping = aes(x = x1, y = y, size = early_af, fill = logprop), shape =21,show.legend = F)+
  geom_point(data = df1_non[df1_non$late_af>0.0,], mapping = aes(x = x2, y = y, size = late_af, fill = logprop), shape=21,show.legend = T)+
  scale_size(limits = c(0,1))+
  scale_x_continuous(expand = expansion(mult = 0.3))+
  geom_text_repel(data = df1_non[df1_non$early_af>0.0,], mapping = aes(x = x1, y = y, label = early), force_pull = 0, 
                  nudge_x = -7,direction = "y", hjust = 0,angle = 0,size=2.5, segment.size = 0.25, max.iter = 1e4)+
  geom_text_repel(data = df1_non[df1_non$late_af>0.0,], mapping = aes(x = x2, y = y, label = late), force_pull = 0, 
                  nudge_x = 3,direction = "y", hjust = 0,angle = 0,size=2.5, segment.size = 0.25, max.iter = 1e4)+
  theme_minimal()+theme(panel.grid = element_blank(), axis.text = element_blank(),plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', title = 'IDHmut-noncodel', color = 'Positive\npredictor')


##IDH mutant, 1p/19q codel #
gns_cod = gns[!gns%in% c('IDH1.2',"codel","g7l10",'MDM2amp','MDM4amp','EGFRamp','EGFR','PDGFRA','PDGFRAamp','METamp','CNLOH_17p','MYC_gain','NF1del','TP53del','CDK4amp','ATRX','PTEN','NF1')]
gns_cod = c(gns_cod,'histoGrade','HM')


df = expand.grid(early = gns_cod, late = gns_cod)
df$early = as.character(df$early)
df$late = as.character(df$late)

df$pvalue = NA; df$Coef = NA

df$early_mut = 0; df$early_cod = 0
df$late_mut= 0; df$late_cod =0
df$early_af = 0; df$late_af = 0

ini_cod = ini[which(ini$IDH =="IDHmut" & ini$X1p19qcodel=="codel"),gns_cod]
rec_cod = rec[which(rec$IDH =="IDHmut" & rec$X1p19qcodel=="codel"),gns_cod]

for (i in 1:nrow(df)){
  x = ini_cod[,df$early[i]]
  y = rec_cod[,df$late[i]]
  tb = table(x,y)
  if (dim(tb)[1] ==2 & dim(tb)[2]==2){
    df$pvalue[i] = fisher.test(tb)$p.value
    df$Coef[i] = log(fisher.test(tb)$estimate)
  }
  #md = glm(y ~ x, data = data.frame(x,y), family = binomial)
  #coef = summary(md)$coef
  # if (nrow(coef)==2){
  #   df$Coef[i] = coef[2,1]
  #   df$pvalue[i] = coef[2,4]
  # }else{
  #   df$Coef[i] = NA
  #   df$pvalue[i] = NA
  # }
  
  df$early_mut[i] = sum(x==1, na.rm = T); df$early_cod[i] =sum(x==0,na.rm = T)
  df$late_mut[i] = sum(y==1, na.rm = T); df$late_cod[i] =sum(y==0,na.rm = T)
  
  df$early_af[i] = sum(x==1,na.rm = T)/sum(!is.na(x))
  df$late_af[i] = sum(y==1,na.rm = T)/sum(!is.na(y))
}

idx = 1:length(gns_cod); names(idx)=gns_cod
df$early_y = idx[df$early]
df$late_y = idx[df$late]
library(ggplot2)
library(ggrepel)

df_pos = df#[df$late %in% c('histologicalgrade','highestgrade2021','hmstatus'),]
df_pos$lab1 = df_pos$early; df_pos$lab1[duplicated(df_pos$lab1)]=NA
df_pos$lab2 = df_pos$late; df_pos$lab2[duplicated(df_pos$lab2)]=NA
#df_pos = df_pos[df_pos$early!=df_pos$late,]

#change the shape
dfacc1 = myf(n = length(idx), direction = 'left'); dfacc2 = myf(n = length(idx), direction = 'right')
dfacc1$name = names(idx); dfacc2$name = names(idx)

df_pos$x1 = dfacc1$x[match(df_pos$early, dfacc1$name)]
df_pos$early_y = dfacc1$y[match(df_pos$early, dfacc1$name)]
df_pos$x2 = dfacc2$x[match(df_pos$late, dfacc2$name)]
df_pos$late_y = dfacc2$y[match(df_pos$late, dfacc2$name)]
#
ggplot(df_pos)+
  geom_point(aes(x = x1, y = early_y, size = early_af, fill = early), shape =21,show.legend = F)+
  geom_point(aes(x = x2, y = late_y, size = late_af, fill = late), shape=21,show.legend = F)+
  geom_segment(mapping = aes(x = x1, y = early_y, xend = x2, yend = late_y, alpha = -log10(pvalue),color = Coef>0), 
               data = df_pos[which(df_pos$pvalue<0.05),], show.legend = F)+
  scale_alpha(range = c(0.25,1))+#scale_color_manual(values = c('blue','red'))+
  scale_x_continuous(expand = expansion(mult = 0.4))+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  geom_text_repel(data = df_pos, mapping = aes(x = x1, y = early_y, label = lab1), 
                  force_pull = 0, nudge_x = -15,direction = "y",
                  hjust = 0,angle = 0,size=2, segment.size = 0.25, max.iter = 1e4)+
  geom_text_repel(data = df_pos, mapping = aes(x = x2, y = late_y, label = lab2), 
                  force_pull = 0, nudge_x = 15,direction = "y",
                  hjust = 0,angle = 0,size=2, segment.size = 0.25, max.iter = 1e4)+
  theme_minimal()+theme(panel.grid = element_blank(), axis.text = element_blank(),
                        plot.title = element_text(hjust = 0.5))+labs(x = '', y = '', title = 'IDHmut-codel')

df1_cod = df[df$early==df$late,]
df1_cod$alt_af = df1_cod$late_af - df1_cod$early_af
df1_cod$lab = ifelse(df1_cod$early %in% c('histoGrade','HM'), df1_cod$early, NA)
#df1_cod = df1_cod[df1_cod$early_af>0.05 | df1_cod$late_af>0.05,]
df1_cod$prop.p = NA
for (i in 1:nrow(df1_cod)){
  df1_cod$prop.p[i] = prop.test(x = c(df1_cod$early_mut[i], df1_cod$late_mut[i]),
                                n = c((df1_cod$early_mut[i] + df1_cod$early_cod[i]),(df1_cod$late_mut[i] + df1_cod$late_cod[i])))$p.value
}
df1_cod$prop.p[is.na(df1_cod$prop.p)]=1 #fix the NA issue
df1_cod$sigornot = ifelse(df1_cod$prop.p>=0.05,  'nosig', ifelse(df1_cod$late_af>df1_cod$early_af,'pos','neg'))
df1_cod$logp = -log10(df1_cod$prop.p); df1_cod$logp[df1_cod$logp>3]=3
ggplot(df1_cod, aes(x = 100*early_af, y = 100*late_af, label = early))+
  geom_point(aes(fill = logp, size = -log10(prop.p)),pch=21)+
  scale_fill_viridis_c(guide = "none")+
  geom_abline( intercept = 0,slope =1,lty=2)+
  #geom_text_repel(size= 3, segment.size=0.25,nudge_y = 5,min.segment.length = 0, max.overlaps = 100)+
  theme_classic()+theme(axis.text = element_text(colour = 'black'))+
  lims(x = c(0,105), y =c(0,105))+
  labs(x = 'Frequency at initial (%)', y = 'Frequency at recurrence (%)', size = '-log10(P)')

df1_cod$n11 = df1_cod$n10 = df1_cod$n01 = df1_cod$n00 = 0
for (i in 1:nrow(df1_cod)){
  df1_cod$n00[i] = length(which(ini_cod[,df1_cod$early[i]]==0 &rec_cod[,df1_cod$early[i]]==0))
  df1_cod$n01[i] = length(which(ini_cod[,df1_cod$early[i]]==0 &rec_cod[,df1_cod$early[i]]==1))
  df1_cod$n10[i] = length(which(ini_cod[,df1_cod$early[i]]==1 &rec_cod[,df1_cod$early[i]]==0))
  df1_cod$n11[i] = length(which(ini_cod[,df1_cod$early[i]]==1 &rec_cod[,df1_cod$early[i]]==1))
}


for (i in 1:nrow(df1_cod)){ #use the df1_cod as an example
  gn_i = df1_cod$late[i]
  fts = df_pos$early[which(df_pos$late==gn_i & df_pos$pvalue<0.05)]
  data4mod = ini_cod
  data4mod$target = rec_cod[,gn_i]
  if (length(fts)>1){
    dfi = data4mod[,c(fts,'target')]
    #drop the features that are highly correlated with the feature itself, as this will be meaningless
    if (gn_i %in% fts){
      fts_selected = c(gn_i,'target')
      fts_cad = fts[!fts %in% c(gn_i, 'target')]
      for (ft in fts_cad){
        if (fisher.test(table(dfi[,c(ft,gn_i)]))$p.value>0.1){
          fts_selected = append(fts_selected, ft)
        }
      }
    }else{
      fts_selected=c(fts,'target')
    }
    dfi = dfi[,fts_selected]
    #
    if (length(fts_selected)>2){
      mvfit_IDHcod <- glm(target ~ ., data = dfi, family = binomial)#
      as.data.frame(summary(mvfit_IDHcod)$coef) -> res_from_mod
      res_from_mod = res_from_mod[-1,] #discard the intercept
      res_from_mod$Feature = rownames(res_from_mod)
      names(res_from_mod)=c('Estimate','Stderr','Zvalue','Pvalue','Feature')
      res_from_mod$Target = gn_i
    }else{
      tmp = fisher.test(table(dfi))
      res_from_mod = data.frame(Estimate = tmp$estimate, Stderr = NA, Zvalue = NA, Pvalue = tmp$p.value, Feature = fts, Target = gn_i)
    }
    
  }else if (length(fts)==1){
    dfi =  data4mod[,c(fts,'target')]
    tmp = fisher.test(table(dfi))
    res_from_mod = data.frame(Estimate = tmp$estimate, Stderr = NA, Zvalue = NA, Pvalue = tmp$p.value, Feature = fts, Target = gn_i)
  }else{
    res_from_mod = data.frame(Estimate = NA, Stderr = NA, Zvalue = NA, Pvalue = NA, Feature = NA, Target = NA)
  }
  if (i ==1){
    res_from_mod_cod = res_from_mod
  }else{
    res_from_mod_cod = rbind(res_from_mod_cod, res_from_mod)
  }
}
res_from_mod_cod = res_from_mod_cod[!is.na(res_from_mod_cod$Feature),]

res_from_mod_cod0 = res_from_mod_cod
res_from_mod_cod0$univar_pval = 1
for (i in 1:nrow(res_from_mod_cod0)){
  res_from_mod_cod0$univar_pval[i] = df_pos$pvalue[which(df_pos$early==res_from_mod_cod0$Feature[i] & df_pos$late==res_from_mod_cod0$Target[i])]
}

df_pos_cod = df_pos
df_pos_cod$multivar_pval = NA
for (i in 1:nrow(df_pos_cod)){
  if (length(res_from_mod_cod0$Pvalue[which(res_from_mod_cod0$Feature==df_pos_cod$early[i] & res_from_mod_cod0$Target==df_pos_cod$late[i])])>0){
    df_pos_cod$multivar_pval[i] =res_from_mod_cod0$Pvalue[which(res_from_mod_cod0$Feature==df_pos_cod$early[i] & res_from_mod_cod0$Target==df_pos_cod$late[i])]
  }
}

res_from_mod_cod = res_from_mod_cod[res_from_mod_cod$Pvalue<0.05,]
res_from_mod_cod$direction  = ifelse(res_from_mod_cod$Estimate>0,'positive','negative')

res_from_mod_cod$early_af = df1_cod$early_af[match(res_from_mod_cod$Feature, df1_cod$early)]
res_from_mod_cod$late_af = df1_cod$late_af[match(res_from_mod_cod$Target, df1_cod$late)]

#write.table(df_pos[which(df_pos$pvalue<0.05),], file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/res_from_univar_mod_cod.txt',
#            row.names = F, quote = F, sep = "\t")

#write.table(res_from_mod_cod, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/res_from_multivar_mod_cod.txt',
#            row.names = F, quote = F, sep = "\t")
#write.table(df1_cod, file = '~/Dropbox/2020-PanGliomaEvolution/Manuscript/Nature\ Genetics/code/node_feature_cod.txt',
#            row.names = F, quote = F, sep = "\t")

fts = df1_cod$early[order(-df1_cod$prop.p, -df1_cod$early_af)]
idx = 1:length(fts); names(idx) = fts
res_from_mod_cod$x1 = idx[res_from_mod_cod$Feature]; res_from_mod_cod$x2 = idx[res_from_mod_cod$Target]

#change the shape
dfacc1 = myf(n = length(idx), direction = 'left'); dfacc2 = myf(n = length(idx), direction = 'right')
dfacc1$name = names(idx); dfacc2$name = names(idx)

df1_cod$x1 = dfacc1$x[match(df1_cod$early, dfacc1$name)]
df1_cod$x2 = dfacc2$x[match(df1_cod$early, dfacc2$name)]
df1_cod$y = dfacc1$y[match(df1_cod$early, dfacc1$name)]

res_from_mod_cod$x1 = dfacc1$x[match(res_from_mod_cod$Feature, dfacc1$name)];
res_from_mod_cod$y1 = dfacc1$y[match(res_from_mod_cod$Feature, dfacc1$name)];
res_from_mod_cod$x2 = dfacc2$x[match(res_from_mod_cod$Target, dfacc2$name)];
res_from_mod_cod$y2 = dfacc2$y[match(res_from_mod_cod$Target, dfacc2$name)];
res_from_mod_cod$type = ifelse(res_from_mod_cod$Feature==res_from_mod_cod$Target,'self','cross')
#
library(ggarchery)
df1_cod$logprop = ifelse(df1_cod$prop.p<0.001, 3,-log10(df1_cod$prop.p))



uni_cod = df_pos[df_pos$pvalue<0.05,]
uni_cod = uni_cod[uni_cod$early %in% df1_cod$early & uni_cod$late %in% df1_cod$late,]
uni_cod$x1 = df1_cod$x1[match(uni_cod$early, df1_cod$early)]
uni_cod$x2 = df1_cod$x2[match(uni_cod$late, df1_cod$late)]
uni_cod$y1 = df1_cod$y[match(uni_cod$early, df1_cod$early)]
uni_cod$y2 = df1_cod$y[match(uni_cod$late, df1_cod$late)]

multi_cod = res_from_mod_cod[,c('Feature','Target','Pvalue','Estimate','x1','x2','y1','y2')]
uni_cod = uni_cod[,c('early','late','pvalue','Coef','x1','x2','y1','y2')]
names(uni_cod) = names(multi_cod)
multi_cod$uid = paste(multi_cod$Feature, multi_cod$Target, sep = "--")
uni_cod$uid = paste(uni_cod$Feature, uni_cod$Target, sep = "--")

IDHcod_uni_multi = rbind(uni_cod,multi_cod)
IDHcod_uni_multi$type = #ifelse(IDHcod_uni_multi$Feature==IDHcod_uni_multi$Target,'self',
  ifelse(IDHcod_uni_multi$uid %in% multi_cod$uid,'multi', 'uni')#)
IDHcod_uni_multi$type[IDHcod_uni_multi$Feature==IDHcod_uni_multi$Target & IDHcod_uni_multi$type=='multi']='self'
ggplot()+
  # geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.25,
  #                   size=0.5,data = IDHcod_uni_multi[IDHcod_uni_multi$type=='uni',], show.legend = T,arrow_positions = 0.975,color= '#999999',arrow_fills='#999999',
  #                   arrows = arrow(angle = ifelse(IDHcod_uni_multi$Estimate[IDHcod_uni_multi$type=='uni']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
  #                   position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.5,
                    size=0.5,data = IDHcod_uni_multi[IDHcod_uni_multi$type=='self',], show.legend = T,arrow_positions = 0.975,color= '#e5c494',arrow_fills='#e5c494',
                    arrows = arrow(angle = ifelse(IDHcod_uni_multi$Estimate[IDHcod_uni_multi$type=='self']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
                    position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),alpha = 0.75,
                    size=0.5,data = IDHcod_uni_multi[IDHcod_uni_multi$type=='multi',], show.legend = T,arrow_positions = 0.95,color= 'indianred3',arrow_fills='indianred3',
                    arrows = arrow(angle = ifelse(IDHcod_uni_multi$Estimate[IDHcod_uni_multi$type=='multi']>0,30,90), type = 'closed',length = unit(0.05, "inches")),
                    position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  # scale_alpha_manual(values = c(0.7,0.35,0.14))+
  # geom_arrowsegment(mapping = aes(x = x1, y = y1, xend = x2, yend = y2, alpha = type),
  #                           size=0.5,data = res_from_mod_cod, show.legend = T,arrow_positions = 0.975,color= 'indianred3',arrow_fills='indianred3',
  #                     arrows = arrow(angle = ifelse(res_from_mod_cod$Estimate>0,30,90), type = 'closed',length = unit(0.05, "inches")),
  #                     position = position_attractsegment(start_shave = 0, end_shave = 0.025, type_shave = "distance"))+
  #scale_alpha_manual(values = c(0.7,0.35))+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  scale_fill_viridis_c()+
  geom_point(data = df1_cod[df1_cod$early_af>0.0,], mapping = aes(x = x1, y = y, size = early_af, fill = logprop), shape =21,show.legend = F)+
  geom_point(data = df1_cod[df1_cod$late_af>0.0,], mapping = aes(x = x2, y = y, size = late_af, fill = logprop), shape=21,show.legend = T)+
  #scale_size_continuous(guide = "none")+
  scale_x_continuous(expand = expansion(mult = 0.3))+
  geom_text_repel(data = df1_cod[df1_cod$early_af>0.0,], mapping = aes(x = x1, y = y, label = early), force_pull = 0, 
                  nudge_x = -5,direction = "y", hjust = 0,angle = 0,size=2.5, segment.size = 0.25, max.iter = 1e4)+
  geom_text_repel(data = df1_cod[df1_cod$late_af>0.0,], mapping = aes(x = x2, y = y, label = late), force_pull = 0, 
                  nudge_x = 3,direction = "y", hjust = 0,angle = 0,size=2.5, segment.size = 0.25, max.iter = 1e4)+
  theme_minimal()+theme(panel.grid = element_blank(), axis.text = element_blank(),plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', title = 'IDHmut-codel', color = 'Positive\npredictor')

#Alluvials 
library(alluvial)
dt = as.data.frame(table(ini_cod$HM, rec_cod$HM))
alluvial(dt[,1:2], freq=dt$Freq,gap.width = .5,
         #col = ifelse(dt$Recurrence == 2, "#c6dbef",ifelse(dt$Recurrence== 3, "#6baed6","#2171b5")),
         #col = ifelse( dt$prog=='yes', "#010078", "#fed976" ),
         #col = ifelse(dt$Subtype=="IDHwt","#107050", ifelse(dt$Subtype== "IDHcodel","#c7d34d","#55d9c0")),
         col = ifelse(dt$Var2==1,'red','#999999'),
         border = NA,alpha=0.75,hide = dt$Freq<1,
         # ordering = list(
         #   #order(dt$Grade_I,dt$Grade_R,decreasing = T),
         #   NULL,NULL,
         #   NULL,NULL
         #   #order(dt$Grade_I,dt$Grade_R, decreasing = T)
         # ), 
         blocks = T,cex = 0.001,cw = 0.3, cex.axis = 0.6)
#
dt = as.data.frame(table(ini$MYC_gain[ini$IDH=='IDHwt'],ini$TMZ_R[ini$IDH=='IDHwt'], rec$HM[rec$IDH=='IDHwt']))
alluvial(dt[,1:3], freq=dt$Freq,gap.width = .05,
         #col = ifelse(dt$Recurrence == 2, "#c6dbef",ifelse(dt$Recurrence== 3, "#6baed6","#2171b5")),
         #col = ifelse( dt$prog=='yes', "#010078", "#fed976" ),
         #col = ifelse(dt$Subtype=="IDHwt","#107050", ifelse(dt$Subtype== "IDHcodel","#c7d34d","#55d9c0")),
         col = ifelse(dt$Var3==1,'red','#999999'),
         border = NA,alpha=0.75,hide = dt$Freq<1,
         # ordering = list(
         #   #order(dt$Grade_I,dt$Grade_R,decreasing = T),
         #   NULL,NULL,
         #   NULL,NULL
         #   #order(dt$Grade_I,dt$Grade_R, decreasing = T)
         # ), 
         blocks = T,cex = 0.001,cw = 0.15, cex.axis = 0.9)

dt = as.data.frame(table(ini_wt$ATRX, rec_wt$HM))
alluvial(dt[,1:2], freq=dt$Freq,gap.width = .5,
         #col = ifelse(dt$Recurrence == 2, "#c6dbef",ifelse(dt$Recurrence== 3, "#6baed6","#2171b5")),
         #col = ifelse( dt$prog=='yes', "#010078", "#fed976" ),
         #col = ifelse(dt$Subtype=="IDHwt","#107050", ifelse(dt$Subtype== "IDHcodel","#c7d34d","#55d9c0")),
         col = ifelse(dt$Var2==1,'red','#999999'),
         border = NA,alpha=0.75,hide = dt$Freq<1,
         # ordering = list(
         #   #order(dt$Grade_I,dt$Grade_R,decreasing = T),
         #   NULL,NULL,
         #   NULL,NULL
         #   #order(dt$Grade_I,dt$Grade_R, decreasing = T)
         # ), 
         blocks = T,cex = 0.001,cw = 0.3, cex.axis = 0.6)

dt = as.data.frame(table(ini_wt$ATRX,ini_wt$MYC_gain, rec_wt$HM))
alluvial(dt[,1:3], freq=dt$Freq,gap.width = .1,
         #col = ifelse(dt$Recurrence == 2, "#c6dbef",ifelse(dt$Recurrence== 3, "#6baed6","#2171b5")),
         #col = ifelse( dt$prog=='yes', "#010078", "#fed976" ),
         #col = ifelse(dt$Subtype=="IDHwt","#107050", ifelse(dt$Subtype== "IDHcodel","#c7d34d","#55d9c0")),
         col = ifelse(dt$Var3==1,'red','#999999'),
         border = NA,alpha=0.75,hide = dt$Freq<1,
         # ordering = list(
         #   #order(dt$Grade_I,dt$Grade_R,decreasing = T),
         #   NULL,NULL,
         #   NULL,NULL
         #   #order(dt$Grade_I,dt$Grade_R, decreasing = T)
         # ), 
         blocks = T,cex = 0.001,cw = 0.15, cex.axis = 0.6)
#

dt = as.data.frame(table(ini$MYC_gain[ini$IDH=='IDHmut' & ini$X1p19qcodel=='noncodel'],
                         ini$TMZ_R[ini$IDH=='IDHmut' & ini$X1p19qcodel=='noncodel'],
                         rec$HM[rec$IDH=='IDHmut' & rec$X1p19qcodel=='noncodel']))
alluvial(dt[,1:3], freq=dt$Freq,gap.width = .1,
         #col = ifelse(dt$Recurrence == 2, "#c6dbef",ifelse(dt$Recurrence== 3, "#6baed6","#2171b5")),
         #col = ifelse( dt$prog=='yes', "#010078", "#fed976" ),
         #col = ifelse(dt$Subtype=="IDHwt","#107050", ifelse(dt$Subtype== "IDHcodel","#c7d34d","#55d9c0")),
         col = ifelse(dt$Var3==1,'red','#999999'),
         border = NA,alpha=0.75,hide = dt$Freq<1,
         # ordering = list(
         #   #order(dt$Grade_I,dt$Grade_R,decreasing = T),
         #   NULL,NULL,
         #   NULL,NULL
         #   #order(dt$Grade_I,dt$Grade_R, decreasing = T)
         # ), 
         blocks = T,cex = 0.001,cw = 0.15, cex.axis = 0.6)


tmpdf = data.frame(PDGFRAamp_ini = ini_non$PDGFRAamp,TP53_ini = ini_non$TP53,
                   PDGFRAamp_rec= rec_non$PDGFRAamp,TP53_rec = rec_non$TP53)
library(pheatmap)
tmpdf = tmpdf[apply(tmpdf[,c(1,4)], 1, function(x) sum(is.na(x)))==0,]
pheatmap(t(tmpdf[,c(1,4)]), show_colnames = F, legend = F)
