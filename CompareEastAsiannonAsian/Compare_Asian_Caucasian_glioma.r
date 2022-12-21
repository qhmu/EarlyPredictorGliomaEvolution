setwd('~/Dropbox/2020-PanGliomaEvolution/Figures/mutationLandscape/')
options(stringsAsFactors = F)
library(ggplot2)
library(cowplot)
ini = read.delim('glioma.initial.mutations_0723.txt', na.strings = c('NA','#N/A'))
ini$X1p19qcodel[ini$Patient_ID=="PS128" |ini$Patient_ID=="PS148"]="noncodel" #!!!!
ini$codel[ini$Patient_ID=="PS128" |ini$Patient_ID=="PS148"]=0
rec = read.delim('glioma.recurrence.mutations_0723.txt', na.strings = c('NA','#N/A'))
rec$X1p19qcodel[rec$Patient_ID=="PS128"|rec$Patient_ID=="PS148"]="noncodel" #!!!!
rec$codel[rec$Patient_ID=="PS128"|rec$Patient_ID=="PS148"]=0
ini$TP53del[ini$Cohorts=='Yale'] = ini$MDM2amp[ini$Cohorts=='Yale'] = ini$MDM4amp[ini$Cohorts=='Yale'] =ini$CDKN2Adel[ini$Cohorts=='Yale'] =ini$CDK4amp[ini$Cohorts=='Yale'] =ini$RB1del[ini$Cohorts=='Yale'] =ini$EGFRamp[ini$Cohorts=='Yale'] =ini$PDGFRAamp[ini$Cohorts=='Yale'] =ini$METamp[ini$Cohorts=='Yale'] =ini$NF1del[ini$Cohorts=='Yale'] =ini$PTENdel[ini$Cohorts=='Yale'] =NA
rec$TP53del[rec$Cohorts=='Yale'] = rec$MDM2amp[rec$Cohorts=='Yale'] = rec$MDM4amp[rec$Cohorts=='Yale'] =rec$CDKN2Adel[rec$Cohorts=='Yale'] =rec$CDK4amp[rec$Cohorts=='Yale'] =rec$RB1del[rec$Cohorts=='Yale'] =rec$EGFRamp[rec$Cohorts=='Yale'] =rec$PDGFRAamp[rec$Cohorts=='Yale'] =rec$METamp[rec$Cohorts=='Yale'] =rec$NF1del[rec$Cohorts=='Yale'] =rec$PTENdel[rec$Cohorts=='Yale'] =NA
stopifnot(identical(ini$Patient_ID, rec$Patient_ID))
names(ini)[names(ini)=="chr17p_NLOH"] = "chr17p_CNLOH"
names(rec)[names(rec)=="chr17p_NLOH"] = "chr17p_CNLOH"
gns = c("IDH1.2","TP53","chr17p_CNLOH","ATRX","codel","CIC","FUBP1", 
        "MDM2amp","MDM4amp","CDKN2Adel","CDK4amp","RB1","RB1del",
        "EGFR","EGFRamp","PDGFRA","PDGFRAamp","METamp","F3T3","ZM","METex14","EGFRvIII","NF1","NF1del",
        "PTEN","PIK3CA","PIK3CG","PIK3R1",
        "chr7p_gain","chr7q_gain","chr10p_loss","chr10q_loss", "MYC_gain")#
#br=read.delim('~/Documents/pangliomaevolution/Rev1/Novelty/TreatmentEffect/glioma.branch.evolution.txt',sep=' ')
inim = ini[,gns]
recm = rec[,gns]
p1m = inim
for (i in 1:nrow(p1m)){
  for (j in 1:ncol(p1m)){
    a = inim[i,j]; b = recm[i,j]
    if (is.na(a) | is.na(b)){
      p1m[i,j]=NA
    }else if (a==1 & b ==1){
      p1m[i,j] = 3
    } else if (a==1 & b ==0){
      p1m[i,j] = 1
    }else if (a==0 & b ==1){
      p1m[i,j] = 2
    }else if (a==0 & b ==0) {
      p1m[i,j] = 0
    }else{
      print(paste('WARNING: weired value detected at coordinate',i,j))
    }
  }
}
p1a = data.frame(Subtype = paste0(ini$IDH,ini$X1p19qcodel),
                 Grade_I = ini$Grade_1,
                 Grade_R = ini$Grade_2,
                 TMZ_R = ifelse(is.na(rec$TMZ_R),NA,ifelse(rec$TMZ_R==1,4,0)),
                 MGMTfusion = ifelse(is.na(rec$MGMTfusion),NA,ifelse(rec$MGMTfusion==1,4,0)),
                 Hypermutation = ifelse(is.na(rec$HM_R),NA,ifelse(rec$HM_R=="YES"|rec$HM_R=="Yes",4,0)),
                 MMR = ifelse(is.na(rec$MMR),NA,ifelse(rec$MMR==1,4,0)),
                 #MYC_gain = ifelse(is.na(ini$MYC_gain),NA,ifelse(ini$MYC_gain==1,4,0)),
                 stringsAsFactors = F)

rownames(p1a) = ini$Patient_ID
p1a$Subtype[p1a$Subtype=='IDHwtnoncodel']='IDHwt'
p1a$Subtype[p1a$Subtype=='IDHmutnoncodel']='IDHnon'
p1a$Subtype[p1a$Subtype=='IDHmutcodel']='IDHcodel'

clin=read.delim('../dataHub/glioma.clinical.txt', stringsAsFactors = F,row.names = 1)
p1a$Cohort = clin$Cohorts[match(rownames(p1a), rownames(clin))]
p1a$New =ifelse(p1a$Cohort %in% c('CGGA','Tiantan','CUHK','SMC'),"New",'Published')
p1a$Race = clin$Race[match(rownames(p1a), rownames(clin))]

#        1. compare surgical interval and survival        #
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
clin$Subtype = p1a$Subtype[match(rownames(clin), rownames(p1a))]
clin$Subtype = factor(clin$Subtype, levels = c('IDHwt','IDHnon','IDHcodel'))
ggplot(aes(x = Subtype, y = Age, color = Race), data = clin[!is.na(clin$Race),])+
  #geom_violin(trim = F, width = 0.75) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15),cex=0.5,alpha=0.7)+
  stat_compare_means(label = "p.format",paired = F,show.legend = F, method = 'wilcox.test',label.y.npc = 0.95)+
  scale_y_continuous(limits = c(15,85))+
  scale_color_manual(values = c("#d7191c","#2b83ba"))+
  theme_classic()+labs(x = '', y = 'Age (years)')

# ggplot(aes(x = Subtype, y = Surgerical_interval.month., color = Race), data = clin[!is.na(clin$Race),])+
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position=position_jitterdodge(jitter.width = 0.15),cex=0.5,alpha=0.7)+
#   stat_compare_means(label = "p.format",paired = F,show.legend = F,label.y.npc = 0.95, method = 't.test')+
#   scale_color_manual(values = c("#d7191c","#2b83ba"))+
#   theme_classic()+labs(x = '', y = 'Surgical Interval (months)')

library(survival)
library(survminer)
library(weights)
clin$censor_PFS=1
for (sbt in c('IDHwt','IDHnon','IDHcodel')){
  clin_sub <- clin[which(!is.na(clin$Race) &clin$Subtype==sbt),]
  my.Surv <- with(clin_sub,Surv(time = Surgerical_interval.month., event = censor_PFS == 1)) 
  expr.surv <- survminer::surv_fit(my.Surv ~ Race, data = clin_sub)
  smax <- max(clin_sub[ ,"Surgerical_interval.month."], na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (90*tmax)/100
  log.rank <- survdiff(my.Surv ~ Race, rho = 0, data = clin_sub)
  mantle.cox <- survdiff(my.Surv ~ Race, rho = 1, data = clin_sub)
  surv <- data.frame(summary(expr.surv)$table)
  model <- summary(coxph(my.Surv ~ Race, data=clin_sub))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
  mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
  star.log <- weights::starmaker(log.rank.p)
  star.log <- weights::starmaker(log.rank.p)
  star.mcox <- weights::starmaker(mantle.cox.p)
  legend.labs = c(sprintf("Asian (n=%s, events=%s, median=%s)", surv$records[1], surv$events[1], surv$median[1]),
                  sprintf("Non-Asian (n=%s, events=%s, median=%s)",  surv$records[2], surv$events[2], surv$median[2]))
  figure_s1f <- survminer::ggsurvplot(fit = expr.surv, censor = T, conf.int = F, legend = c(0.7,0.9),surv.median.line='hv',
                                      surv.scale = "percent", ylab = "Surviving probability", xlab = "Progression-free survival (Months)",
                                      legend.labs = legend.labs, legend.title = "", xlim = c(0,smax),
                                      font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(),
                                      data = clin_sub)
  figure_s1f$plot <-  figure_s1f$plot + annotate("text", x = xmax, y = c(0.7,0.625,0.55), size = 8/3,
                                                 label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                                           sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                                           sprintf("%s Mantle-Cox p value= %s", star.mcox, mantle.cox.p)))
  plt <-figure_s1f + labs(title = sbt); print(plt)
}
for (sbt in c('IDHwt','IDHnon','IDHcodel')){
  clin_sub <- clin[which(!is.na(clin$Race) &clin$Subtype==sbt ),] #& clin$TMZ_R==0
  my.Surv <- with(clin_sub,Surv(time = OS.month., event = Censor_OS..0.alive.1.dead. == 1)) 
  expr.surv <- survminer::surv_fit(my.Surv ~ Race, data = clin_sub)
  smax <- max(clin_sub[ ,"OS.month."], na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (90*tmax)/100
  log.rank <- survdiff(my.Surv ~ Race, rho = 0, data = clin_sub)
  mantle.cox <- survdiff(my.Surv ~ Race, rho = 1, data = clin_sub)
  surv <- data.frame(summary(expr.surv)$table)
  model <- summary(coxph(my.Surv ~ Race, data=clin_sub))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
  mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
  star.log <- weights::starmaker(log.rank.p)
  star.log <- weights::starmaker(log.rank.p)
  star.mcox <- weights::starmaker(mantle.cox.p)
  legend.labs = c(sprintf("Asian (n=%s, events=%s, median=%s)", surv$records[1], surv$events[1], surv$median[1]),
                  sprintf("Non-Asian (n=%s, events=%s, median=%s)",  surv$records[2], surv$events[2], surv$median[2]))
  figure_s1f <- survminer::ggsurvplot(fit = expr.surv, censor = T, conf.int = F, legend = c(0.7,0.9),surv.median.line='hv',
                                      surv.scale = "percent", ylab = "Surviving probability", xlab = "Overall survival (Months)",
                                      legend.labs = legend.labs, legend.title = "", xlim = c(0,smax),
                                      font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(),
                                      data = clin_sub)
  figure_s1f$plot <-  figure_s1f$plot + annotate("text", x = xmax, y = c(0.7,0.625), size = 8/3,
                                                 label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                                           sprintf("%s Log-rank p value= %s", star.log, log.rank.p)))
  plt <-figure_s1f + labs(title = sbt); print(plt)
}
clin$PPS = clin$OS.month. - clin$Surgerical_interval.month.
for (sbt in c('IDHwt','IDHnon','IDHcodel')){
  clin_sub <- clin[which(!is.na(clin$Race) &clin$Subtype==sbt),]
  my.Surv <- with(clin_sub,Surv(time = PPS, event = Censor_OS..0.alive.1.dead. == 1)) 
  expr.surv <- survminer::surv_fit(my.Surv ~ Race, data = clin_sub)
  smax <- max(clin_sub[ ,"PPS"], na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (90*tmax)/100
  log.rank <- survdiff(my.Surv ~ Race, rho = 0, data = clin_sub)
  mantle.cox <- survdiff(my.Surv ~ Race, rho = 1, data = clin_sub)
  surv <- data.frame(summary(expr.surv)$table)
  model <- summary(coxph(my.Surv ~ Race, data=clin_sub))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
  mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
  star.log <- weights::starmaker(log.rank.p)
  star.log <- weights::starmaker(log.rank.p)
  star.mcox <- weights::starmaker(mantle.cox.p)
  legend.labs = c(sprintf("Asian (n=%s, events=%s, median=%s)", surv$records[1], surv$events[1], signif(surv$median[1],2)),
                  sprintf("Non-Asian (n=%s, events=%s, median=%s)",  surv$records[2], surv$events[2], signif(surv$median[2],2)))
  figure_s1f <- survminer::ggsurvplot(fit = expr.surv, censor = T, conf.int = F, legend = c(0.7,0.9),surv.median.line='hv',
                                      surv.scale = "percent", ylab = "Surviving probability", xlab = "Post-progression survival (Months)",
                                      legend.labs = legend.labs, legend.title = "", xlim = c(0,smax),
                                      font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(),
                                      data = clin_sub)
  figure_s1f$plot <-  figure_s1f$plot + annotate("text", x = xmax, y = c(0.7,0.625), size = 8/3,
                                                 label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                                           sprintf("%s Log-rank p value= %s", star.log, log.rank.p)))
  plt <-figure_s1f + labs(title = sbt); print(plt)
}

# cox models, on OS
sbt = 'IDHwt'
clin_sub <- clin[which(!is.na(clin$Race) &clin$Subtype==sbt),]
clin_sub$Age_I = ifelse(clin_sub$Age <45,'young','old')
clin_sub$Grade_I = ifelse(clin_sub$Grade_1=="IV","GBM","LGG")
covariates <- c("Age_I", "Gender",  "Race", "Grade_I", "TMZ_R")
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(OS.month., Censor_OS..0.alive.1.dead.)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data =clin_sub )})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         #HR <- paste0(HR, " (",              HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR","HR_lo","HR_up", "wald.test", "pval")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res = as.data.frame(res); res$Var = rownames(res)
res$group = ifelse(res$pval>0.05,'nosig',ifelse(res$HR<1,'fav','inf'))
ggplot(res, aes(x = Var, y = HR))+
  geom_pointrange(aes(ymin = HR_lo, ymax = HR_up, color = group),pch=15, show.legend = F)+
  geom_hline(yintercept = 1,lty=3,col='black')+
  scale_color_manual(values = c('#5ab4ac','tomato','#cccccc'))+
  scale_y_continuous(limits = c(0,5),breaks = 0:3,labels = 0:3)+
  coord_flip()+theme_classic()+theme(axis.text = element_text(color = 'black'))+
  labs(x = '', y = 'Hazard ratio (95% CI)                   ', title = paste0(sbt,', univariate Cox model'))+
  geom_text(y = 3.5, label = paste('P =',signif(res$pval,2)), hjust = 0,size=3)
#
res.cox <- coxph(Surv(OS.month., Censor_OS..0.alive.1.dead.) ~ Age_I + Race +Gender+ Grade_I + TMZ_R, data =  clin_sub)
res = as.data.frame(summary(res.cox)$conf.int)
names(res) = c('HR','xxx','HR_lo','HR_up')
res$Var = rownames(res)
res$pval = summary(res.cox)$coefficients[,5]
res$group = ifelse(res$pval>0.05,'nosig',ifelse(res$HR<1,'fav','inf'))

ggplot(res, aes(x = Var, y = HR))+
  geom_pointrange(aes(ymin = HR_lo, ymax = HR_up, color = group),pch=15, show.legend = F)+
  geom_hline(yintercept = 1,lty=3,col='black')+
  scale_color_manual(values = c('#5ab4ac','tomato','#cccccc'))+
  scale_y_continuous(limits = c(0,5),breaks = 0:3,labels = 0:3)+
  coord_flip()+theme_classic()+theme(axis.text = element_text(color = 'black'))+
  labs(x = '', y = 'Hazard ratio (95% CI)                   ', title = paste0(sbt,', multivariate Cox model'))+
  geom_text(y = 3.5, label = paste('P =',signif(res$pval,2)), hjust = 0,size=3)

#batch effects? seems no
sbt = 'IDHwt'
clin_sub <- clin[which(!is.na(clin$Race) &clin$Subtype==sbt),]
clin_sub$Grade_I = ifelse(clin_sub$Grade_1=="IV","GBM","LGG")
clin_sub$Cohorts[clin_sub$Cohorts=="Tiantan"] = "CGGA"
clin_sub$Cohorts[clin_sub$Cohorts=="SMC_old"] = "SMC"
clin_sub$Cohorts[clin_sub$Cohorts=="TCGA-GBM"|clin_sub$Cohorts=="TCGA-LGG"] = "TCGA"
clin_sub = clin_sub[clin_sub$Grade_I=='GBM' & clin_sub$Cohorts %in% c("CGGA","SMC"),]
my.Surv <- with(clin_sub,Surv(time = PPS, event = Censor_OS..0.alive.1.dead. == 1)) 
expr.surv <- survminer::surv_fit(my.Surv ~ Cohorts, data = clin_sub)
smax <- max(clin_sub[ ,"PPS"], na.rm = TRUE)
tmax <- smax-(25*smax)/100
xmax <- (90*tmax)/100
log.rank <- survdiff(my.Surv ~ Cohorts, rho = 0, data = clin_sub)
mantle.cox <- survdiff(my.Surv ~ Cohorts, rho = 1, data = clin_sub)
surv <- data.frame(summary(expr.surv)$table)
model <- summary(coxph(my.Surv ~ Cohorts, data=clin_sub))
HR <- round(model$conf.int[1],2)
HR.lower <- round(model$conf.int[3],2)
HR.upper <- round(model$conf.int[4],2)
log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
star.log <- weights::starmaker(log.rank.p)
star.log <- weights::starmaker(log.rank.p)
star.mcox <- weights::starmaker(mantle.cox.p)
legend.labs = c(sprintf("CGGA (n=%s, events=%s, median=%s)", surv$records[1], surv$events[1], surv$median[1]),
                sprintf("SMC (n=%s, events=%s, median=%s)",  surv$records[2], surv$events[2], surv$median[2]))
figure_s1f <- survminer::ggsurvplot(fit = expr.surv, censor = T, conf.int = F, legend = c(0.7,0.95),surv.median.line='hv',
                                    surv.scale = "percent", ylab = "Surviving probability", xlab = "Post-progression survival (Months)",
                                    legend.labs = legend.labs, legend.title = "", #xlim = c(0,40),
                                    font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(),
                                    data = clin_sub)
figure_s1f$plot <-  figure_s1f$plot + annotate("text", x = xmax, y = c(0.775,0.7,0.625), size = 8/3,
                                               label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                                         sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                                         sprintf("%s Mantle-Cox p value= %s", star.mcox, mantle.cox.p)))
plt <-figure_s1f + labs(title = sbt)+lims(x = c(0,110)); print(plt)


#cptac 
cptac = read.delim('~/Dropbox/communter/Rev1/CPTAC.IDHwt.clinical.txt', stringsAsFactors = F)
my.Surv <- with(cptac,Surv(time = OS, event = CensorOS)) 
expr.surv <- survminer::surv_fit(my.Surv ~ ethnicity_group, data = cptac)
smax <- max(cptac[ ,"OS"], na.rm = TRUE)
tmax <- smax-(25*smax)/100
xmax <- (90*tmax)/100
log.rank <- survdiff(my.Surv ~ ethnicity_group, rho = 0, data = cptac)
mantle.cox <- survdiff(my.Surv ~ ethnicity_group, rho = 1, data = cptac)
surv <- data.frame(summary(expr.surv)$table)
model <- summary(coxph(my.Surv ~ ethnicity_group, data=cptac))
HR <- round(model$conf.int[1],2)
HR.lower <- round(model$conf.int[3],2)
HR.upper <- round(model$conf.int[4],2)
log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
star.log <- weights::starmaker(log.rank.p)
star.log <- weights::starmaker(log.rank.p)
star.mcox <- weights::starmaker(mantle.cox.p)
legend.labs = c(sprintf("Asian (n=%s, events=%s, median=%s)", surv$records[1], surv$events[1], signif(surv$median[1],2)),
                sprintf("non-Asian (n=%s, events=%s, median=%s)",  surv$records[2], surv$events[2], signif(surv$median[2],2)))
figure_s1f <- ggsurvplot(fit = expr.surv, censor = T, conf.int = F, legend = c(0.7,0.95),surv.median.line='hv',
                                    surv.scale = "percent", ylab = "Surviving probability", xlab = "Overall survival (Months)",
                                    legend.labs = legend.labs, legend.title = "", #xlim = c(0,smax),
                                    font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(),
                                    data = cptac)
figure_s1f$plot <-  figure_s1f$plot + annotate("text", x = xmax*1.1, y = c(0.775,0.7,0.625), size = 8/3,
                                               label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                                         sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                                         sprintf("%s Mantle-Cox p value= %s", star.mcox, mantle.cox.p)))
plt <-figure_s1f + labs(title = sbt); print(plt)

# cox models, on PPS
sbt = 'IDHwt'
clin_sub <- clin[which(!is.na(clin$Race) &clin$Subtype==sbt),]
clin_sub$Age_I = ifelse(clin_sub$Age <45,'young','old')
clin_sub$Grade_I = ifelse(clin_sub$Grade_1=="IV","GBM","LGG")
covariates <- c("Age_I", "Gender",  "Race", "Grade_I", "TMZ_R")
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(PPS, Censor_OS..0.alive.1.dead.)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data =clin_sub )})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         #HR <- paste0(HR, " (",              HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR","HR_lo","HR_up", "wald.test", "pval")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res = as.data.frame(res); res$Var = rownames(res)
res$group = ifelse(res$pval>0.05,'nosig',ifelse(res$HR<1,'fav','inf'))
ggplot(res, aes(x = Var, y = HR))+
  geom_pointrange(aes(ymin = HR_lo, ymax = HR_up, color = group),pch=15, show.legend = F)+
  geom_hline(yintercept = 1,lty=3,col='black')+
  scale_color_manual(values = c('#5ab4ac','tomato','#cccccc'))+
  scale_y_continuous(limits = c(0,5),breaks = 0:3,labels = 0:3)+
  coord_flip()+theme_classic()+theme(axis.text = element_text(color = 'black'))+
  labs(x = '', y = 'Hazard ratio (95% CI)                   ', title = paste0(sbt,', univariate Cox model'))+
  geom_text(y = 3.5, label = paste('P =',signif(res$pval,2)), hjust = 0,size=3)
#
res.cox <- coxph(Surv(OS.month., Censor_OS..0.alive.1.dead.) ~ Age_I + Race +Gender+ Grade_I + TMZ_R, data =  clin_sub)
res = as.data.frame(summary(res.cox)$conf.int)
names(res) = c('HR','xxx','HR_lo','HR_up')
res$Var = rownames(res)
res$pval = summary(res.cox)$coefficients[,5]
res$group = ifelse(res$pval>0.05,'nosig',ifelse(res$HR<1,'fav','inf'))

ggplot(res, aes(x = Var, y = HR))+
  geom_pointrange(aes(ymin = HR_lo, ymax = HR_up, color = group),pch=15, show.legend = F)+
  geom_hline(yintercept = 1,lty=3,col='black')+
  scale_color_manual(values = c('#5ab4ac','tomato','#cccccc'))+
  scale_y_continuous(limits = c(0,5),breaks = 0:3,labels = 0:3)+
  coord_flip()+theme_classic()+theme(axis.text = element_text(color = 'black'))+
  labs(x = '', y = 'Hazard ratio (95% CI)                   ', title = paste0(sbt,', multivariate Cox model'))+
  geom_text(y = 3.5, label = paste('P =',signif(res$pval,2)), hjust = 0,size=3)





#             2. compare driver mutations in Asian and non-Asian by subtype         #

lowPurity = c("PS003","PS006","PS011","PS017","PS019","PS028","PS065","PS069","PS095","PS097","PS103",
              "PS104","PS119","PS123","PS269","PS274","PS278","PS284","PS287")

idx = which((!(is.na(ini$IDH1.2))) & (!(is.na(p1a$Race)))); p1a = p1a[idx,]; p1m = p1m[idx,] #all Asian

#idx = which(ini$Cohorts=='GLSS'); p1a = p1a[idx,]; p1m = p1m[idx,] #GLASS
#idx = which(startsWith(rownames(p1a),"PSX")); p1a = p1a[idx,]; p1m = p1m[idx,] #highqual203



p1m = as.data.frame(t(p1m)) 
idx = order(p1a$Subtype,p1a$Grade_R, p1a$Grade_I, decreasing = T)
p1a = p1a[idx,]
p1m = p1m[,idx]
#p1mmyc = as.data.frame(t(p1m));# p1m = p1m[rownames(p1m)!="MYC_gain",] #MYC_gain is moved to annotation

#
scdf = read.delim('~/Desktop/Asian_nonAsian_score.txt')
scdf = scdf[order(scdf$Asian_Cauca_Dif_Founding.Events),]
scdf = scdf[!scdf$gene %in% c("IDH1.2","codel"),]
##
# IDHwt
##
table( p1a$Subtype,p1a$Race)
p1mt0 = as.data.frame(t(p1m))
gns = c("CIC","FUBP1","ATRX","TP53","chr17p_CNLOH","MDM2amp","MDM4amp",#"IDH1.2","codel",
        "CDKN2Adel","CDK4amp","RB1del","RB1",
        "EGFRamp","EGFR","EGFRvIII","PDGFRAamp","PDGFRA","METamp","ZM","METex14","F3T3",
        "NF1","NF1del","PTEN","PIK3CA","PIK3CG","PIK3R1","chr7p_gain","chr7q_gain","chr10p_loss","chr10q_loss",
        "MYC_gain","Hypermutation","MMR")
gns2 = c("CIC","FUBP1","ATRX","TP53","chr17p_CNLOH","MDM2amp","MDM4amp","p53_pathway", #"IDH1.2","codel",
        "CDKN2Adel","CDK4amp","RB1del","RB1","cellcycle_pathway",
        "EGFRamp","EGFR","EGFRvIII","PDGFRAamp","PDGFRA","METamp","ZM","METex14","F3T3","NF1","NF1del","RTK_pathway",
        "PTEN","PIK3CA","PIK3CG","PIK3R1", "PI3K_pathway",
        "chr7p_gain","chr7q_gain","chr10p_loss","chr10q_loss",
        "MYC_gain","Hypermutation","MMR")
rc = "Asian";sbt = "IDHwt"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] 
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_Asian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_Asian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_Asian) = names(p1mt)
for (i in 1:nrow(p1st_Asian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_Asian[i,1] = sum(x==0)
  p1st_Asian[i,2] = sum(x==1)
  p1st_Asian[i,3] = sum(x==2)
  p1st_Asian[i,4] = sum(x==3)
}

p1st_Asian$gene = rownames(p1st_Asian)
p1st_Asian$race = "Asian"
#hm: only consider TMZ treated cases
p1st_Asian$wt[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_Asian$wt[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_Asian$alt_ini = (p1st_Asian$ini+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)
p1st_Asian$alt_rec = (p1st_Asian$rec+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)

#
rc = "non-Asian";sbt = "IDHwt"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] #ignore MYC gain
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_nonAsian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_nonAsian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_nonAsian) = names(p1mt)
for (i in 1:nrow(p1st_nonAsian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_nonAsian[i,1] = sum(x==0)
  p1st_nonAsian[i,2] = sum(x==1)
  p1st_nonAsian[i,3] = sum(x==2)
  p1st_nonAsian[i,4] = sum(x==3)
}


p1st_nonAsian$gene = rownames(p1st_nonAsian)
p1st_nonAsian$race = "non-Asian"
#hm: only consider TMZ treated cases
p1st_nonAsian$wt[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_nonAsian$wt[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_nonAsian$alt_ini = (p1st_nonAsian$ini+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)
p1st_nonAsian$alt_rec = (p1st_nonAsian$rec+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)

p1st = rbind(p1st_Asian, p1st_nonAsian)
p1st$inip = p1st$ini/(p1st$wt+p1st$shr+p1st$ini)
p1st$recp = p1st$rec/(p1st$wt+p1st$shr+p1st$rec)
p1st$shrp = p1st$shr/(p1st$wt+p1st$shr+p1st$ini)
p1st$wtp = p1st$wt/(p1st$wt+p1st$shr+p1st$ini)
p1st_idhwt = p1st
#
library(reshape2)
tmp = melt(p1st[,c('wtp','inip','recp','shrp','gene','race')], id.vars = c('gene','race'),variable.name = 'status')
tmp$xlab = paste(tmp$gene,tmp$race, sep = '.') #lost during evolution

tmp$value[tmp$race=='non-Asian']=-tmp$value[tmp$race=='non-Asian']
tmp$race = factor(tmp$race, levels = c("non-Asian",'Asian'))

tmp_idhwt = tmp
gnorder = scdf$gene#p1st$gene[order(p1st$alt_rec[p1st$race=='non-Asian'] - p1st$alt_rec[p1st$race=='Asian'],decreasing = T)]
tmp$gene = factor(tmp$gene, levels = gnorder)
gnorder0 = gnorder
#tmp = tmp[!tmp$gene %in% c('IDH1.2','codel','CIC','FUBP1','METex14','ZM','F3T3','TERTp'),]
#tmp$gene = factor(tmp$gene)

#
df2 = data.frame(gene = rev(levels(tmp$gene)), p_alt = 1, p_evol = 1,p_overall = 1, stringsAsFactors = F)
for (ix in 1:nrow(df2)){ #comapre Asian and non-Asian
  tb = p1st[p1st$gene==df2$gene[ix],1:4]
  tb$alt = tb[,2]+tb[,3]+tb[,4]
  p = fisher.test(tb[,c(1,5)])$p.value 
  df2$p_alt[ix] = signif(p,2)
  tb$statuschange = tb[,2]+tb[,3]
  tb$nochange = tb[,1]+tb[,4]
  q = fisher.test(tb[,c(7,6)])$p.value
  df2$p_evol[ix] = signif(q,2)
  r = fisher.test(tb[,c(1,4,6)])$p.value
  r = fisher.test(tb[,1:4])$p.value
  df2$p_overall[ix] = signif(r,2)
}
df2$gene = factor(df2$gene, levels = levels(tmp$gene))
df2_idhwt = df2
df2 = reshape2::melt(df2, id.vars = c('gene'))
library(ggplot2)
p1<-ggplot(mapping = aes(x = gene, y = 100*value),data = tmp[!is.na(tmp$gene),])+
  geom_bar(mapping = aes(fill = status),stat = 'identity', show.legend = F)+
  geom_hline(yintercept = 0, size= 0.25)+
  scale_fill_manual(values = c(NA,"#d02a7c","#010078","#e2b449"))+
  scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = paste0(sbt,'') )+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0, 0, 0), "cm"))

p1_2<- ggplot(mapping = aes(x = gene, y = 90,label = value),data = df2[df2$variable=='p_overall',])+
  geom_text(size=2,position = 'identity',aes(color = value<0.1), show.legend = F)+
  scale_color_manual(values = c('black','red'))+
  coord_flip()+#facet_grid( ~ variable,scales = 'free_x')+
  theme_bw()+labs(y = '', x = '' )+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),
        panel.grid.minor.x  = element_blank(),panel.grid.major.x  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

cowplot::plot_grid(p1,p1_2,align = 'h',rel_widths = c(0.8,0.2),nrow = 1)
write.table(p1st_idhwt, file = 'p1st_idhwt.txt',quote = F,sep = "\t")
#
# now IDHmut-noncodel
#
rc = "Asian";sbt = "IDHnon"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] 
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_Asian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_Asian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_Asian) = names(p1mt)
for (i in 1:nrow(p1st_Asian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_Asian[i,1] = sum(x==0)
  p1st_Asian[i,2] = sum(x==1)
  p1st_Asian[i,3] = sum(x==2)
  p1st_Asian[i,4] = sum(x==3)
}

p1st_Asian$gene = rownames(p1st_Asian)
p1st_Asian$race = "Asian"
#hm: only consider TMZ treated cases
p1st_Asian$wt[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_Asian$wt[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_Asian$alt_ini = (p1st_Asian$ini+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)
p1st_Asian$alt_rec = (p1st_Asian$rec+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)

#
rc = "non-Asian";sbt = "IDHnon"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] #ignore MYC gain
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_nonAsian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_nonAsian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_nonAsian) = names(p1mt)
for (i in 1:nrow(p1st_nonAsian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_nonAsian[i,1] = sum(x==0)
  p1st_nonAsian[i,2] = sum(x==1)
  p1st_nonAsian[i,3] = sum(x==2)
  p1st_nonAsian[i,4] = sum(x==3)
}

p1st_nonAsian$gene = rownames(p1st_nonAsian)
p1st_nonAsian$race = "non-Asian"
#hm: only consider TMZ treated cases
p1st_nonAsian$wt[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_nonAsian$wt[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_nonAsian$alt_ini = (p1st_nonAsian$ini+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)
p1st_nonAsian$alt_rec = (p1st_nonAsian$rec+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)


p1st = rbind(p1st_Asian, p1st_nonAsian)
p1st$inip = p1st$ini/(p1st$wt+p1st$shr+p1st$ini)
p1st$recp = p1st$rec/(p1st$wt+p1st$shr+p1st$rec)
p1st$shrp = p1st$shr/(p1st$wt+p1st$shr+p1st$ini)
p1st$wtp = p1st$wt/(p1st$wt+p1st$shr+p1st$ini)
p1st_idhnon = p1st
#

tmp = melt(p1st[,c('wtp','inip','recp','shrp','gene','race')], id.vars = c('gene','race'),variable.name = 'status')
tmp$xlab = paste(tmp$gene,tmp$race, sep = '.') #lost during evolution

tmp$value[tmp$race=='non-Asian']=-tmp$value[tmp$race=='non-Asian']
tmp$race = factor(tmp$race, levels = c("non-Asian",'Asian'))
tmp_idhnon = tmp
gnorder = scdf$gene #p1st$gene[order(p1st$alt_rec[p1st$race=='non-Asian'] - p1st$alt_rec[p1st$race=='Asian'],decreasing = T)]
tmp$gene = factor(tmp$gene, levels = gnorder)
#gnorder1 = gnorder
#tmp = tmp[!tmp$gene %in% c('IDH1.2','codel','CIC','FUBP1','METex14','ZM','F3T3','TERTp'),]
#tmp$gene = factor(tmp$gene)

#
df2 = data.frame(gene = rev(levels(tmp$gene)), p_alt = 1, p_evol = 1,p_overall = 1, stringsAsFactors = F)
for (ix in 1:nrow(df2)){ #comapre Asian and non-Asian
  tb = p1st[p1st$gene==df2$gene[ix],1:4]
  tb$alt = tb[,2]+tb[,3]+tb[,4]
  p = fisher.test(tb[,c(1,5)])$p.value 
  df2$p_alt[ix] = signif(p,2)
  tb$statuschange = tb[,2]+tb[,3]
  tb$nochange = tb[,1]+tb[,4]
  q = fisher.test(tb[,c(7,6)])$p.value
  df2$p_evol[ix] = signif(q,2)
  r = fisher.test(tb[,c(1,4,6)])$p.value
  r = fisher.test(tb[,1:4])$p.value
  df2$p_overall[ix] = signif(r,2)
}
df2$gene = factor(df2$gene, levels = levels(tmp$gene))
df2_idhnon = df2
df2 = reshape2::melt(df2, id.vars = c('gene'))

p2<-ggplot(mapping = aes(x = gene, y = 100*value),data = tmp[!is.na(tmp$gene),])+
  geom_bar(mapping = aes(fill = status),stat = 'identity', show.legend = F)+
  geom_hline(yintercept = 0, size= 0.25)+
  scale_fill_manual(values = c(NA,"#d02a7c","#010078","#e2b449"))+
  scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = paste0(sbt,'') )+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0, 0, 0), "cm"))

p2_2<- ggplot(mapping = aes(x = gene, y = 90,label = value),data = df2[df2$variable=='p_overall',])+
  geom_text(size=2,position = 'identity',aes(color = value<0.1), show.legend = F)+
  scale_color_manual(values = c('black','red'))+
  coord_flip()+#facet_grid( ~ variable,scales = 'free_x')+
  theme_bw()+labs(y = '', x = '' )+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),
        panel.grid.minor.x  = element_blank(),panel.grid.major.x  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

cowplot::plot_grid(p2,p2_2,align = 'h',rel_widths = c(0.8,0.2),nrow = 1)
write.table(p1st_idhnon, file = 'p1st_idhnon.txt',quote = F, sep = "\t")
#
# now IDHmut-codel
#
rc = "Asian";sbt = "IDHcodel"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] 
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_Asian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_Asian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_Asian) = names(p1mt)
for (i in 1:nrow(p1st_Asian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_Asian[i,1] = sum(x==0)
  p1st_Asian[i,2] = sum(x==1)
  p1st_Asian[i,3] = sum(x==2)
  p1st_Asian[i,4] = sum(x==3)
}

p1st_Asian$gene = rownames(p1st_Asian)
p1st_Asian$race = "Asian"

#hm: only consider TMZ treated cases
p1st_Asian$wt[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_Asian$wt[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#

p1st_Asian$alt_ini = (p1st_Asian$ini+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)
p1st_Asian$alt_rec = (p1st_Asian$rec+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)

#
rc = "non-Asian";sbt = "IDHcodel"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] #ignore MYC gain
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_nonAsian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_nonAsian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_nonAsian) = names(p1mt)
for (i in 1:nrow(p1st_nonAsian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_nonAsian[i,1] = sum(x==0)
  p1st_nonAsian[i,2] = sum(x==1)
  p1st_nonAsian[i,3] = sum(x==2)
  p1st_nonAsian[i,4] = sum(x==3)
}

p1st_nonAsian$gene = rownames(p1st_nonAsian)
p1st_nonAsian$race = "non-Asian"
#hm: only consider TMZ treated cases
p1st_nonAsian$wt[p1st_nonAsian$gene=='Hypermutation'] = 0  #table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_nonAsian$wt[p1st_nonAsian$gene=='MMR'] = 0#table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
#
p1st_nonAsian$alt_ini = (p1st_nonAsian$ini+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)
p1st_nonAsian$alt_rec = (p1st_nonAsian$rec+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)

p1st = rbind(p1st_Asian, p1st_nonAsian)
p1st$inip = p1st$ini/(p1st$wt+p1st$shr+p1st$ini)
p1st$recp = p1st$rec/(p1st$wt+p1st$shr+p1st$rec)
p1st$shrp = p1st$shr/(p1st$wt+p1st$shr+p1st$ini)
p1st$wtp = p1st$wt/(p1st$wt+p1st$shr+p1st$ini)
p1st[is.na(p1st)]=0
p1st_idhcod = p1st
#

tmp = melt(p1st[,c('wtp','inip','recp','shrp','gene','race')], id.vars = c('gene','race'),variable.name = 'status')
tmp$xlab = paste(tmp$gene,tmp$race, sep = '.') #lost during evolution

tmp$value[tmp$race=='non-Asian']=-tmp$value[tmp$race=='non-Asian']
tmp$race = factor(tmp$race, levels = c("non-Asian",'Asian'))
tmp_idhcod = tmp
gnorder = scdf$gene #p1st$gene[order(p1st$alt_rec[p1st$race=='non-Asian'] - p1st$alt_rec[p1st$race=='Asian'],decreasing = T)]
tmp$gene = factor(tmp$gene, levels = gnorder)
#gnorder2 = gnorder
#tmp = tmp[!tmp$gene %in% c('IDH1.2','codel','CIC','FUBP1','METex14','ZM','F3T3','TERTp'),]
#tmp$gene = factor(tmp$gene)

#
df2 = data.frame(gene = rev(levels(tmp$gene)), p_alt = 1, p_evol = 1,p_overall = 1, stringsAsFactors = F)
for (ix in 1:nrow(df2)){ #comapre Asian and non-Asian
  tb = p1st[p1st$gene==df2$gene[ix],1:4]
  tb$alt = tb[,2]+tb[,3]+tb[,4]
  p = fisher.test(tb[,c(1,5)])$p.value 
  df2$p_alt[ix] = signif(p,2)
  tb$statuschange = tb[,2]+tb[,3]
  tb$nochange = tb[,1]+tb[,4]
  q = fisher.test(tb[,c(7,6)])$p.value
  df2$p_evol[ix] = signif(q,2)
  r = fisher.test(tb[,c(1,4,6)])$p.value
  r = fisher.test(tb[,1:4])$p.value
  df2$p_overall[ix] = signif(r,2)
}
df2$gene = factor(df2$gene, levels = levels(tmp$gene))
df2_idhcod=df2
df2 = reshape2::melt(df2, id.vars = c('gene'))

p3<-ggplot(mapping = aes(x = gene, y = 100*value),data = tmp[!is.na(tmp$gene),])+
  geom_bar(mapping = aes(fill = status),stat = 'identity', show.legend = F)+
  geom_hline(yintercept = 0, size= 0.25)+
  scale_fill_manual(values = c(NA,"#d02a7c","#010078","#e2b449"))+
  scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = paste0(sbt,'') )+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0, 0, 0), "cm"))

p3_2<- ggplot(mapping = aes(x = gene, y = 90,label = value),data = df2[df2$variable=='p_overall',])+
  geom_text(size=2,position = 'identity',aes(color = value<0.1), show.legend = F)+
  scale_color_manual(values = c('black','red'))+
  coord_flip()+#facet_grid( ~ variable,scales = 'free_x')+
  theme_bw()+labs(y = '', x = '' )+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),
        panel.grid.minor.x  = element_blank(),panel.grid.major.x  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

cowplot::plot_grid(p3,p3_2,align = 'h',rel_widths = c(0.8,0.2),nrow = 1)
write.table(p1st_idhcod, file = 'p1st_idhcod.txt',quote =F, sep = "\t")
#


cowplot::plot_grid(p1,p1_2,
                   p2+theme(axis.text.y = element_blank()),p2_2,
                   p3+theme(axis.text.y = element_blank()),p3_2,
                   align = 'h',nrow = 1, rel_widths = c(0.35,0.05,0.25,0.05,0.25,0.05))
cowplot::plot_grid(p1, p2+theme(axis.text.y = element_blank()),
                   p3+theme(axis.text.y = element_blank()),
                   align = 'h',nrow = 1, rel_widths = c(0.35,0.25,0.25))
scdf$gene = factor(scdf$gene, levels = scdf$gene)
p4<-ggplot(mapping = aes(x = gene, y = Asian_Cauca_Dif_Founding.Events),data = scdf)+
  geom_bar(stat = 'identity', show.legend = F,fill = 'black')+
  geom_hline(yintercept = 0, size= 0.25)+
  #scale_fill_manual(values = c(NA,"#d02a7c","#010078","#e2b449"))+
  #scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = 'Diff. in founding')+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0, 0, 0), "cm"))

p5<-ggplot(mapping = aes(x = gene, y = Dif_Evo),data = scdf)+
  geom_bar(stat = 'identity', show.legend = F,fill = '#010078')+
  geom_hline(yintercept = 0, size= 0.25)+
  #scale_fill_manual(values = c(NA,"#d02a7c","#010078","#e2b449"))+
  #scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = 'Evolutionary' )+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0, 0, 0), "cm"))

cowplot::plot_grid(p1, p2+theme(axis.text.y = element_blank()),
                   p3+theme(axis.text.y = element_blank()),
                   #p5+theme(axis.text.y = element_blank()),
                   p4+theme(axis.text.y = element_blank()),
                   align = 'h',nrow = 1, rel_widths = c(0.25,0.18,0.18,0.12,0.12))

#
#genomic feature and branched evolution
#
br=read.delim('~/Documents/pangliomaevolution/Rev1/Novelty/TreatmentEffect/glioma.branch.evolution.txt',sep = " ")
ini$br = br$branched[match(ini$Patient_ID, br$ID)]
ini$Subtype = paste(ini$IDH, ini$X1p19qcodel, sep = "-")
sbt = "IDHmut-codel"
fisher.test(table(ini$Grade_1[ini$Subtype==sbt], ini$br[ini$Subtype==sbt]))
fisher.test(table(ini$Grade_2[ini$Subtype==sbt], ini$br[ini$Subtype==sbt]))
fisher.test(table(ini$disease_progression[ini$Subtype==sbt], ini$br[ini$Subtype==sbt]))
wilcox.test(ini$Age[ini$Subtype==sbt]~ini$br[ini$Subtype==sbt])
fisher.test(table(ini$Race[ini$Subtype==sbt],ini$br[ini$Subtype==sbt]))
fisher.test(table(ini$HM_R[ini$Subtype==sbt],ini$br[ini$Subtype==sbt]))

gn = NULL; pv = NULL;subtype = NULL
for (sbt in c("IDHwt-noncodel","IDHmut-noncodel","IDHmut-codel")){
  for ( i in 26:61){
    #print(names(ini)[i])
    if (nrow(table(ini[ini$Subtype==sbt,i],ini$br[ini$Subtype==sbt]))>1){
      tt = fisher.test(table(ini[ini$Subtype==sbt,i],ini$br[ini$Subtype==sbt]))
      if (tt$p.value <= 1){
        gn = append(gn,names(ini)[i] )
        pv = append(pv ,tt$p.value )
        subtype = append(subtype,sbt)
      }
    }
    
  }
  
}
df = data.frame(gn,pv, subtype)
df$subtype = factor(df$subtype, levels = c("IDHwt-noncodel","IDHmut-noncodel","IDHmut-codel"))
df$padj = p.adjust(df$pv)
ggplot(df, aes(x = gn, y = padj, fill = subtype))+
  geom_bar(stat = 'identity',position = position_dodge())+
  coord_flip()+scale_fill_manual(values = c("#107050","#c7d34d","#55d9c0"))+theme_classic()+
  labs(y = "adjusted p-value", x = '')+theme(axis.text = element_text(color = 'black'))


