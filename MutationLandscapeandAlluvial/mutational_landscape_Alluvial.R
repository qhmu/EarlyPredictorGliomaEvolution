setwd('~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/')
options(stringsAsFactors = F)

ini = read.delim('~/Dropbox/2020-PanGliomaEvolution/clinicalmutationCNAs/Pairedgliomas.genomics_initial_20220204.txt', na.strings = c('NA','#N/A'))
rec = read.delim('~/Dropbox/2020-PanGliomaEvolution/clinicalmutationCNAs/Pairedgliomas.genomics_recurrence_20220204.txt', na.strings = c('NA','#N/A'))
stopifnot(identical(ini$Patient_ID, rec$Patient_ID))

c710 = read.delim('~/Downloads/panglioma.chr710.new.tsv')
ini$chr7gain10loss = c710$new_chr710_ini[match(ini$Patient_ID, c710$Patient_ID)]
rec$chr7gain10loss = c710$new_chr710_rec[match(rec$Patient_ID, c710$Patient_ID)]
# #only show East Asian
# dt = read.delim("~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/seq_data_availability_0723.txt",stringsAsFactors = F)
#
# ini = ini[ini$Patient_ID %in% dt$Patient_ID[which(dt$Asian=='yes')],]
# rec = rec[rec$Patient_ID %in% dt$Patient_ID[which(dt$Asian=='yes')],]
# ##

gns = c("IDH1.2","TP53","TP53del","ATRX","TERTp","CIC","FUBP1",
        "MDM2amp","MDM4amp","CDKN2Adel","CDK4amp","RB1","RB1del",
        "EGFR","EGFRamp","PDGFRA","PDGFRAamp","METamp","F3T3","ZM",
        "METex14","EGFRvIII","NF1","NF1del",
        "PTEN","PTENdel","PIK3CA","PIK3CG","PIK3R1", "MYC_gain","chr7gain10loss","chr17p_NLOH","codel")#

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

# ini$Grade2021_Ini = as.numeric(substr(ini$WHO2021Classification_initial,
#                                       start = nchar(ini$WHO2021Classification_initial),
#                                       stop = nchar(ini$WHO2021Classification_initial)))
# ini$Grade2021_Rec = as.numeric(substr(ini$WHO2021Classification_recurrence,
#                                       start = nchar(ini$WHO2021Classification_recurrence),
#                                       stop = nchar(ini$WHO2021Classification_recurrence)))
grd21new = read.delim('~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/pairedGlioma.who2021grade.txt')
ini$Grade2021_Ini =grd21new$grade2021_ini[match(ini$Patient_ID, grd21new$Patient_ID)]
ini$Grade2021_Rec =grd21new$grade2021_rec[match(ini$Patient_ID, grd21new$Patient_ID)]

p1a = data.frame(Subtype = paste0(ini$IDH,ini$X1p19qcodel),
                 #Cohort = ifelse(ini$Cohorts %in% c('CGGA','Tiantan','CUHK','SMC'),"New",'Published'),
                 Grade16_I = ini$WHO2016_grade_initial,
                 Grade16_R = ini$WHO2016_grade_recurrence,
                 Grade21_I = ini$Grade2021_Ini,
                 Grade21_R = ini$Grade2021_Rec,
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

lowPurity = c("PS003","PS006","PS011","PS017","PS019","PS028","PS065","PS069","PS095","PS097","PS103",
              "PS104","PS119","PS123","PS269","PS274","PS278","PS284","PS287")


#idx = which(!(ini$Cohorts == "MSKCC"|is.na(ini$IDH1.2)| ini$Cohorts == "GLSS"  |ini$Patient_ID %in% lowPurity )); p1a = p1a[idx,]; p1m = p1m[idx,] #highqual203
idx = which(!(is.na(ini$IDH1.2))  ); p1a = p1a[idx,]; p1m = p1m[idx,] #all DNA
#idx = which((!(is.na(ini$IDH1.2)) & (p1a$New=="New"))); p1a = p1a[idx,]; p1m = p1m[idx,] #all new
#idx = which((!(is.na(ini$IDH1.2)) & (!(is.na(p1a$Race)) & (p1a$Race=="Asian"))  )); p1a = p1a[idx,]; p1m = p1m[idx,] #all Asian
#idx = which(ini$Cohorts=='GLSS'); p1a = p1a[idx,]; p1m = p1m[idx,] #GLASS
#idx = which(startsWith(rownames(p1a),"PSX")); p1a = p1a[idx,]; p1m = p1m[idx,] #highqual203



p1m = as.data.frame(t(p1m))
idx = order(p1a$Subtype,p1a$Grade21_R,p1a$Grade16_R, p1a$Grade21_I,p1a$Grade16_I, decreasing = T)
p1a = p1a[idx,]
p1m = p1m[,idx]
p1mmyc = as.data.frame(t(p1m)); p1m = p1m#[rownames(p1m)!="MYC_gain",]

##################### ############# ######
#                     ************  ******
#     1 mutational     landscape    \\\\\\
##@@@########@@##### @@@@@@@@@@$@   #  $ #
#################### ############# #######

cls = c("#d02a7c","#010078","#e2b449","#000000","coral4")
value2color = function(x){ #
  if (is.na(x)){
    cl = "#DDDDDD"
  } else if (x ==0){
    cl = "white"
  } else if (x ==1){
    cl = cls[1]
  } else if (x ==2){
    cl = cls[2]
  } else if (x==3){
    cl = cls[3]
  }else{
    cl = cls[4]
  }
  return(cl)
}

#set overall margins and initilize plot
par(mar = c(2.5,0.5,0.5,0.1), fig = c(0,1,0,1),las=2,mgp = c(3,0.5,0))
offset =2;  #offset is dist of the top to 0, i.e. the top is at -2, so there is space for some title/label above the top
npanels = 7 #the number of chunks in row
w = sum(p1a$Subtype=="IDHwt"); n = sum(p1a$Subtype=="IDHnon"); d = sum(p1a$Subtype=="IDHcodel");o = nrow(p1a);
npatient = ncol(p1m); nfeature= nrow(p1m) + ncol(p1a)  +0.5*npanels+ offset
pxi = ifelse(npatient >200,1,0.85) #width of each cell. max 1, <1 values will lead to small gaps between cells. recommend to use 0.85 when number of samples <=200, and 1 if >=250
pxj = ifelse(nfeature >50,1,0.9) #height of each cell. max 1, <1 values will lead to small gaps between row. recommend to use 0.9 when number of features <50, and 1 if >50

gapri = max(1,as.integer(0.01*npatient)) #gap width in row, i.e.e the gap between different subtypes
pdf('landscape_20220616.pdf',width = 8.5, height = 6.5)
plot(-npatient, -nfeature,  xlim = c(-(npatient + 2*gapri +10),0), ylim = c(-nfeature,0), type = "n",  bty = "n",xaxt = "n", yaxt = "n",  ann=FALSE,frame.plot = T,yaxs="i") #+2 for spaces, +10 for labels
#idx = which(rownames(p1a)%in%c("PS157","PS137","PS146","PS138","PS135","PS143","PS140","PS132"))
#par(las=2);axis(side = 1, at = -(npatient+2*gapri)+idx,labels = rownames(p1a)[idx] )
#abline(v = -(npatient+2*gapri)+idx)
#text(x = -npatient + w/2, -1, labels = paste0("IDHwt (n = ",w,")"), col = '#107050',cex=1,font=2)
#text(x = -npatient + w + n/2, -1, labels = paste0("IDHnon (n = ",n,")"), col = '#55d9c0',cex=1,font=2)
#text(x = -npatient + w + n+ d/2, -1, labels = paste0("IDHcodel (n = ",d,")"), col = '#c7d34d',cex=1,font=2)
#Map(axis, side=1, at=c(-224:-118,-116:-35,-33:-1)+0.5, col.axis=as.factor(p1$cohort[idx]), labels=rownames(p1a), lwd=0, las=2,cex.axis = 0.45)
#axis(side = 1,at = c(-224:-118,-116:-35,-33:-1)+0.5,labels = F,tck=-0.01)
#chr level alterations, 4 rows.
n1=4  #chr7p chr7q chr10p chr10q
for (i in 1:w){
  for (j in 1:n1){
    rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -(nfeature+0.5) + j, ytop = -(nfeature+0.5) + j + 1 -(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
    # v=p1m[(nrow(p1m)+1-j),i]
    # if (is.na(v)){}
  }
}
for (i in (w+1):(w+n)){
  for (j in 1:n1){
    rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -(nfeature+0.5) + j, ytop = -(nfeature+0.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+n+1):npatient){
  for (j in 1:n1){
    rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -(nfeature+0.5) + j, ytop = -(nfeature+0.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
lbs = c("MYC gain","+7/-10","17p_CNLOH","1p/19q codel")
lbs = rev(lbs)
for (j in 1:n1){
  text(-(npatient + 2 +10),-(nfeature+0.5) + j + 0.5, labels = lbs[j], cex = 0.6, col = ifelse(j>0,"black","black"), adj=1) #rownames(p1m)[nrow(p1m)+1-j]
}
#pi3k pathway. five rows.
n2=5
for (i in 1:w){
  for (j in (n1+1):(n1+n2)){
    rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -nfeature + j, ytop = -nfeature + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+1):(w+n)){
  for (j in (n1+1):(n1+n2)){
    rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -nfeature + j, ytop = -nfeature + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+n+1):npatient){
  for (j in (n1+1):(n1+n2)){
    rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -nfeature + j, ytop = -nfeature + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
lbs = c(expression(italic("PTEN")),substitute(paste(italic('PTEN'), " del")),expression(italic("PIK3CA")),expression(italic("PIK3CG")),expression(italic("PIK3R1")))
lbs = rev(lbs)
for (j in 1:n2){
  text(-(npatient + 2 +10),-nfeature + n1+j + 0.5, labels =lbs[j], cex = 0.6, col = ifelse(j>0,"black","black"), adj=1) #add label
}
#rtk pathway. 11 rows
n3= 11
for (i in 1:w){
  for (j in (n1+n2+1):(n1+n2+n3)){
    rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -(nfeature-0.5) + j, ytop = -(nfeature-0.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+1):(w+n)){
  for (j in (n1+n2+1):(n1+n2+n3)){
    rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -(nfeature-0.5) + j, ytop = -(nfeature-0.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+n+1):npatient){
  for (j in (n1+n2+1):(n1+n2+n3)){
    rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -(nfeature-0.5) + j, ytop = -(nfeature-0.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
lbs = c(expression(italic("EGFR")),substitute(paste(italic('EGFR'), " amp")),expression(italic("PDGFRA")),substitute(paste(italic('PDGFRA'), " amp")),substitute(paste(italic('MET'), " amp")),
        expression(italic("FGFR3-TACC3")),expression(italic("PTPRZ1-MET")),substitute(paste(italic('MET'), "ex14")),substitute(paste(italic('EGFR'), "vIII")),expression(italic("NF1")),substitute(paste(italic('NF1'), " del")))
lbs = rev(lbs)
for (j in 1:n3){
  text(-(npatient + 2 +10),-(nfeature-0.5) + n1+n2+j + 0.5, labels = lbs[j], cex = 0.6, col = "black", adj=1)
}
#cell cycle. 6 rows
n4=6
for (i in 1:w){
  for (j in (n1+n2+n3+1):(n1+n2+n3+n4)){
    rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -(nfeature-1) + j, ytop = -(nfeature-1) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+1):(w+n)){
  for (j in (n1+n2+n3+1):(n1+n2+n3+n4)){
    rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -(nfeature-1) + j, ytop = -(nfeature-1) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+n+1):npatient){
  for (j in (n1+n2+n3+1):(n1+n2+n3+n4)){
    rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -(nfeature-1) + j, ytop = -(nfeature-1) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
lbs = c(substitute(paste(italic('MDM2'), " amp")),substitute(paste(italic('MDM4'), " amp")),substitute(paste(italic('CDKN2A'), " del")),
        substitute(paste(italic('CDK4'), " amp")),expression(italic("RB1")),substitute(paste(italic('RB1'), " del")))#c("MDM2amp","MDM4amp","CDKN2Adel","CDK4amp","RB1","RB1del")
lbs = rev(lbs)
for (j in 1:n4){
  text(-(npatient + 2 +10),-(nfeature-1) + n1+n2+n3+j + 0.5, labels = lbs[j], cex = 0.6, col = "black", adj=1)
}

#key markers, n =9
n5=7
for (i in 1:w){
  for (j in (n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)){
    rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -(nfeature-1.5) + j, ytop = -(nfeature-1.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+1):(w+n)){
  for (j in (n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)){
    rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i +pxi, ybottom = -(nfeature-1.5) + j, ytop = -(nfeature-1.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
for (i in (w+n+1):npatient){
  for (j in (n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)){
    rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -(nfeature-1.5) + j, ytop = -(nfeature-1.5) + j + 1-(1-pxj), col = value2color(p1m[(nrow(p1m)+1-j),i]),border=NA)
  }
}
lbs = c(expression(italic('IDH1')),expression(italic('TP53')),substitute(paste(italic('TP53'), " del")),
        expression(italic('ATRX')),substitute(paste(italic('TERT'), " promoter")),expression(italic('CIC')),expression(italic('FUBP1')))
lbs = rev(lbs)
for (j in 1:n5){
  text(-(npatient + 2 +10),-(nfeature-1.5) + n1+n2+n3+n4+j + 0.5, labels = lbs[j], cex = 0.6, col = "black", adj=1)
}

#hypermutation, n=5
n6=4
p1ar = p1a[,c("MMR", "Hypermutation","MGMTfusion","TMZ_R")]
for (i in 1:w){
  for (j in (n1+n2+n3+n4+n5+1):(n1+n2+n3+n4+n5+n6)){
    rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -(nfeature-2) + j, ytop = -(nfeature-2) + j + 1-(1-pxj), col = value2color(p1ar[i,j-(n1+n2+n3+n4+n5)]),border=NA)
  }
}
for (i in (w+1):(w+n)){
  for (j in (n1+n2+n3+n4+n5+1):(n1+n2+n3+n4+n5+n6)){
    rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -(nfeature-2) + j, ytop = -(nfeature-2) + j + 1-(1-pxj), col = value2color(p1ar[i,j-(n1+n2+n3+n4+n5)]),border=NA)
  }
}
for (i in (w+n+1):npatient){
  for (j in (n1+n2+n3+n4+n5+1):(n1+n2+n3+n4+n5+n6)){
    rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -(nfeature-2) + j, ytop = -(nfeature-2) + j + 1-(1-pxj), col = value2color(p1ar[i,j-(n1+n2+n3+n4+n5)]),border=NA)
  }
}
lbs = c("Alkylator_Rec","MGMT fusion","Hypermutation","MMR")
lbs = rev(lbs)
for (j in 1:n6){
  text(-(npatient + 2 +10),-(nfeature-2) + n1+n2+n3+n4+n5+j + 0.5, labels = lbs[j], cex = 0.6, col = "black", adj=1)
}
#clinical annotations
# for (i in 1:86){
#   text(-87 + i,-35, labels = names(p1m)[i], cex = 0.5, srt = 90)
# }
n7=5;
for (i in 1:w){
  rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -offset-1, ytop = -offset-(1-pxj),
       col = ifelse(p1a$Subtype[i] == "IDHwt","#107050", ifelse(p1a$Subtype[i] == "IDHcodel","#c7d34d","#55d9c0")),border=NA)
}
for (i in (w+1):(w+n)){
  rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -offset-1, ytop = -offset-(1-pxj),
       col = ifelse(p1a$Subtype[i] == "IDHwt","#107050", ifelse(p1a$Subtype[i] == "IDHcodel","#c7d34d","#55d9c0")),border=NA)
}
for (i in (w+n+1):npatient){
  rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -offset-1, ytop = -offset-(1-pxj),
       col = ifelse(p1a$Subtype[i] == "IDHwt","#107050", ifelse(p1a$Subtype[i] == "IDHcodel","#c7d34d","#55d9c0")),border=NA)
}
text(-(npatient + 2 +10), -2.55, labels = "Subtype", cex = 0.6, adj=1)
for (i in 1:w){
  rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -offset-2, ytop = -offset-1-(1-pxj),
       col = ifelse(is.na(p1a$Grade16_I[i]),"#DDDDDD",ifelse(p1a$Grade16_I[i] == "II", "#c6dbef", ifelse(p1a$Grade16_I[i] == "III", "#6baed6","#2171b5"))),
       border=NA)

}
for (i in (w+1):(w+n)){
  rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -offset-2, ytop = -offset-1-(1-pxj),
       col = ifelse(is.na(p1a$Grade16_I[i]),"#DDDDDD",ifelse(p1a$Grade16_I[i] == "II", "#c6dbef", ifelse(p1a$Grade16_I[i] == "III", "#6baed6","#2171b5"))),
       border=NA)
}
for (i in (w+n+1):npatient){
  rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -offset-2, ytop = -offset-1-(1-pxj),
       col = ifelse(is.na(p1a$Grade16_I[i]),"#DDDDDD",ifelse(p1a$Grade16_I[i] == "II", "#c6dbef", ifelse(p1a$Grade16_I[i] == "III", "#6baed6","#2171b5"))),
       border=NA)
}
text(-(npatient + 2 +10), -3.55, labels = "WHO16Grade_Ini", cex = 0.6, adj=1)
for (i in 1:w){
  rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -offset-3, ytop = -offset-2-(1-pxj),
       col = ifelse(is.na(p1a$Grade16_R[i]),"#DDDDDD",ifelse(p1a$Grade16_R[i] == "II", "#c6dbef", ifelse(p1a$Grade16_R[i] == "III", "#6baed6","#2171b5"))),
       border=NA)
}
for (i in (w+1):(w+n)){
  rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -offset-3, ytop = -offset-2-(1-pxj),
       col = ifelse(is.na(p1a$Grade16_R[i]),"#DDDDDD",ifelse(p1a$Grade16_R[i] == "II", "#c6dbef", ifelse(p1a$Grade16_R[i] == "III", "#6baed6","#2171b5"))),
       border=NA)
}
for (i in (w+n+1):npatient){
  rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -offset-3, ytop = -offset-2-(1-pxj),
       col = ifelse(is.na(p1a$Grade16_R[i]),"#DDDDDD",ifelse(p1a$Grade16_R[i] == "II", "#c6dbef", ifelse(p1a$Grade16_R[i] == "III", "#6baed6","#2171b5"))),
       border=NA)
}
text(-(npatient + 2 +10), -4.55, labels = "WHO16Grade_Rec", cex = 0.6, adj=1)

for (i in 1:w){
  rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -offset-4, ytop = -offset-3-(1-pxj),
       col = ifelse(is.na(p1a$Grade21_I[i]),"#DDDDDD",ifelse(p1a$Grade21_I[i] == 2, "#c6dbef", ifelse(p1a$Grade21_I[i] == 3, "#6baed6","#2171b5"))),
       border=NA)

}
for (i in (w+1):(w+n)){
  rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -offset-4, ytop = -offset-3-(1-pxj),
       col = ifelse(is.na(p1a$Grade21_I[i]),"#DDDDDD",ifelse(p1a$Grade21_I[i] == 2, "#c6dbef", ifelse(p1a$Grade21_I[i] == 3, "#6baed6","#2171b5"))),
       border=NA)
}
for (i in (w+n+1):npatient){
  rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -offset-4, ytop = -offset-3-(1-pxj),
       col = ifelse(is.na(p1a$Grade21_I[i]),"#DDDDDD",ifelse(p1a$Grade21_I[i] == 2, "#c6dbef", ifelse(p1a$Grade21_I[i] == 3, "#6baed6","#2171b5"))),
       border=NA)
}
text(-(npatient + 2 +10), -5.55, labels = "WHO21Grade_Ini", cex = 0.6, adj=1)
for (i in 1:w){
  rect(xleft = -(npatient+2*gapri+1)+ i, xright = -(npatient+2*gapri+1) + i + pxi, ybottom = -offset-5, ytop = -offset-4-(1-pxj),
       col = ifelse(is.na(p1a$Grade21_R[i]),"#DDDDDD",ifelse(p1a$Grade21_R[i] == 2, "#c6dbef", ifelse(p1a$Grade21_R[i] == 3, "#6baed6","#2171b5"))),
       border=NA)
}
for (i in (w+1):(w+n)){
  rect(xleft = -(npatient+1*gapri+1)+ i, xright = -(npatient+1*gapri+1) + i + pxi, ybottom = -offset-5, ytop = -offset-4-(1-pxj),
       col = ifelse(is.na(p1a$Grade21_R[i]),"#DDDDDD",ifelse(p1a$Grade21_R[i] == 2, "#c6dbef", ifelse(p1a$Grade21_R[i] == 3, "#6baed6","#2171b5"))),
       border=NA)
}
for (i in (w+n+1):npatient){
  rect(xleft = -(npatient+1)+ i, xright = -(npatient+1) + i + pxi, ybottom = -offset-5, ytop = -offset-4-(1-pxj),
       col = ifelse(is.na(p1a$Grade21_R[i]),"#DDDDDD",ifelse(p1a$Grade21_R[i] == 2, "#c6dbef", ifelse(p1a$Grade21_R[i] == 3, "#6baed6","#2171b5"))),
       border=NA)
}
text(-(npatient + 2 +10), -6.55, labels = "WHO21Grade_Rec", cex = 0.6, adj=1)
#now add the frames outside

rect(xleft = -(npatient+2*gapri+1)+1,xright = -(npatient+2*gapri+1)+w+pxi,ybottom = -(n7+offset),ytop = -offset-(1-pxj),lwd=0.7) #clinical
rect(xleft = -(npatient+1*gapri+1)+w+1,xright = -(npatient+1*gapri+1)+w+n+pxi,ybottom = -(n7+offset),ytop = -offset-(1-pxj),lwd=0.7)
rect(xleft = -(npatient+1)+w+n+1,xright = -1+pxi,ybottom = -(n7+offset),ytop = -offset-(1-pxj),lwd=0.7)

n7=5
rect(xleft = -(npatient+2*gapri+1)+1,xright = -(npatient+2*gapri+1)+w+pxi,ybottom = -(n6+n7+offset+0.5),ytop = -(n7+offset+0.5)-(1-pxj),lwd=0.7) #TMZ resistance
rect(xleft = -(npatient+1*gapri+1)+w+1,xright = -(npatient+1*gapri+1)+w+n+pxi,ybottom = -(n6+n7+offset+0.5),ytop = -(n7+offset+0.5)-(1-pxj),lwd=0.7)
rect(xleft = -(npatient+1)+w+n+1,xright = -1+pxi,ybottom = -(n6+n7+offset+0.5),ytop = -(n7+offset+0.5)-(1-pxj),lwd=0.7)

rect(xleft = -(npatient+2*gapri+1)+1,xright = -(npatient+2*gapri+1)+w+pxi,ybottom = -(n5+n6+n7+offset+1),ytop = -(n6+n7+offset+1)-(1-pxj),lwd=0.7) #IDH TERT
rect(xleft = -(npatient+1*gapri+1)+w+1,xright = -(npatient+1*gapri+1)+w+n+pxi,ybottom = -(n5+n6+n7+offset+1),ytop = -(n6+n7+offset+1)-(1-pxj),lwd=0.7)
rect(xleft = -(npatient+1)+w+n+1,xright = -1+pxi,ybottom = -(n5+n6+n7+offset+1),ytop = -(n6+n7+offset+1)-(1-pxj),lwd=0.7)

rect(xleft = -(npatient+2*gapri+1)+1,xright = -(npatient+2*gapri+1)+w+pxi,ybottom = -(n4+n5+n6+n7+offset+1.5),ytop = -(n5+n6+n7+offset+1.5)-(1-pxj),lwd=0.7) #IDH TERT
rect(xleft = -(npatient+1*gapri+1)+w+1,xright = -(npatient+1*gapri+1)+w+n+pxi,ybottom = -(n4+n5+n6+n7+offset+1.5),ytop = -(n5+n6+n7+offset+1.5)-(1-pxj),lwd=0.7)
rect(xleft = -(npatient+1)+w+n+1,xright = -1+pxi,ybottom = -(n4+n5+n6+n7+offset+1.5),ytop = -(n5+n6+n7+offset+1.5)-(1-pxj),lwd=0.7)

rect(xleft = -(npatient+2*gapri+1)+1,xright = -(npatient+2*gapri+1)+w+pxi,ybottom = -(n3+n4+n5+n6+n7+offset+2),ytop = -(n4+n5+n6+n7+offset+2)-(1-pxj),lwd=0.7) #MAPK
rect(xleft = -(npatient+1*gapri+1)+w+1,xright = -(npatient+1*gapri+1)+w+n+pxi,ybottom = -(n3+n4+n5+n6+n7+offset+2),ytop = -(n4+n5+n6+n7+offset+2)-(1-pxj),lwd=0.7)
rect(xleft = -(npatient+1)+w+n+1,xright = -1+pxi,ybottom = -(n3+n4+n5+n6+n7+offset+2),ytop = -(n4+n5+n6+n7+offset+2)-(1-pxj),lwd=0.7)

rect(xleft = -(npatient+2*gapri+1)+1,xright = -(npatient+2*gapri+1)+w+pxi,ybottom = -(n2+n3+n4+n5+n6+n7+offset+2.5),ytop = -(n3+n4+n5+n6+n7+offset+2.5)-(1-pxj),lwd=0.7) #PI3K
rect(xleft = -(npatient+1*gapri+1)+w+1,xright = -(npatient+1*gapri+1)+w+n+pxi,ybottom = -(n2+n3+n4+n5+n6+n7+offset+2.5),ytop = -(n3+n4+n5+n6+n7+offset+2.5)-(1-pxj),lwd=0.7)
rect(xleft = -(npatient+1)+w+n+1,xright = -1+pxi,ybottom = -(n2+n3+n4+n5+n6+n7+offset+2.5),ytop = -(n3+n4+n5+n6+n7+offset+2.5)-(1-pxj),lwd=0.7)

rect(xleft = -(npatient+2*gapri+1)+1,xright = -(npatient+2*gapri+1)+w+pxi,ybottom = -(n1+n2+n3+n4+n5+n6+n7+offset+3),ytop = -(n2+n3+n4+n5+n6+n7+offset+3)-(1-pxj),lwd=0.7) #arms
rect(xleft = -(npatient+1*gapri+1)+w+1,xright = -(npatient+1*gapri+1)+w+n+pxi,ybottom = -(n1+n2+n3+n4+n5+n6+n7+offset+3),ytop = -(n2+n3+n4+n5+n6+n7+offset+3)-(1-pxj),lwd=0.7)
rect(xleft = -(npatient+1)+w+n+1,xright = -1+pxi,ybottom = -(n1+n2+n3+n4+n5+n6+n7+offset+3),ytop = -(n2+n3+n4+n5+n6+n7+offset+3)-(1-pxj),lwd=0.7)

dev.off()
print('Finished ploting mutational landscape!')




###
library(alluvial)
library(dplyr)
f = data.frame(Initial = p1a$Grade21_I, Recurrence = p1a$Grade21_R, Subtype = p1a$Subtype)
f %>% group_by(Subtype,Initial, Recurrence) %>% summarise(n = n()) ->dt
dt
dt = dt[!(is.na(dt$Initial)|is.na(dt$Recurrence)),]
dt$prog = ifelse(dt$Recurrence==4,'yes','no')
dt$prog[which(dt$Subtype=='IDHcodel' & dt$Initial==2 & dt$Recurrence==3)] = 'yes'
#"IDHwt","#107050", ifelse(p1a$Subtype[i] == "IDHcodel","#c7d34d","#55d9c0"
alluvial(dt[,1:3], freq=dt$n,gap.width = .5,
         #col = ifelse(dt$Recurrence == 2, "#c6dbef",ifelse(dt$Recurrence== 3, "#6baed6","#2171b5")),
         #col = ifelse( dt$prog=='yes', "#010078", "#fed976" ),
         col = ifelse(dt$Subtype=="IDHwt","#107050", ifelse(dt$Subtype== "IDHcodel","#c7d34d","#55d9c0")),
         border = NA,alpha=0.75,hide = dt$n<2,
         # ordering = list(
         #   #order(dt$Grade_I,dt$Grade_R,decreasing = T),
         #   NULL,NULL,
         #   NULL,NULL
         #   #order(dt$Grade_I,dt$Grade_R, decreasing = T)
         # ),
         blocks = T,cex = 0.6,cw = 0.2, cex.axis = 0.6)

ini$HM_I[which(ini$HM_I=='NO')] = "No"
f = data.frame(Initial = ini$HM_I, Recurrence = ini$HM_R)
f %>% group_by(Initial, Recurrence) %>% summarise(n = n()) ->dt
dt = dt[!is.na(dt$Initial),]
dt
alluvial(dt[,1:2], freq=dt$n,gap.width = .5,
         col = ifelse(dt$Recurrence == "Yes", "#010078", "#fed976"),

         border = NA,alpha=0.8,hide = dt$n<2,
         # ordering = list(
         #   #order(dt$Grade_I,dt$Grade_R,decreasing = T),
         #   NULL,NULL,
         #   NULL,NULL
         #   #order(dt$Grade_I,dt$Grade_R, decreasing = T)
         # ),
         blocks = T,cex = 0.01,cw = 0.2, cex.axis = 0.6)
