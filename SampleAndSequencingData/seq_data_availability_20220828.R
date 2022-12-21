dt = read.delim("~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/seq_data_availability_0723.txt",stringsAsFactors = F)
ini = read.delim('~/Dropbox/2020-PanGliomaEvolution/clinicalmutationCNAs/Pairedgliomas.genomics_initial_20220204.txt', na.strings = c('NA','#N/A'))
ini$Subtype = ifelse(ini$IDH=="IDHwt", "IDHwt", ifelse(ini$X1p19qcodel=='codel',"IDHcodel","IDHnoncod"))
dt$Subtype = ini$Subtype[match(dt$Patient_ID, ini$Patient_ID)]

#Add new RNAseq data from GLASS,2022
osnd = read.delim('~/Documents/Research/Longitudinal/GLASS/glass_cell2022_oldsamplenewRNAseq.txt')
dt$RNAseq_Initial[dt$Patient_ID %in% osnd$ID_panglioma] = 'YES'
dt$RNAseq_Recurrence[dt$Patient_ID %in% osnd$ID_panglioma] = 'YES'

#Add new samples from GLASS, 2022
nsnd = read.delim('~/Documents/Research/Longitudinal/GLASS/glass_cell2022_newSamplewithRNAseq.txt')
pl = read.delim('~/Documents/Research/Longitudinal/GLASS/GLASS-DNAseq-samples-hmmissing.hmstatusdetermined.txt')
nsnd$plt = pl$platform[match(nsnd$Patient, pl$Patient)]
nsnd$plt[nsnd$Patient=='GLSS-MD-0011']='WGS'
nsnd$plt[nsnd$Patient=='GLSS-HF-4B46']='WXS'
nsnd = nsnd[order(nsnd$Patient),]
dta = data.frame(Patient_ID = c(paste0('PSA0',1:9),paste0('PSA',10:33)), Cohorts = 'GLSS',Asian = 'no',
                 WGS_normal = c(rep('NO',14),rep('YES',19)),WGS_Initial = c(rep('NO',14),rep('YES',19)), WGS_Recurrence=c(rep('NO',14),rep('YES',19)),
                 WES_Normal = c(rep('YES',14),rep('NO',19)),WES_Initial = c(rep('YES',14),rep('NO',19)), WES_Recurrence=c(rep('YES',14),rep('NO',19)),
                 RNAseq_Initial = rep('YES',33),RNAseq_Recurrence = rep('YES',33),GLASS='yes', Subtype = nsnd$Subtype[!duplicated(nsnd$Patient)])
dta$Subtype[dta$Subtype=='IDHmut-noncodel'] = 'IDHnoncod'
dta$Subtype[dta$Subtype=='IDHmut-codel'] = 'IDHcodel'

dt= rbind(dt, dta)


dt$Subtype = factor(dt$Subtype, levels = c('IDHwt',"IDHnoncod","IDHcodel"))
dt = dt[,c("Cohorts","Asian","GLASS","Subtype","WGS_normal","WGS_Initial","WGS_Recurrence","WES_Normal","WES_Initial","WES_Recurrence","RNAseq_Initial","RNAseq_Recurrence")]
dt$New = ifelse(dt$Cohorts %in% c("Tiantan","SMCn","CUHK"),"ANew",NA)
dt$EAGLE = ifelse(dt$Asian=='yes',"Asian",NA)
dt$GLASS = ifelse(dt$GLASS=='yes',"GLASS",NA)

dt$WGS_normal[dt$WGS_normal=="YES"]="WGS"
dt$WGS_Initial[dt$WGS_Initial=="YES"]="WGS"
dt$WGS_Recurrence[dt$WGS_Recurrence=="YES"]="WGS"

dt$WES_Normal[dt$WES_Normal=="YES"&dt$Cohorts!= "MSKCC"]="WES"
dt$WES_Normal[dt$WES_Normal=="panel"&dt$Cohorts== "MSKCC"]="XPanel"
dt$WES_Initial[dt$WES_Initial=="YES"&dt$Cohorts!= "MSKCC"]="WES"
dt$WES_Initial[dt$WES_Initial=="panel"&dt$Cohorts== "MSKCC"]="XPanel"
dt$WES_Recurrence[dt$WES_Recurrence=="YES"&dt$Cohorts!= "MSKCC"]="WES"
dt$WES_Recurrence[dt$WES_Recurrence=="panel"&dt$Cohorts== "MSKCC"]="XPanel"

dt$RNAseq_Initial[dt$RNAseq_Initial=="YES"]="RNAseq"
dt$RNAseq_Recurrence[dt$RNAseq_Recurrence=="YES"]="RNAseq"


dt[dt=="NO"]=NA
library(ComplexHeatmap)
alter_fun = list(
  background=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "grey90", col = "grey90"))
  },
  Asian=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = "red"))
  },
  ANew=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "black", col = "black"))
  },
  GLASS=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "black", col = "black"))
  },
  IDHwt=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#107050", col = "#107050"))
  },
  IDHnoncod=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#55d9c0", col = "#55d9c0"))
  },
  IDHcodel=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#c7d34d", col = "#c7d34d"))
  },
  WGS=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#B4322C", col = "#B4322C"))
  },
  WES=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "orange", col = "orange"))
  },
  XPanel=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#b2df8a", col = "#b2df8a")) #
  },
  RNAseq=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#5A7DB1", col = "#5A7DB1"))
  }
)

col = c("WGS" = "#B4322C","WES" = "orange","RNAseq" = "#5A7DB1","XPanel"="#b2df8a","ANew"= "black","GLASS"='black',"Asian"= "black","IDHwt"="#107050","IDHnoncod"="#55d9c0", "IDHcodel"="#c7d34d")
dt0 = dt
#dt = dt0[!dt0$Cohorts %in% c("MSKCC","GLSS"),c("New","WES_Normal","WES_Initial","WES_Recurrence","WGS_normal","WGS_Initial","WGS_Recurrence","RNAseq_Initial","RNAseq_Recurrence")]
dt = dt0[,c("Subtype","EAGLE","New",'GLASS',"WES_Recurrence","WES_Initial","WES_Normal","WGS_Recurrence","WGS_Initial","WGS_normal","RNAseq_Recurrence","RNAseq_Initial")]
dt = dt[order(dt$Subtype,dt$New,dt$EAGLE, dt$GLASS),]
oncoPrint(mat = t(dt), pct_gp = gpar(fontsize = 9),get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, remove_empty_columns = F,row_order =  names(dt), column_order = rownames(dt),
          #right_annotation = NULL,
          row_names_side = "left",show_column_names =F,show_pct = T,pct_side = 'right',
          top_annotation = NULL, row_names_gp = gpar(fontsize = 9)
          #column_order = NULL,column_names_gp = gpar(fontsize = 7),
           )


for (i in 1:ncol(dt)) {print(table(dt[,i]))}
#width = 12inch is good

# dt = dt0[which(dt0$New=="New"),c("New","WES_Normal","WES_Initial","WES_Recurrence","WGS_normal","WGS_Initial","WGS_Recurrence","RNAseq_Initial","RNAseq_Recurrence")]
# dt = dt[order(dt$New,dt$WES_Initial,dt$WES_Recurrence,dt$WES_Normal,dt$WGS_Initial,dt$WGS_Recurrence,dt$WGS_normal,dt$RNAseq_Initial,dt$RNAseq_Recurrence),]
# oncoPrint(mat = t(dt), pct_gp = gpar(fontsize = 9),get_type = function(x) strsplit(x, ";")[[1]],
#           alter_fun = alter_fun, col = col, remove_empty_columns = F,row_order = NULL, column_order = NULL,
#           row_names_side = "left",pct_side = "right",#show_column_names =F,#show_row_barplot = T,
#           top_annotation = NULL, row_names_gp = gpar(fontsize = 9)
#           #column_order = NULL,column_names_gp = gpar(fontsize = 7),
# )
# 


