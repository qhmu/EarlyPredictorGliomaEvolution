
setwd('~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/TEDG/')
#1. Define function
mutDirectedGraph <- function(mutation_gene_table){
  library(igraph)
  library(reshape2)
  
  input.table <- mutation_gene_table
  selected_geneList <- as.character(colnames(input.table))
  
  #calculating TEDG edge table
  temp <- rep(0,length(selected_geneList))
  edge.matrix <- temp
  for( i in 2:length(selected_geneList)){
    edge.matrix <- cbind(edge.matrix,temp)
  }
  
  edge.table <- c('geneA','geneB','weight','label')
  
  end <- length(selected_geneList)-1
  for( i in 1:end){
    start <- i+1
    for( j in start:length(selected_geneList) ){
      
      edge.matrix[i,j] <- length( which( (input.table[,i] == 'C') & (input.table[,j] %in% c('P','R')) ) )
      edge.matrix[j,i] <- length( which( (input.table[,j] == 'C') & (input.table[,i] %in% c('P','R')) ) )
      labelA <- paste( rownames(input.table)[which( (input.table[,i] == 'C') & (input.table[,j] %in% c('P','R')) )], collapse = ";")
      labelB <- paste( rownames(input.table)[which( (input.table[,j] == 'C') & (input.table[,i] %in% c('P','R')) )], collapse = ";")
      
      if(edge.matrix[i,j] < edge.matrix[j,i]){
        edge.matrix[i,j] <- 0
        edge <- c(selected_geneList[j],selected_geneList[i],edge.matrix[j,i],labelB)
        edge.table <- rbind(edge.table,edge)
      }
      else if(edge.matrix[i,j] > edge.matrix[j,i]){
        edge.matrix[j,i] <- 0
        edge <- c(selected_geneList[i],selected_geneList[j],edge.matrix[i,j],labelA)
        edge.table <- rbind(edge.table,edge)
      }
      else{
        edge.matrix[i,j] <- 0
        edge.matrix[j,i] <- 0
      }
      
    }
  }
  
  rownames(edge.matrix) <- selected_geneList
  colnames(edge.matrix) <- selected_geneList
  
  edge.table <- edge.table[-1,]
  colnames(edge.table) <- c('geneA','geneB','weight','label')
  rownames(edge.table) <- c(1:nrow(edge.table))
  
  #write.table(edge.table,"TEDGedge.txt",row.names = F,quote = F,sep = '\t')
  
  #calculating TEDG node table
  Mut.freq <- cbind(rep(0,length(selected_geneList)),rep(0,length(selected_geneList)),rep(0,length(selected_geneList)),rep(0,length(selected_geneList)))
  for(i in 1:length(selected_geneList)){
    Mut.freq[i,1] <- length(which(input.table[,i] == 'P'))
    Mut.freq[i,2] <- length(which(input.table[,i] == 'R'))
    Mut.freq[i,3] <- length(which(input.table[,i] == 'C')) * 2
    Mut.freq[i,4] <- length(which(input.table[,i] == 'N')) * 2
  }
  
  sample.size <- rep(0,length(selected_geneList))
  for(i in 1:length(selected_geneList)){
    sample.size[i] <- (Mut.freq[i,1] + Mut.freq[i,2] + Mut.freq[i,3])/(Mut.freq[i,1] + Mut.freq[i,2] + Mut.freq[i,3]+ Mut.freq[i,4])
  }
  
  ins <- rep(0, length(selected_geneList))
  outs <- rep(0,length(selected_geneList))
  for(i in 1:length(selected_geneList)){
    ins[i] <- length(which(edge.matrix[,i]>0))
    outs[i] <- length(which(edge.matrix[i,]>0))
  }
  
  
  pcdf <- rep(1,length(selected_geneList))
  for(i in 1:length(pcdf)){
    if(ins[i] < outs[i]){
      #y = binocdf(x,N,p) computes (x,y) that follow binomial dist (N,p)
      pcdf[i] <- binom.test(ins[i], ins[i]+outs[i], 0.5)$p.value
    }
    else{
      pcdf[i] <- binom.test(outs[i], ins[i]+outs[i], 0.5)$p.value
    }
  }
  
  fc = log2((outs+1) / (ins+1)) # positive = early; negative = late.
  
  node.table <- data.frame(selected_geneList,pcdf,fc,sample.size)
  colnames(node.table) <- c('Gene','P_CDF',	'FC',	'Occurrence')
  #write.table(node.table,"TEDGnode.txt",row.names = F,quote = F,sep = '\t')
  colnames(edge.table) <- c('source','target','weight','label')
  
  node <- data.frame(node.table);edge <- data.frame(edge.table)
  net <- graph_from_data_frame(d=edge[which(as.numeric(edge$weight) > 0),], vertices=node, directed=T)
  
  #network deconvolution
  E(net)$weight2 = as.numeric(E(net)$weight)
  A0 = as_adjacency_matrix(net,attr =  "weight2",sparse = T)
  rho = 1
  A = as.matrix(A0 )
  B = rho*A+diag(dim(A)[1])
  Bi = solve(B)
  C = A%*%Bi
  
  C[C<0]=0;  diag(C) = 0
  eg2 <- melt(C)
  eg2 = eg2[eg2$value>0,]
  names(eg2)=c('source','target','weight')
  edge.table = as.data.frame(edge.table)
  
  eg2$weight = edge.table$weight[match(paste(eg2$source,eg2$target),paste(edge.table$source, edge.table$target))]
  eg2$weight[is.na(eg2$weight)]=0
  eg2 = eg2[eg2$weight>2,]
  #return values
  returnList <- list("edge.table" = edge.table,'edge.deconv' = eg2, "node.table" = node.table,"network"=net)
  return(returnList)
}

#example Data Format:
data.frame(gene1 = c('N','C','P','R'),gene2 = c('C','R','N','P'), row.names = c('Pt1','Pt2','Pt3','Pt4'))
#2. read in data
gi = read.delim('~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/panglioma.initial.clinicalanddrivers.0320.txt',na.strings = c("NA","#N/A"))
gr = read.delim('~/Dropbox/2020-PanGliomaEvolution/Figures/Figure_1/panglioma.recurrence.clinicalanddrivers.0320.txt',na.strings = c("NA","#N/A"))
gi = gi [order(gi$Patient_ID),]; gr = gr[order(gr$Patient_ID),]
stopifnot(identical(gi$Patient_ID, gr$Patient_ID))
gi$HM = 0
gr$HM = ifelse(is.na(gr$HM_R),NA,ifelse(gi$HM_R=='Yes',1,0))

c710 = read.delim('~/Downloads/panglioma.chr710.new.tsv')
gi$chr710 = c710$new_chr710_ini[match(gi$Patient_ID, c710$Patient_ID)]
gr$chr710 = c710$new_chr710_rec[match(gi$Patient_ID, c710$Patient_ID)]

#3. calculate IDHwt
gns_wt = c('chr710',"MYC_gain","EGFRamp","CDKN2Adel",
           "TERTp","TP53","PTEN","PIK3CA","PIK3R1",
           "EGFR","RB1","NF1","EGFRvIII","HM")
giwt = gi[gi$X=="IDHwt_noncodel",c('Patient_ID',gns_wt)]
grwt = gr[gi$X=="IDHwt_noncodel",c('Patient_ID',gns_wt)]

gidhwt = giwt[,-1] + 2*grwt[,-1]
gidhwt_mat = gidhwt
for (i in 1:nrow(gidhwt)){
  for (j in 1:ncol(gidhwt)){
    x = gidhwt[i,j]
    gidhwt_mat[i,j]=ifelse(is.na(x),NA,ifelse(x==0,'N',ifelse(x==1,'P',ifelse(x==2,'R','C'))))
  }
}
idhwt_tedg = mutDirectedGraph(gidhwt_mat)

plot(idhwt_tedg$network)
write.table(idhwt_tedg$edge.table,"IDHwt.202206082.TEDGedge.txt",row.names = F,quote = F,sep = '\t')
write.table(idhwt_tedg$edge.deconv,"IDHwt.202206082.deconTEDGedge.txt",row.names = F,quote = F,sep = '\t')
write.table(idhwt_tedg$node.table,"IDHwt.202206082.TEDGnode.txt",row.names = F,quote = F,sep = '\t')

#4. calculate IDHnoncodel
gns_non = c('IDH1.2',"chr17p_NLOH","TP53","ATRX","PTEN",
            "MYC_gain","CDKN2Adel",
            "METex14","RB1","NF1","PIK3CA","HM")
ginon = gi[gi$X=="IDHmut_noncodel",gns_non]
grnon = gr[gi$X=="IDHmut_noncodel",gns_non]

gidhnon = ginon + 2*grnon
gidhnon_mat = gidhnon
for (i in 1:nrow(gidhnon)){
  for (j in 1:ncol(gidhnon)){
    x = gidhnon[i,j]
    gidhnon_mat[i,j]=ifelse(is.na(x),NA,ifelse(x==0,'N',ifelse(x==1,'P',ifelse(x==2,'R','C'))))
  }
}
idhnon_tedg = mutDirectedGraph(gidhnon_mat)
write.table(idhnon_tedg$edge.table,"IDHnon.202206082.TEDGedge.txt",row.names = F,quote = F,sep = '\t')
write.table(idhnon_tedg$edge.deconv,"IDHnon.202206082.deconTEDGedge.txt",row.names = F,quote = F,sep = '\t')
write.table(idhnon_tedg$node.table,"IDHnon.202206082.TEDGnode.txt",row.names = F,quote = F,sep = '\t')

#5. calculate IDHcod
gns_cod = c('IDH1.2',"codel","MYC_gain","TERTp","CIC","FUBP1","PIK3CA","TP53","PIK3R1","HM")
gicod = gi[gi$X=="IDHmut_codel",gns_cod]
grcod = gr[gi$X=="IDHmut_codel",gns_cod]

gidhcod = gicod + 2*grcod
gidhcod_mat = gidhcod
for (i in 1:nrow(gidhcod)){
  for (j in 1:ncol(gidhcod)){
    x = gidhcod[i,j]
    gidhcod_mat[i,j]=ifelse(is.na(x),NA,ifelse(x==0,'N',ifelse(x==1,'P',ifelse(x==2,'R','C'))))
  }
}
idhcod_tedg = mutDirectedGraph(gidhcod_mat)
write.table(idhcod_tedg$edge.table,"IDHcod.20220608.TEDGedge.txt",row.names = F,quote = F,sep = '\t')
write.table(idhcod_tedg$edge.deconv,"IDHcod.20220608.deconTEDGedge.txt",row.names = F,quote = F,sep = '\t')
write.table(idhcod_tedg$node.table,"IDHcod.20220608.TEDGnode.txt",row.names = F,quote = F,sep = '\t')


##
library(igraph)

E(idhcod_tedg$network)$weight2 = as.numeric(E(idhcod_tedg$network)$weight)
A0 = as_adjacency_matrix(idhcod_tedg$network,attr =  "weight2",sparse = T)
#
rho = 1
A = as.matrix(A0 )
B = rho*A+diag(dim(A)[1])
Bi = solve(B)
C0 = A%*%Bi
# C

A[A<0]=0; A[A>0]=1
g1 = graph_from_adjacency_matrix(A)
plot(g1,layout= layout_on_grid,main = paste0('Rho = ',rho,', before'))

library(reshape2)
DF <- melt(C0)
DF = DF[DF$value>0,]
C=C0; C[C<0]=0; C[C>0.5]=1; diag(C) = 0
g2 =  graph_from_adjacency_matrix(C)
plot(g2, layout= layout_on_grid, main = paste0('Rho = ',rho,', after'))

eg2 = as.data.frame(get.edgelist(g2))
names(eg2)=c('source','target')
g3 = mst(g2)
plot(g3, layout= layout_on_grid, main = paste0('Rho = ',rho,', MST'))

