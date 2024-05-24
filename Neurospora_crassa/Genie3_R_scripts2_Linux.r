conda activate m.peng_R
R
library("GENIE3");

setwd("/home/m.peng_cbs-niob.local/regulator/Genie3/Neurospora");
## exprMatr=read.table("microarray_data_paper_transporters_noMean_removeLowExpressionMax20.txt",sep="\t",head=T);
exprMatr=read.table("neurospora_transcriptome_data_noMutant.xls",sep="\t",head=T);


row.names(exprMatr)=exprMatr[,1]
exprMatr=exprMatr[,-1]

## exprMatr=exprMatr[exprMatr$max>20,]
## a=ncol(exprMatr)
## exprMatr=exprMatr[,-c(a,a-1)]


exprMatr=as.matrix(exprMatr)
set.seed(123) # For reproducibility of results
## weightMat <- GENIE3(exprMatr)

TFomes=read.table("JGI_neucr2_TFs.txt",head=T,sep="\t");
TFs=as.character(TFomes$gene)
regulators <- TFs


genes=rownames(exprMatr);
tfs=intersect(genes,regulators)
weightMat <- GENIE3(exprMatr, regulators=tfs, nCores=22, verbose=TRUE)



## linkList <- getLinkList(weightMat)
linkList <- getLinkList(weightMat, threshold=0.05)
write.table(linkList,"predicted_regulation_005.txt",sep="\t",col.names=NA)



linkList <- getLinkList(weightMat, threshold=0.025)
write.table(linkList,"predicted_regulation_0025.txt",sep="\t",col.names=NA)

linkList <- getLinkList(weightMat, threshold=0.01)
write.table(linkList,"predicted_regulation_001.txt",sep="\t",col.names=NA)

linkList <- getLinkList(weightMat, threshold=0.005)
write.table(linkList,"predicted_regulation_0005.txt",sep="\t",col.names=NA)


pbd=read.table("neurospora_PBC_genes.txt",head=T,sep="\t");

TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y=TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"PBD_related_regulation_prediction_0025.xls",sep="\t",col.names=NA)



linkList=read.table("predicted_regulation_0025.txt",head=T,sep="\t");
TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"Regulation_prediction_with_PBD_annotation_0025.xls",sep="\t",col.names=NA)


linkList=read.table("predicted_regulation_001.txt",head=T,sep="\t");
TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"Regulation_prediction_with_PBD_annotation_001.xls",sep="\t",col.names=NA)


linkList=read.table("predicted_regulation_0005.txt",head=T,sep="\t");
TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"Regulation_prediction_with_PBD_annotation_0005.xls",sep="\t",col.names=NA)





linkList <- getLinkList(weightMat, threshold=0.015)
write.table(linkList,"predicted_regulation_0015.txt",sep="\t",col.names=NA)

linkList=read.table("predicted_regulation_0015.txt",head=T,sep="\t");
TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"Regulation_prediction_with_PBD_annotation_0015.xls",sep="\t",col.names=NA)


linkList <- getLinkList(weightMat, threshold=0.02)
write.table(linkList,"predicted_regulation_002.txt",sep="\t",col.names=NA)

linkList=read.table("predicted_regulation_002.txt",head=T,sep="\t");
TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"Regulation_prediction_with_PBD_annotation_002.xls",sep="\t",col.names=NA)

