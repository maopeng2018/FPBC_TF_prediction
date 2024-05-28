
## Linux vrie1   /data2/regulator_prediction


library("GENIE3");

## setwd("E:\\Projects\\creative_ideas\\gene_prioritization\\new_TF\\RNA_seq\\raw_counts_select_samples\\no_mutant");
exprMatr=read.table("DEseq2_normalized_counts.xls",sep="\t",head=T);

row.names(exprMatr)=exprMatr[,1]
exprMatr=exprMatr[,-1]



## exprMatr=as.matrix(exprMatr)
exprMatr=as.matrix(exprMatr)
set.seed(123) # For reproducibility of results
## weightMat <- GENIE3(exprMatr)

TFomes=read.table("TFome_NRRL3_JGI_list.xls",head=T,sep="\t");
TFs=as.character(TFomes$ids)
regulators <- TFs

genes=rownames(exprMatr);
tfs=intersect(genes,regulators)
## weightMat <- GENIE3(exprMatr, regulators=tfs)   ## nCores=10
weightMat <- GENIE3(exprMatr,nCores=25,regulators=tfs)   ## nCores=10
## weightMat <- GENIE3(exprMatr,regulators=tfs)   ## nCores=10

## linkList <- getLinkList(weightMat)
linkList <- getLinkList(weightMat, threshold=0.0001)
write.table(linkList,"predicted_regulation_Genie3_treshold0.0001_FPKM10_noMutant2.txt",sep="\t",col.names=NA)


## linkList=read.table("predicted_regulation_Genie3_treshold0.0001_FPKM10_noMutant.txt",head=T,sep="\t");
pbd=read.table("PBC_crucial_genes_nrrl3.xls",head=T,sep="\t");

TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"PBD_related_regulation_prediction_Genie3_treshold0.0001_FPKM10_noMutant2.xls",sep="\t",col.names=NA)


linkList <- getLinkList(weightMat, threshold=0.005)
pbd=read.table("PBC_crucial_genes_nrrl3.xls",head=T,sep="\t");

TFs1=merge(pbd,linkList,by.x="ids",by.y="targetGene",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulatoryGene",all.y=TRUE)
write.table(TFs,"PBD_related_regulation_prediction_Genie3_treshold0.005_FPKM10_noMutant2.xls",sep="\t",col.names=NA)

