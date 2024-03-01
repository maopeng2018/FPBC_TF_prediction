setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_TF\\RNA_seq\\raw_counts_select_samples\\no_mutant");



net=read.table("Regulation_prediction_with_PBD_annotation_AracneE1.xls",head=T,sep="\t");
net=net[net$pvalue<1e-1,]
net <- net[order(net$pvalue), ]
net2=net[,c(2,4,6)]
colnames(net2)=c("regulator","target","weight")
aracne=net2
len=nrow(aracne)


net=read.table("Regulation_prediction_with_PBD_annotation_Genie3_w0.005.xls",head=T,sep="\t");
net <- net[order(-net$weight), ]
##net=net[1:len,]
net2=net[,c(2,4,7)]
colnames(net2)=c("regulator","target","weight")
len2=nrow(net2)
genie3=net2[1:len,]



net=read.table("Regulation_kboostProb00002_with_PBD_annotation.xls",head=T,sep="\t");
net <- net[order(-net$score), ]
##net=net[1:len,]
net2=net[,c(2,4,6)]
colnames(net2)=c("regulator","target","weight")
len3=nrow(net2)
kboost=net2[1:len,]


grn_list<-list(aracne,genie3,kboost)
combine <- Reduce(function(x, y) merge(x, y, by = c("regulator","target"),all.x=TRUE,all.y=TRUE), grn_list)
naCount<-function(x){length(x)-sum(is.na(x))}
combine$topRanked <- apply(combine[, 3:ncol(combine)], 1, naCount)

pbd=read.table("PBC_crucial_genes_nrrl3.xls",head=T,sep="\t");
TFs1=merge(pbd,combine,by.x="ids",by.y="target",all.y = TRUE)
TFs=merge(pbd,TFs1,by.x="ids",by.y="regulator",all.y=TRUE)
TFs2=TFs[TFs$topRanked>1,]
write.table(TFs2,"Regulation_prediction_with_PBD_annotation_3ToolsConsensus2_E1fixLength.xls",sep="\t",col.names=NA)

TFs3=TFs[TFs$topRanked>0,]
write.table(TFs3,"Regulation_prediction_with_PBD_annotation_3ToolsCombine_E1fixLength.xls",sep="\t",col.names=NA)

