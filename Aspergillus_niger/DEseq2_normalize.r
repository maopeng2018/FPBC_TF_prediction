setwd("E:\\Projects\\creative_ideas\\gene_prioritization\\new_TF\\RNA_seq\\raw_counts_select_samples\\no_mutant\\remove_mutant")

library(corrplot)

fpkms <- read.table("combine_counts_published_papers_selected.xls", header=T, sep="\t")
fpkm3<-fpkms[, -1]
rownames(fpkm3)=fpkms[,1]
fpkm3=fpkm3[rowSums(fpkm3)>0,];
fpkm3=as.matrix(fpkm3)
M <- cor(fpkm3)
## corrplot(M, method="circle")
corrplot.mixed(M)
## pdf("FPKM_correlation.pdf")
## corrplot.mixed(M,tl.cex=0.3,cl.cex=0.3,number.cex=0.3,pch.cex=0.3,pointsize=15,，bg="white",align=c("c"))
## dev.off();
write.table(M,file="samples_correlation.matrix.xls",sep="\t",col.names=NA)

library(matrixStats)
fpkm2=fpkm3[rowMaxs(fpkm3)>10,];


library("DESeq2");

lineas=fpkm2


countMatrix <- as.matrix(lineas)
countMatrix=round(countMatrix)
### countMatrix=apply(countMatrix, 1:ncol(countMatrix), round)
head(countMatrix)

## Create colData to run the comparisons
comp<- matrix(ncol=2, nrow=0)
for (n in colnames(countMatrix)){
    new<-unlist(strsplit(n, "\\.1$|\\.2$|\\.3$"))[1]  ##### to format the name of duplicate experiment
  line<-cbind(n, new)
  #print (line)
  comp<-rbind(comp, line)
  #print (lineas)
}
expt_design <- data.frame(comp)
print(expt_design)


dds <- DESeqDataSetFromMatrix(countMatrix, colData=expt_design, design= ~ new)
#dds$condition <- relevel(dds$condition, levels=expt_design$new)
##dds2 <- DESeq(dds) 
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalizedCounts=counts(dds, normalized=TRUE)


write.table(normalizedCounts,file="DEseq2_normalized_counts2.xls",sep="\t",col.names=NA)


