setwd("F:\\Projects\\creative_ideas\\gene_prioritization\\new_TF\\RNA_seq\\raw_counts_select_samples\\no_mutant");
data=read.table("pvalues_TFs_PBCclassEnrichment_3toolsConsensus_E1fixLength_all.xls",head=T,sep="\t");
tfs=data[,c("JGI_geneID","PBC_pvalues","known")]
tfs2=tfs[order(-tfs$PBC_pvalues),]

tfs2$orders=1:nrow(tfs2)

tfs2$PBC_pvalues=-log10(tfs2$PBC_pvalues)

tfs2$PBC_pvalues[tfs2$PBC_pvalues>10]=10

library(ggplot2)

ggplot(tfs2, aes(x=orders, y=PBC_pvalues)) + 
   geom_point(color='grey',shape=21,size=3,aes(fill=factor(known)))+
   scale_fill_manual(values=c('cyan', 'red'))
 scale_y_continuous(limits = c(0, 10),breaks = c(0,2, 4, 6,8,10))+
 scale_x_continuous(limits = c(0, 900),breaks = c(0,200, 400, 600,800))+ 
 geom_hline(yintercept=2,color = "darkgrey")
 
 
ggplot(tfs2, aes(x=orders, y=PBC_pvalues)) + 
   geom_point(color='cyan',shape=21,size=3,aes(fill=factor(known)),alpha=0.4)+
   scale_fill_manual(values=c('cyan', 'darkred'))+
 scale_y_continuous(limits = c(0, 10),breaks = c(0,2, 4, 6,8,10))+
 scale_x_continuous(limits = c(0, 900),breaks = c(0,200, 400, 600,800))+ 
 geom_hline(yintercept=2,color = "darkgrey")
 	
	
	
	
	