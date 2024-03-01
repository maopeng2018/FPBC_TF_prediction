setwd("F:\\Projects\\creative_ideas\\gene_prioritization\\new_TF\\Neurospora\\network_analysis\\new_analysis_noMutant");

exps=read.table("Neurospora_pvalues_3toolsCensus_E1fixLength_sig0.05.xls",head=T,sep="\t");

## select=exps
select=exps[exps$PBC_pvalues<0.05,]


print(colnames(select))

library("ComplexHeatmap")

mat = select[,c(24:ncol(select))]
mat[mat<0.00001]=0.00001
mat=-log10(mat)
rownames(mat) = select$deleted_gene

## select$deleted_gene=gsub("jgi\\|Aspni_NRRL3_1\\|","NRRL3_",select$deleted_gene)
ha = rowAnnotation(bar = select[,8],col = list(bar = c("sugar_transporter" = "blue", "sugar_catabolism" = "yellow", "PBD_CAZy" = "green")))

my_matrix=as.matrix(mat)

## library("colorRamp2") ; col_fun = colorRamp2(c(0,1,2,3,4,5), c("lightblue", "white","#FFCCCB", "red","darkred"));  colors=col_fun(seq(-3, 3))
library("RColorBrewer")
colors=rev(brewer.pal(n = 9, name = "RdBu"))[c(1,4:9)]

ht  <- Heatmap(my_matrix,column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE,rect_gp = gpar(col = "white", lwd = 1.5),col = colors)

ha_row <- rowAnnotation(names = row_anno_text(x = select$deleted_gene, just = "left", gp = gpar(fontsize = 8) ), show_annotation_name = FALSE)
ha_row2 <- rowAnnotation(names = row_anno_text(x = select$growthConditions, just = "left", gp = gpar(fontsize = 8) ), show_annotation_name = FALSE)
## ha_row3 <- rowAnnotation(names = row_anno_text(x = select$geneName, just = "left", gp = gpar(fontsize = 8) ), show_annotation_name = FALSE)

pdf("heatmap_predicted_TFs_3tools_E1_p005.pdf",  width=8, height=7,useDingbats=FALSE)

draw(ht + ha_row + ha_row2)
##draw(ht)
dev.off()



