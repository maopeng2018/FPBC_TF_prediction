setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_TF\\RNA_seq\\raw_counts_select_samples\\no_mutant");


DrawnBall_ID=2
wtDrawnID=1
wtTotal=440 ### CBS strain is 464, NRRL3 is 451 (only 420 with FPKM>10),  440 PBC genes included in the GRN network
TotalBall=11653    ### CBS strain is 13998, NRRL3 is 11846,  11654 genes are used for building GRN network， 11466 genes included in the GRN network

pValues_cal<- function (x) { 
blackBall=TotalBall-as.numeric(wtTotal);
##pval<-return(tryCatch(phyper(x[wtDrawnID],x[wtTotalID],(TotalBall-x[wtTotalID]),DrawnBall),error=function(e)NA))
## return(1-phyper(x[wtDrawnID],wtTotal,blackBall,x[DrawnBall_ID]))  ### this is probability the expected white balled larger than the observed onee 
return(1-phyper(x[wtDrawnID]-1,wtTotal,blackBall,x[DrawnBall_ID]))  ### this is probability the expected white balled larger or equal than the observed onee 
}




## tfs=read.table("TFome_NRRL3_JGI_list.xls",sep="\t",head=T)
reg=read.table("Regulation_prediction_with_PBD_annotation_3ToolsConsensus2_E1fixLength.xls",sep="\t",head=T)
pbcs=read.table("PBC_annotation_aniger.txt",sep="\t",head=T)

df=pbcs
df=df[df$group %in% c("PBD_cazy"),]  ###  %in% will help us to select rows based on a list of matches
df1=df[grep("Starch|starch",df$Substrate),]; df1$Substrate=rep("Starch",nrow(df1));
df2=df[grep("Cellulose|cellulose",df$Substrate),]; df2$Substrate=rep("Cellulose",nrow(df2));
df3=df[grep("Mannan|mannan",df$Substrate),]; df3$Substrate=rep("Mannan",nrow(df3));
df4=df[grep("Xyloglucan|xyloglucan",df$Substrate),]; df4$Substrate=rep("Xyloglucan",nrow(df4));
df5=df[grep("Xylan|xylan",df$Substrate),]; df5$Substrate=rep("Xylan",nrow(df5));
df6=df[grep("Pectin|pectin",df$Substrate),]; df6$Substrate=rep("Pectin",nrow(df6));

##gsets<-c("CAZy: Starch","CAZy: Cellulose","CAZy: Mannan","CAZy: Xyloglucan","CAZy: Xylan","CAZy: Pectin", "Metabolism: galacturonic acid","Metabolism: Rhamnose", "Metabolism: PCP","Metabolism: Galactose","Metabolism: PPP","Metabolism: Glycolysis","Metabolism: Glycerol",  "Metabolism: TCA")

df=pbcs
df=df[df$group %in% c("sugar_metabolism"),]  ###  %in% will help us to select rows based on a list of matches
df7=df[grep("galacturonic",df$Enzyme_family.pathway),];           df7$Enzyme_family.pathway=rep("Galacturonic acid",nrow(df7));
df8=df[grep("rhamnose|Rhamnose",df$Enzyme_family.pathway),];      df8$Enzyme_family.pathway=rep("Rhamnose",nrow(df8));
df9=df[grep("PCP",df$Enzyme_family.pathway),];                    df9$Enzyme_family.pathway=rep("PCP",nrow(df9));
df10=df[grep("oxidoreductive|Oxidoreductive|Leloir|leloir|Doudoroff|doudoroff",df$Enzyme_family.pathway),]; df10$Enzyme_family.pathway=rep("Galactose",nrow(df10));
df11=df[grep("Glycerol|glycerol",df$Enzyme_family.pathway),];      df11$Enzyme_family.pathway=rep("Glycerol",nrow(df11));
df12=df[grep("PPP",df$Enzyme_family.pathway),];                    df12$Enzyme_family.pathway=rep("PPP",nrow(df12));
df13=df[grep("Glycolysis|glycolysis",df$Enzyme_family.pathway),];  df13$Enzyme_family.pathway=rep("Glycolysis",nrow(df13));
df14=df[grep("TCA|glyoxylate",df$Enzyme_family.pathway),];         df14$Enzyme_family.pathway=rep("TCA,glyoxylate cycle",nrow(df14));

df=pbcs
df15=df[df$group %in% c("sugar_transporter"),]
df16=df[df$group %in% c("transcription_factor"),]
df17=df[df$group %in% c("PBD_cazy"),]
df18=df[df$group %in% c("sugar_metabolism"),]



gs<-list(df$GeneID,df15$GeneID,df16$GeneID,df17$GeneID,df18$GeneID,df1$GeneID,df2$GeneID,df3$GeneID,df4$GeneID,df5$GeneID,df6$GeneID,df7$GeneID,df8$GeneID,df9$GeneID,df10$GeneID,df11$GeneID,df12$GeneID,df13$GeneID,df14$GeneID)
gsets<-c("PBC","transporter","tfs","cazy","metabolism","CAZy_starch","CAZy_cellulose","CAZy_mannan","CAZy_xyloglucan","CAZy_xylan","CAZy_pectin", "Metabolism_galacturonicAcid","Metabolism_rhamnose", "Metabolism_PCP","Metabolism_galactose","Metabolism_glycerol","Metabolism_PPP", "Metabolism_glycolysis", "Metabolism_TCA")
names(gs) <- gsets
pbc=length(gs$PBC);
cazy=length(gs$cazy)
metabolism=length(gs$metabolism)
transporter=length(gs$transporter);
regulator=length(gs$tfs);


starch=length(gs$CAZy_starch);
cellulose=length(gs$CAZy_cellulose);
mannan=length(gs$CAZy_mannan);
xyloglucan=length(gs$CAZy_xyloglucan);
xylan=length(gs$CAZy_xylan);
pectin=length(gs$CAZy_pectin);
galacturonicAcid=length(gs$Metabolism_galacturonicAcid);
rhamnose=length(gs$Metabolism_rhamnose);
PCP=length(gs$Metabolism_PCP);
galactose=length(gs$Metabolism_galactose);
glycerol=length(gs$Metabolism_glycerol);
PPP=length(gs$Metabolism_PPP);
glycolysis=length(gs$Metabolism_glycolysis);
TCA=length(gs$Metabolism_TCA);


tfs=unique(as.character(reg$ids))
## gaaR="jgi|Aspni_NRRL3_1|8195";
## araR="jgi|Aspni_NRRL3_1|7564";
## hapX="jgi|Aspni_NRRL3_1|511";
## creA="jgi|Aspni_NRRL3_1|5946";
StatisValues <- data.frame(matrix(NA, nrow = length(tfs), ncol = 0));
genes=c();


for(i in 1:length(tfs)){
tfr=reg[reg$ids %in% tfs[i],]
totalTarget=nrow(tfr)
pbc_t=length(intersect(as.character(tfr$ids.y),gs$PBC))
cazy_t=length(intersect(as.character(tfr$ids.y),gs$cazy))
metabolism_t=length(intersect(as.character(tfr$ids.y),gs$metabolism))
transporter_t=length(intersect(as.character(tfr$ids.y),gs$transporter))
regulator_t=length(intersect(as.character(tfr$ids.y),gs$tfs))
starch_t=length(intersect(as.character(tfr$ids.y),gs$CAZy_starch))
cellulose_t=length(intersect(as.character(tfr$ids.y),gs$CAZy_cellulose))
mannan_t=length(intersect(as.character(tfr$ids.y),gs$CAZy_mannan))
xyloglucan_t=length(intersect(as.character(tfr$ids.y),gs$CAZy_xyloglucan))
xylan_t=length(intersect(as.character(tfr$ids.y),gs$CAZy_xylan))
pectin_t=length(intersect(as.character(tfr$ids.y),gs$CAZy_pectin))
galacturonicAcid_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_galacturonicAcid))
rhamnose_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_rhamnose))
PCP_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_PCP))
galactose_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_galactose))
glycerol_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_glycerol))
PPP_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_PPP))
glycolysis_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_glycolysis))
TCA_t=length(intersect(as.character(tfr$ids.y),gs$Metabolism_TCA))


wtTotal=pbc;          p1=pValues_cal(c(pbc_t,totalTarget));
wtTotal=cazy;         p2=pValues_cal(c(cazy_t,totalTarget));
wtTotal=metabolism;   p3=pValues_cal(c(metabolism_t,totalTarget));
wtTotal=transporter;  p4=pValues_cal(c(transporter_t,totalTarget));
wtTotal=regulator;    p5=pValues_cal(c(regulator_t,totalTarget));
wtTotal=starch;          p6=pValues_cal(c(starch_t,totalTarget));
wtTotal=cellulose;          p7=pValues_cal(c(cellulose_t,totalTarget));
wtTotal=mannan;          p8=pValues_cal(c(mannan_t,totalTarget));
wtTotal=xyloglucan;          p9=pValues_cal(c(xyloglucan_t,totalTarget));
wtTotal=xylan;          p10=pValues_cal(c(xylan_t,totalTarget));
wtTotal=pectin;          p11=pValues_cal(c(pectin_t,totalTarget));
wtTotal=galacturonicAcid;          p12=pValues_cal(c(galacturonicAcid_t,totalTarget));
wtTotal=rhamnose;          p13=pValues_cal(c(rhamnose_t,totalTarget));
wtTotal=PCP;          p14=pValues_cal(c(PCP_t,totalTarget));
wtTotal=galactose;          p15=pValues_cal(c(galactose_t,totalTarget));
wtTotal=glycerol;          p16=pValues_cal(c(glycerol_t,totalTarget));
wtTotal=PPP;          p17=pValues_cal(c(PPP_t,totalTarget));
wtTotal=glycolysis;          p18=pValues_cal(c(glycolysis_t,totalTarget));
wtTotal=TCA;          p19=pValues_cal(c(TCA_t,totalTarget));


sig=c(totalTarget,pbc_t,cazy_t,metabolism_t,transporter_t,regulator_t,starch_t,cellulose_t,mannan_t,xyloglucan_t,xylan_t,pectin_t,galacturonicAcid_t,
    rhamnose_t,PCP_t,galactose_t,glycerol_t,PPP_t,glycolysis_t,TCA_t,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19);
StatisValues=rbind(StatisValues,sig);
genes=c(genes,tfs[i])
}

StatisValues=cbind(genes,StatisValues)



## padj=p.adjust(pvalues,method="BH")
## result=data.frame(compare2,pvalues,padj)


knownTF=read.table("knownTF.xls",head=T,sep="\t");
TFs2=merge(knownTF,StatisValues,by.x="genes",by.y="genes",all.y=TRUE)
colnames(TFs2)=c("JGI_geneID","TFs_function","Regulon_genes","PBC_genes","Cazy_genes","metabolism_genes","transporter_genes","regulator_genes","starch_genes","cellulose_genes","mannan_genes","xyloglucan_genes","xylan_genes","pectin_genes","galacturonicAcid_genes","rhamnose_genes","PCP_genes","galactose_genes","glycerol_genes","PPP_genes","glycolysis_genes","TCA_genes",
"PBC_pvalues","Cazy_pvalues","metabolism_pvalues","transporter_pvalues","regulator_pvalues","starch_pvalues","cellulose_pvalues","mannan_pvalues","xyloglucan_pvalues","xylan_pvalues","pectin_pvalues","galacturonicAcid_pvalues","rhamnose_pvalues","PCP_pvalues","galactose_pvalues","glycerol_pvalues","PPP_pvalues","glycolysis_pvalues","TCA_pvalues"
)

sigs=c();
for (k in 23:ncol(TFs2)){
TFs2[,k]=p.adjust(TFs2[,k],method="BH")
sig=TFs2[TFs2[,k]<0.05,1]
sigs=union(sigs,as.character(sig))
}

TFs3=TFs2[TFs2$JGI_geneID %in% sigs,]

blast=read.table("best_hits_between_Aspergillus_CBS_NRRL3.txt",head=T,sep="\t");
others=read.table("known_other_TF.txt",head=T,sep="\t");
cbs=merge(blast,others,by.x="query",by.y="anNumber",all.x=TRUE)
cbs2=cbs[,c("ids","query","geneName","Annotation","Reference")]

TFs4=merge(TFs3,cbs2,by.x="JGI_geneID",by.y="ids",all.x=TRUE)
write.table(TFs4,file="pvalues_TFs_PBCclassEnrichment_3toolsConsensus_E1fixLength_significanceP005_v2.xls",sep="\t",col.names=NA);

TFs22=merge(TFs2,cbs2,by.x="JGI_geneID",by.y="ids",all.x=TRUE)
write.table(TFs22,file="pvalues_TFs_PBCclassEnrichment_3toolsConsensus_E1fixLength_all_v2.xls",sep="\t",col.names=NA);
