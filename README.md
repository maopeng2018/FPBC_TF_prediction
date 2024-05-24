This project is about "Discovery of novel plant biomass conversion associated fungal transcription factors using a network-based approach"

Following steps are required to perfrom the prediction of transcription factors (TFs) related to fungal plant biomass conversion (FPBC):

Step 1: Collection a large set of transcriptome data. In this case, the related transcriptome data for Aspergillus niger and Neurospora crassa have been prepared and normalized.


Step 2: Inferencing gene regulatory networks (GRNs) from transcriptome data:
    In this project, we tried three different methods: 
    Method 1, Aracne   see https://pubmed.ncbi.nlm.nih.gov/27153652/, you can find related command to run this analysis is attached in corresponding folders.           

    Method 2, Genie3, see https://pubmed.ncbi.nlm.nih.gov/20927193/ . 
              It's a well-documented R packages, https://bioconductor.org/packages/release/bioc/html/GENIE3.html. You can also find our R scripts in corresponding folders. 

    Method 3, KBoost, see https://pubmed.ncbi.nlm.nih.gov/34326402/.
              It's also a well-documented R packages, https://www.bioconductor.org/packages/release/bioc/html/KBoost.html. You can also find our R scripts in corresponding folders. 

Step 3, Combine the consensus network prediction results of above three algorithms using R script "mergemerge_table - Aracne_Genie3_Kboost_predictions.r".   

Step 4, Calculate the enrichment function for regulons of each TF using the R script "hypergeometric test_3ToolsConsensus_v5_pbcClass_E1fixlength.R".

Step 5, You can visualizing your results using heatmap by R script, "heatmap_prediction_results_v2.r".
