# Install packages

BiocManager::install("recount")
BiocManager::install("apeglm")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("clusterProfiler")
install.packages("msigdbr")
#DESeq2 and apeglm: Core packages for DGE analysis and LFC shrinkage.
#recount: Used to easily download public RNA-seq data from the recount2 database.
#EnhancedVolcano: For advanced visualization of results (Volcano plot).
#clusterProfiler and msigdbr: Used for downstream pathway and gene set enrichment analysis (though not used in the executed script).
#tidyverse: A collection of packages for data manipulation and visualization (e.g., dplyr,ggplot2).
#EnsDb.Hsapiens.v86: A database package for gene annotation.

# Load libraries
library(tidyverse)
library(recount)
library(DESeq2)

# Problem: Effects of EWS-FLI1 KD in Ewing sarcoma

# RQ: What is the effect of knocking down EWS-FLI1 in Ewing sarcoma? How does it change cellular behavior?



#Ewing sarcoma (ES) is a highly aggressive malignant bone and soft tissue tumor, primarily affecting children and young adults.
#The tumor is defined and driven by a characteristic chromosomal rearrangement that creates a potent fusion oncogene.
#EWS-FLI1 is the primary oncogenic driver for approximately 85% of Ewing sarcoma (ES) cases.
#It is an abnormal, chimeric transcription factor resulting from a chromosomal translocation, typically t(11;22).
#EWS-FLI1 is the primary oncogenic driver for approximately 85% of Ewing sarcoma (ES) cases. It is an abnormal, chimeric transcription factor resulting from a chromosomal translocation, typically t(11;22).
#The observed effects of EWS-FLI1 knockdown are essentially the reversal of the oncogenic phenotype that the protein causes.
#The knockdown (KD) of EWS-FLI1 is a critical experimental manipulation used to confirm its role and identify its downstream targets.
#The results of the differential gene expression (DGE) analysis we are running will identify the full list of genes whose expression is controlled (activated or repressed) by this fusion protein.




# The dataset we use here is downloaded from recount2
# Find a project of interest (note: it's hit or miss, the database is not limitless)
project_info <- abstract_search(query = "Ewing sarcoma")


# select study
download_study("SRP015989")

# load the data
load("SRP015989/rse_gene.Rdata")

# examine the read counts
counts_data <- assay(rse_gene)

# examine the coldata
col_data <- as.data.frame(colData(rse_gene))

# fix the colData to give a column with the appropriate groups
rse_gene$condition <- c(rep("shCTR", 3), rep("ShEF1", 4))

# check colnames of count_data and rownames of col_data
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) %in% rownames(col_data))

# are they in same order
all(colnames(counts_data) == rownames(col_data))

# perform DESeq2 analysis using DESeqDataSetFromMatrix (not suitable for this condition)
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~condition
)

# perform DESeq2 analysis using DESeqDataSet
dds <- DESeqDataSet(rse_gene, design = ~condition)

# normalize the data
rld <- rlog(dds)
plotPCA(rld)

# perform DEG analysis
dds <- DESeq(dds)

# save as results
res <- results(dds)
res

# exploring the results
summary(res)

# working with alpha (significance level)
# if alpha = 0.1 or 10% ~ 90% CI (default)
# if alpha = 0.01 or 1% ~ 99% CI
res_0.01 <- results(dds, alpha = 0.01)
summary(res_0.01)

# if alpha = 0.05 or 5% ~ 95% CI
res_0.05 <- results(dds, alpha = 0.05)
summary(res_0.05)

# results name
resultsNames(dds)

# contrast
contrast_res <- results(dds, contrast = c("condition", "shEF1", "shCTR"))
contrast_res
summary(contrast_res)

# Visualization (What? Why? How?)
# MA Plot (M vs A Plot)
# What: An MA plot represents the log-fold change (M) on the y-axis and the average expression (A) on the x-axis for each gene or feature.
# Why: It is used to visualize differential expression between two conditions.The log-fold change (M) gives an idea of the magnitude of change, and the average expression (A) helps identify if the change is dependent on the expression level.
plotMA(dds)

# Interpretations
# 1. Points on the plot represent genes.
# 2. Genes with significant differential expression are often `found at the extremes` of the plot.
# 3. Upregulated genes are at the top
# 4. Downregulated genes at the bottom
# 5. Non-differentially expressed genes are centered around zero on the y-axis.




