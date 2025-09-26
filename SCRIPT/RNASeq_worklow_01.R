# install required bioconductor packages
BiocManager::install("EnhancedVolcano")


# install required packages
install.packages("ggthemes")

# load required packages
library(tidyverse)
library(airway)
library(DESeq2)

# get info about the data
?airway
#Read counts per gene for airway smooth muscle cell lines RNA-Seq experiment
#The dataset contains four cell lines in two conditions: control and treatment with dexamethasone.
# get the data
data(airway)

#counts_table
#The assay() function is a key accessor function in Bioconductor,
#used to extract the main experimental data matrix (or matrices) from a SummarizedExperiment or related object, such as a DESeqDataSet (which inherits from SummarizedExperiment).
counts_data <- assay(airway)

#The colData() function is a key Bioconductor accessor that extracts the sample metadata from a SummarizedExperiment or DESeqDataSet object. 
#It provides descriptive information about each column (sample) in the experiment.#coldata=column data

airway
head(assay(airway))
#rows=gene_names
#column_names=sample_names
# metadata (coldata)
col_data <- as.data.frame(colData(airway))

view(assay(airway))

view(col_data)

# making sure the row names in `col_data` (metadata) matches to columns names in `counts_data`
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) %in% rownames(col_data))
#%in%: This R operator checks for membership. The expression A %in% B returns a logical vector of the same length as A, indicating whether each element of A is found anywhere in B. 
#In this case, it checks if each sample ID from the count data is present in the metadata table.

## are they in same order?
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) == rownames(col_data))

# construct a DESeqDataSetFromMatrix data object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~dex # condition
)

#counts() → retrieves the gene × sample count matrix from a DESeq2 object (raw or normalized).
# pre-filtering: removing rows with low gene counts
# keep rows that have at least 10 reads total
rowSums(counts(dds))
rowSums(counts(dds)) >= 10

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]

# reference category ~ set the factor level
#relevel(...): Sets the reference factor level for the primary comparison variable (dex).
#By setting the ref = "untrt" (untreated), the subsequent DGE analysis will compare the treated ("trt") group relative to the untreated ("untrt") group.
#The log fold change (LFC) results will be log_2(expression in "trt"/expression in "untrt").

dds$dex <- relevel(dds$dex, ref = "untrt")


## perform differential gene expression analysis
dds <- DESeq(dds)

# save as results
res <- results(dds)
res
#BaseMean:The average of the normalized count values for that gene across all samples in the comparison.
#lfcSE:The standard error of the log 2 FoldChange estimate.Measures the precision of the log 2​FoldChange estimate.
#A smaller value indicates a more precise (less uncertain) estimate.
#Stat:The Wald test statistic used to test the hypothesis that the log 2​FoldChange is equal to zero. 
#Large absolute values of the statistic indicate stronger evidence against the null hypothesis (i.e., that the gene is not differentially expressed).
#padj:This is the main metric for significance.#we have to use this not the raw p value
#For genes to be reliably considered differentially expressed (DE), you should use the adjusted p-value (padj) and require it to be less than a predetermined significance threshold (α value), typically 0.05 or 0.01.
#Standard Cutoff: padj<0.05
#Stricter Cutoff: padj<0.01#Used when a lower false discovery rate is required, for example, if you plan to proceed immediately to costly wet-lab validation.

# exploring the results
summary(res)
#we have now 22389 genes out of 63,677 genes because the deseq2 has filtered out genes  that have no counts
#pdj < 0.1 is used bacause the dataset was noisy

## working with alpha (significance level)
# if alpha = 0.1 or 10% ~ 90% CI (default)
# if alpha = 0.01 or 1% ~ 99% CI
res_0.01 <- results(dds, alpha = 0.01)
summary(res_0.01)

# if alpha = 0.05 or 5% ~ 95% CI
res_0.05 <- results(dds, alpha = 0.05)
summary(res_0.05)



#The alpha parameter sets the significance threshold for the adjusted p-values (padj).
#By default, DESeq2 uses α=0.1
#Here, the we recalculates the results with stricter α values to see how the number of significant genes changes
# the more we decrease the confidence level,the better confidence we get about the results

# results name
resultsNames(dds)

# contrast
contrast_res <- results(dds, contrast = c("dex", "trt", "untrt"))
contrast_res
summary(contrast_res)
#The vector c("factor", "numerator", "denominator") tells DESeq2 to calculate the log_2 fold change for the level "trt" (numerator) relative to the level "untrt" (denominator) for the factor "dex".
#This is the same result as the default res because of the relevel step, but using contrast is a robust practice for specifying exact comparisons.











