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

