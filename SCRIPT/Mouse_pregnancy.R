 #Loading the Librarie
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF) 

seqdata <- read.delim("GSE60450_LactationGenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("SampleInfo_Corrected (1).txt", stringsAsFactors = TRUE)

#Check first 6 lines of the count data
head(seqdata)
#Check the sample information data
sampleinfo


#Set the rowname to geneid
rownames(seqdata) <- seqdata$EntrezGeneID
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)





# using substr, you extract the characters starting at position 1 and stopping at position 7 of
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
head(countdata)
#Check whether the names of the sampleinfo file and the colnames of the countdata are same or not
colnames(countdata)==sampleinfo$SampleName
#The DGEList() function is the core function in the edgeR package (part of Bioconductor) used to create an object specifically designed for Differential Gene Expression (DGE) analysis of RNA-Seq count data.
#It is the standard way to package your raw gene count matrix along with any associated sample information (like experimental group or library size) before proceeding with normalization and statistical testing.
#library size:It is the raw, unnormalized total count of all genes detected in that specific sample.
#It is calculated by simply summing the raw counts across all genes for a single sample (one column in your gene count matrix).
#The primary reason the DGEList() function is applied to count data and not sample information is because the purpose of the DGEList object is to store the raw, measured molecular data that will be used for statistical testing.

#Apply DGEList fucntion on Count Data
y <- DGEList(countdata)
# have a look at y
y

# paste the celltype and status together to make specific types with cell type info and status information
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")

# Convert to factor
group <- factor(group)

# Add the group information into the DGEList
y$samples$group <- group
y$samples
#Load the library for mouse
library(org.Mm.eg.db)
#Annotate the genes using the mouse as reference
#This step is crucial in RNA-Seq analysis because it takes the abstract IDs (which are typically Ensembl or Entrez IDs) used in the gene count matrix and retrieves human-readable, descriptive information for those genes.
ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
#Add the genes to your DGEList
y$genes <- ann
#$genes: This is a specific component within the DGEList object reserved for storing gene-specific annotation. It starts empty or contains minimal data when the DGEList is first created.
#Check the output
head(y$genes)
# Obtain CPMs
myCPM <- cpm(countdata)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))


# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)


#We Filter to keep the selected genes only
#When you remove genes, the total library size for each sample technically changes. Setting this to FALSE tells edgeR to re-calculate the library sizes based only on the counts of the genes that were kept. This is usually the desired behavior after filtering.
y <- y[keep, keep.lib.sizes=FALSE]

#Before Normalisation
# Get log2 counts per million
#log:This stabilizes the variance and makes the data distribution more symmetrical, which is ideal for plotting and clustering.
logcounts <- cpm(y,log=TRUE)
#The code boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2) is used for a critical Quality Control (QC) step in RNA-Seq analysis: visualizing the distribution of gene expression values across all your samples.
#This plot helps you quickly identify any samples that are outliers or have vastly different overall expression profiles before you proceed with the main analysis.
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
#After Nomalisation
# Apply normalisation to DGEList object
y <- calcNormFactors(y)

# Log2 CPM on normalised data
logcounts_norm <- cpm(y, log = TRUE)

# Create Boxplot after normalisation
boxplot(logcounts_norm, xlab = "", ylab = "Log2 counts per million", las = 2) 
abline(h=median(logcounts_norm),col="blue")
title("Boxplots of logCPMs (normalised)")


