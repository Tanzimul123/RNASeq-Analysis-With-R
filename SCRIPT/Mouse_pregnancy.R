 #Loading the Librarie
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF) 

seqdata <- read.delim("GSE60450_LactationGenewiseCounts.txt", stringsAsFactors = FALSE)
View(seqdata)


# preparing metadata THROUGH CODE
meta <- read.delim("GSE60450_series_matrix.txt", 
                   header = FALSE, 
                   comment.char = "!")
head(meta)
# Read all lines
lines <- readLines( 'GSE60450_series_matrix.txt')

# Extract metadata rows
titles <- grep("!Sample_title", lines, value = TRUE)
sources <- grep("!Sample_source_name_ch1", lines, value = TRUE)
characteristics <- grep("!Sample_characteristics_ch1", lines, value = TRUE)

# Convert metadata fields to table
get_values <- function(x) {
  strsplit(sub("!Sample_.*?\t", "", x), "\t")[[1]]
}

sample_titles <- get_values(titles)
sample_sources <- get_values(sources)
sample_characteristics <- get_values(characteristics)

metadata <- data.frame(
  title = sample_titles,
  source = sample_sources,
  condition = sample_characteristics
)

View(metadata)



# Read the sample information into R PROVIDED IN THE BLOG
sampleinfo <- read.delim("SampleInfo_Corrected (GSE60450).txt", stringsAsFactors = TRUE)

#metadata prepared through excel after downloading from GREIN
sampleinfo2 <- read.delim("GSE60450_full_metadata.txt", stringsAsFactors = TRUE)


#Check first 6 lines of the count data
head(seqdata)
#Check the sample information data
sampleinfo
sampleinfo2

#Set the rowname to geneid
rownames(seqdata) <- seqdata$EntrezGeneID
# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)

#



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


## Let's choose purple for basal and orange for luminal
col.cell <- c("purple", "orange")[sampleinfo$CellType]
#Plot MDS plot
plotMDS(y, col = col.cell)
legend("topleft", fill = c("purple", "orange"), legend = levels(sampleinfo$CellType))
title("MDS Plot by Cell Type")


# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
var_genes
View(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)


 
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)


#Load the Libraries
library(RColorBrewer)
library(gplots)
library(NMF)


# 1. Define the palette and color function
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# 2. Open the PNG graphics device
png(file="new_High_var_genes.heatmap.png")

# 3. Create the aheatmap plot (Fix: Missing closing parenthesis)
aheatmap(
  highly_variable_lcpm,
  col = rev(morecols(50)),
  main = "Top 500 most variable genes across samples")  
# 4. Close the graphics device to save the file
dev.off()

getwd()
# preparing design matrix
group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

# VOOM TRANSFORMATION
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)



# Keep running this until the console says "null device" or "Error in dev.off(): cannot shut down device 1 (the null device)"
# Then re-run the voom command:
v <- voom(y, design, plot = TRUE)

#CREATING LINEAR MODEL
fit <- lmFit(v)
names(fit)

#CREATING CONTRAST AND FITTING
# 1. Create the contrast matrix (FIX: Missing closing parenthesis for makeContrasts)
cont.matrix <- makeContrasts(B.PregVsLac = basal.pregnant - basal.lactate, 
                             levels = design)

# 2. Fit the contrasts to the linear model fit
fit.cont <- contrasts.fit(fit, cont.matrix)

# 3. Apply the empirical Bayes smoothing to the contrasts fit
fit.cont <- eBayes(fit.cont)

# 4. Summarize the results of the differential expression tests
summary(decideTests(fit.cont))
# Create the summa.fit object needed for the plotting status argument
summa.fit <- decideTests(fit.cont)

#SAVING RESULT AND VISUALIZATION
# Extract all genes from the contrast fit object
all_genes <- topTable(fit.cont, coef = 1, number = Inf, adjust = "fdr")
# write all genes results in a csv file
write.csv(all_genes,"all_genes.csv")

#CREATING MD PLOT AND VOLCANO PLOT
# We want to highlight the significant genes.
par(mfrow=c(1,2))

# 1. Plot the Mean-Difference (MD) Plot (FIX: Missing closing parenthesis and hl.col value)
plotMD(fit.cont,
       coef=1,
       status=summa.fit[,"B.PregVsLac"], 
       values = c(-1, 1), 
       hl.col=c("blue", "red") # <- FIXED: Added the second color for positive logFC, and closed the parenthesis
)

# 2. Plot the Volcano Plot (FIX: Missing closing parenthesis for volcanoplot)
volcanoplot(fit.cont,
            coef=1,
            highlight=100,
            names=fit.cont$genes$SYMBOL, 
            main="B.PregVsLac") # <- FIXED: Closed the parenthesis


# gene ontology analysis
BiocManager::install("GO.db")
go <- goana(fit.cont, coef = "B.PregVsLac", species = "Mm")
topGO(go, n = 10)


