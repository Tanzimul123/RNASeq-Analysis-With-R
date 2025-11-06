 
#Loading the Libraries
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(GEOquery)


#Read the data into R
Yeast_data <- read.delim("FPKM_NR_CR(YEAST under restricted calorie).txt", stringsAsFactors = FALSE)
View(Yeast_data)


#get metadata
library(GEOquery)
gse_data <- getGEO("GSE53720", GSEMatrix = TRUE)
# The metadata (sample information) is typically stored in the 'pData' slot
 metadata3 <- pData(phenoData(gse_data[[1]]))
View(metadata3)

#subsetting
metadata3=metadata3[,c(2,18,43)]
View(metadata3)


#renamingColumns
colnames(metadata3)<-c('ID','Sample','status')
View(metadata3)

# using substr, you extract the characters starting at position 1 and stopping at position 7 of
colnames(Yeast_data) <- substr(colnames(Yeast_data),start=1,stop=8)
head(Yeast_data)
#Check whether the names of the sampleinfo file and the colnames of the countdata are same or not
colnames(Yeast_data)[4:7]==metadata3$Sample


#Apply DGEList fucntion on Count Data
y <- DGEList(Yeast_data)
# have a look at y
y

# paste the sample and status together 
group <- paste( metadata3$Sample,metadata3$status,sep=".")

# Convert to factor
group <- factor(group)

# Add the group information into the DGEList
y$samples$group <- group
y$samples


# Convert the FPKM matrix (y$counts) to log2(FPKM + pseudocount)
logcounts <- log2(y$counts + 1)

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 fpkm",las=2)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of fpkms (normalised)")


# 1. Convert the 'status' column to a factor
metadata3$status <- as.factor(metadata3$status)

# 2. Check the levels to ensure they are correct
levels(metadata3$status)



# Re-create the color vector based on the fixed factor
col.status <- c("purple", "orange")[metadata3$status]





# Plot MDS plot using the DGEList object and the color vector
plotMDS(y, col = col.status)

# Add the legend using the same colors and the factor levels
legend("topleft", 
       fill = c("purple", "orange"), 
       legend = levels(metadata3$status))

# Add a title
title("MDS Plot by Calorie Restriction Status (GSE53720)")

## Add labels to the points (using the sample column names, if available)
plotMDS(y, col = col.status, labels = colnames(y))

#Cluster 1 (Purple - Calorie-restricted): The samples s_3 CR A and s_4 CR are separated from the other group, primarily along the horizontal axis (Leading logFC dim 1), which accounts for 75% of the variance. This means Calorie-restriction is the major biological difference driving the gene expression patterns.
#Cluster 2 (Orange - Non-restricted): The samples s_1 NR B and s_2 NR B (visible as one single orange label, which likely represents two points overlapping due to high similarity) are clustered tightly together.

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)


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
png(file="new_High_var_genes_YEAST.heatmap.png")

# 3. Create the aheatmap plot (Fix: Missing closing parenthesis)
aheatmap(
  highly_variable_lcpm,
  col = rev(morecols(50)),
  main = "Top 500 most variable genes across samples")  
# 4. Close the graphics device to save the file
dev.off()

getwd()



#CREATING DESIGN MATRIX
group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design


 

# 1. Correctly define the 'group' factor for your 4 samples
#    Adjust the labels based on your actual sample order
group <- factor(c("NR", "NR", "CR", "CR"))

# 2. Recreate the design matrix correctly
design <- model.matrix(~0 + group) 

# 3. Assign clear column names
colnames(design) <- levels(group) 

# 4. Re-run voom
v <- voom(y, design, plot = TRUE)


#FITTING LINEAR MODEL
fit <- lmFit(v)
names(fit)


#CREATING CONTRAST AND FITTING
# Compare the Calorie-Restricted group (CR) against the Non-Restricted group (NR)
# The coefficient for the second group (NR) is subtracted from the first (CR).
cont.matrix <- makeContrasts(CR_vs_NR = CR - NR, levels = design) 

# Re-run the fitting steps with the corrected contrast matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summary(decideTests(fit.cont))

