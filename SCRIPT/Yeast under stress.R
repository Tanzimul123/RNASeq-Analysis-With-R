 
#Loading the Libraries
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(GEOquery)
library(org.Sc.sgd.db)
library(AnnotationDbi)


#Read the data into R
Yeast_data <- read.delim("FPKM_NR_CR(YEAST under restricted calorie).txt", stringsAsFactors = FALSE)
View(Yeast_data)
# Find rows where the gene_id is NOT a valid number (i.e., missing/NA/blank)
Yeast_data_clean <- Yeast_data[!is.na(Yeast_data$gene_id) & Yeast_data$gene_id != "", ]
# Find and remove any explicit duplicate Entrez IDs
Yeast_data_clean <- Yeast_data_clean[!duplicated(Yeast_data_clean$gene_id), ]

# 3. Set the row names using the now-unique 'gene_id' column.
rownames(Yeast_data_clean) <- Yeast_data_clean$gene_id

View(Yeast_data_clean)

 

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
colnames(Yeast_data_clean) <- substr(colnames(Yeast_data_clean),start=1,stop=8)
head(Yeast_data_clean)
#Check whether the names of the sampleinfo file and the colnames of the countdata are same or not
colnames(Yeast_data_clean)[4:7]==metadata3$Sample


#Apply DGEList fucntion on Count Data
y <- DGEList(Yeast_data_clean)
# have a look at y
y

# paste the sample and status together 
group <- paste( metadata3$Sample,metadata3$status,sep=".")

# Convert to factor
group <- factor(group)

# Add the group information into the DGEList
y$samples$group <- group
y$samples

#annotation
library(org.Sc.sgd.db) 
library(AnnotationDbi)

keytypes(org.Sc.sgd.db)
#Annotate the genes using the yeast as reference
#ann <- select(org.Sc.sgd.db,keys=rownames(y$genes$ORF),columns=c("ENTREZID","DESCRIPTION","GENENAME","GO", "ORF")keytype = "ORF")

# Annotate genes using org.Sc.sgd.db
ann <- AnnotationDbi::select(org.Sc.sgd.db,
                             keys = Yeast_data_clean$gene_id,
                             columns = c("ENTREZID","DESCRIPTION","GENENAME","GO", "ORF"),
                             keytype = "ORF")

#Add the genes to your DGEList
y$genes <- ann

#Check the output
head(y$genes)


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
 
select_var



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
cont.matrix
# Re-run the fitting steps with the corrected contrast matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summary(decideTests(fit.cont))

# Extract all genes from the contrast fit object
all_genes <- topTable(fit.cont, coef = 1, number = Inf, adjust = "fdr")
# write all genes results in a csv file
write.csv(all_genes,"all_genes.csv")

# Assuming the contrast name you used was 'CR_vs_NR'
# and the results table from the summary(decideTests) step is named 'summa.fit'
# 1. Store the results of decideTests() into the variable 'summa.fit'
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
par(mfrow=c(1,2))

plotMD(fit.cont, 
       coef = 1, # Use the first (and only) contrast in your model
       status = summa.fit[,"CR_vs_NR"], # CORRECTED: Use the new contrast name
       values = c(-1, 1), # Highlight down-regulated (-1) and up-regulated (1) genes
       hl.col = c("blue", "red") # FIXED: Use clear colors, e.g., blue and red
)



# Highlight the top 100 most differentially expressed (DE) genes
# and label them using the gene symbols.

volcanoplot(
  fit.cont, 
  coef = 1, # Use the first (and only) contrast, CR_vs_NR
  highlight = 100, 
  # CORRECTED: Gene symbols are in 'v$genes$SYMBOL' or 'fit.cont$genes$SYMBOL' if available,
  # but the simplest way to get symbols is often from the original 'y' or 'v' object.
  # Assuming the gene symbols were stored in the 'genes' data frame of the 'v' object:
  names = v$genes$ORF, 
  main = "Calorie Restriction vs. Non-Restricted" # FIXED: Use a descriptive title
)

## gene ontology analysis
library(org.Sc.sgd.db)
library(AnnotationDbi)

fit.cont
View(fit.cont)
# 1. Get the indices of significantly DE genes from the decideTests result
#    'summa.fit' contains 1 (Up), -1 (Down), and 0 (Not DE)
is_de <- which(summa.fit[, "CR_vs_NR"] != 0)
summa.fit
View(summa.fit)

# 2. Extract the SGD ORF IDs for only the DE genes
de_orf_ids <- rownames(fit.cont)[is_de]
de_orf_ids
length(de_orf_ids)

# Check how many DE genes you found
print(paste("Number of Differentially Expressed Genes:", length(de_orf_ids)))


# Columns available in org.Sc.sgd.db (use 'columns(org.Sc.sgd.db)' to see all)
desired_columns <- c(
  "ORF",        # The input key (for confirmation)
  "GENENAME",   # Common Name (Gene Symbol)
  "DESCRIPTION", # Full Description
  "GO"          # Gene Ontology ID
)


# Load the annotation package (already in your session)
library(org.Sc.sgd.db) 
library(AnnotationDbi)

# Retrieve the annotations
de_gene_annotations <- AnnotationDbi::select(
  x = org.Sc.sgd.db,
  keys = de_orf_ids,
  columns = desired_columns,
  keytype = "ORF" # We are searching using the SGD ORF ID
)

# View the structure of the result
head(de_gene_annotations)
dim(de_gene_annotations)
View(all_genes) 
 
de_orf_ids

# 1. Store the Entrez IDs in a clean vector
de_entrez_ids <- de_orf_ids # We established these are Entrez IDs

# 2. Extract the statistical data for only these significant genes from fit.cont
de_stats_table <- topTable(
  fit.cont, 
  coef = "CR_vs_NR", 
  number = length(de_entrez_ids), 
  genelist = de_entrez_ids # Filter to ONLY the significant IDs
)

# 3. Add the Regulation Status column
de_stats_table$Regulation_Status <- ifelse(
  de_stats_table$logFC > 0, 
  "Up-regulated (CR > NR)", 
  "Down-regulated (CR < NR)"
)

# 4. Filter the table to see only the UP-regulated genes
upregulated_genes <- de_stats_table[de_stats_table$Regulation_Status == "Up-regulated (CR > NR)", ]

# 5. Filter the table to see only the DOWN-regulated genes
downregulated_genes <- de_stats_table[de_stats_table$Regulation_Status == "Down-regulated (CR < NR)", ]


View(de_stats_table)
# View the top 10 most up-regulated genes
head(upregulated_genes[order(upregulated_genes$logFC, decreasing = TRUE), 
                       c("logFC", "adj.P.Val", "Regulation_Status")], 10)

# View the top 10 most down-regulated genes
head(downregulated_genes[order(downregulated_genes$logFC, decreasing = FALSE), 
                         c("logFC", "adj.P.Val", "Regulation_Status")], 10)






 



