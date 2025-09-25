countdata <- read.csv('Mouse expression dataset.csv',row.names = 1)
View(countdata)
# 1. Check if ANY duplicates exist in the column:
any(duplicated(countdata$GENE))

# 2. Get the actual duplicate values:
countdata$GeneID[duplicated(countdata$GENE)]

# 3. Get the count of how many unique entries have duplicates:
sum(duplicated(countdata$GENE))

# 1. Identify all values that appear more than once in the "GeneID" column.
#    The 'duplicated()' function marks the 2nd, 3rd, etc. occurrences as TRUE.
#    Reversing the vector and applying duplicated() marks the 1st occurrence as TRUE.
#    The logical OR (|) combines both results to flag ALL duplicate values.
is_duplicate <- duplicated(countdata$GENE) | duplicated(countdata$GENE, fromLast = TRUE)

# 2. Filter the data frame to keep only the rows where 'is_duplicate' is FALSE (i.e., unique entries).
countdata_unique <- countdata[!is_duplicate, ]

# 3. Check the new dimensions
dim(countdata_unique)
View(countdata_unique)
any(duplicated(countdata_unique$GENE))
sum(duplicated(countdata_unique$GENE))

countdata_unique<- countdata_unique[,c(2,3,4,5,6,7)]
head(countdata_unique)
View(countdata_unique)
#Now let me create expermental labels (two conditions)
coldata <- DataFrame(condition = factor(c("ctrl","ctrl","ctrl","treat",
                                          "treat","treat")))
coldata

library(DESeq2) 

#Now lets create a DESeq input matrix
dds <- DESeqDataSetFromMatrix(countData = countdata_unique, colData = coldata,
                              design = ~condition)
dds


##Now lets run DEseq
dds <- DESeq(dds)
#Lets get differentially expressed genes
res <- results(dds)

#Lets look at the table
head(res)


