install.packages('ggpubr')
install.packages('openxlsx')
install.packages('BiocManager')
BiocManager::install('TCGAbiolinks')
library(airway)
data("airway")

library(tidyverse)
library(GEOquery)
library(ggpubr)
library(openxlsx)
library(naniar)
#library(DESeq2)
#library(airway)


#dataManipualtion
#selection,filter,mutate,grouping and summarizing

cancer = read.csv('rna_seq_dataset_30genes.csv')
view(cancer)
cancer$tissue_type="breast tumor"
view(cancer)
dim(cancer)


#subsetting column
cancer=cancer[,c(1,3,5)]
view(cancer)
 

#is there any missing value
is.na(cancer)
sum(is.na(cancer))
miss_var_summary(cancer)
gg_miss_var(cancer)
miss_var_which(cancer)

#column selection
select(cancer,gene)
select(cancer,c(1,2))

#filter
#filtering data using(==)

filter(cancer, metastasis=='no')
filter(cancer,fpkm>10)

#filtering using '&' operator
filter(cancer, fpkm>10 & fpkm< 50)

#filtering using '|' operator
filter(cancer, fpkm>10 | fpkm< 50)


#select and filter
#we can also do this step by step by callin new vaeriable each time
#chaining method using pipe operator
cancer |> 
select(gene,fpkm,metastasis) |> 
         filter(metastasis =='yes')

#multiple filtering criteria
cancer |> 
  filter(gene %in% c('BRCA1','BRCA2','TP53','ALK'))

#MUTATE(CREATING A NEW COLUMN)
cancer |> 
  mutate(fpkm_log=log(fpkm))



