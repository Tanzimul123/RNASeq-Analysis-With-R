library(tidyverse)
library(GEOquery)
library(ggpubr)
library(openxlsx)
library(naniar)


#importing raw counts data
cancer= read.csv("GSE183947_fpkm.csv")
View(cancer)
ncol(cancer)


#get metadata

result = getGEO(GEO = "GSE183947",GSEMatrix = T)

result
class(result)

#metadata#this code is needed to be written in every case
metadata <- pData(phenoData(result[[1]]))
view(metadata)

#subsetting
metadata=metadata[,c(1,10,11,17)]
view(metadata)

#renamingColumns
colnames(metadata)<-c('title','tissue_type',"metastatic_state",'description')
view(metadata)

metadata_new=metadata |> 
  mutate(tissue_type = gsub('tissue: ','',tissue_type) ) 
view(metadata_new)

metadata_new1=metadata_new |> 
  mutate(metastatic_state= gsub('metastasis: ','',metastatic_state))

view(metadata_new1)


#reshaping cancer data from horzontal to vertical 
cancer_new=cancer |> 
  rename(gene= X) |> 
  pivot_longer(-gene,
               names_to = 'samples',
               values_to = 'fpkm')#(-gene excludes that particular column)

view(cancer_new)

#joining data
#actualRNASeqData
counts_final_data= cancer_new |> 
  left_join(metadata_new1, by =c('samples'='description'))

view(counts_final_data) 

#export_data
write.csv(counts_final_data,'GSE183947_counts_matrix.csv',row.names = F)#row.names=F>ENSURES NO NEW NUMBERING OF ROWS ARE DONE

#Grouping and summarizing


#metastatic state is a group
#tissue type is a group
#how the fpkm of a partcular gene  is differing by tissue
#TISSUE_TYPE & GENE BECAUSE WE HAVE TWO VARIABLES
#mean_fpkm is written just to add a column name in the console

counts_final_data |> 
  filter(gene == 'BRCA1' | gene=='TSPAN6') |> 
  group_by(tissue_type,gene)  |> 
  summarise(mean_fpkm = mean(fpkm),median_fpkm=
              median(fpkm))


#arrange

counts_final_data |> 
  filter(gene == 'BRCA1' | gene=='TSPAN6') |> 
  group_by(tissue_type,gene)  |> 
  summarise(mean_fpkm = mean(fpkm),median_fpkm=
              median(fpkm)) |> 
  arrange(desc(mean_fpkm))














