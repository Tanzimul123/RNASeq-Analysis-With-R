library(tidyverse)
#import data

cancer= read.csv('GSE183947_counts_matrix.csv')

#bar chart
cancer |> 
  filter(gene == 'BRCA1') |> 
  ggplot(aes(x= samples ,y= fpkm,fill = tissue_type)) + geom_col()


#boxplot
miss_var_which(cancer)
cancer |> 
  filter(gene == 'BRCA1') |> 
  ggplot(aes(x = metastatic_state,y= fpkm,fill = tissue_type)) + geom_boxplot()

         