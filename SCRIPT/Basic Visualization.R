library(tidyverse)
#import data
library(naniar)
cancer= read.csv('GSE183947_counts_matrix.csv')

#bar chart
cancer |> 
  filter(gene == 'BRCA1') |> 
  ggplot(aes(x= samples ,y= fpkm,fill = tissue_type)) + geom_col()


#boxplot
miss_var_which(cancer)
cancer$metastatic_state <- factor(cancer$metastatic_state, levels = c('yes','no'))
cancer |> 
  filter(gene == 'BRCA1') |> 
  ggplot(aes(x = metastatic_state,y= fpkm,fill = tissue_type)) + geom_boxplot()

#VIOLOIN PLOT(WHEN OBSERVATIONS ARE LARGE)

cancer |> 
  filter(gene == 'BRCA1') |> 
  ggplot(aes(x = metastatic_state,y= fpkm,fill = tissue_type)) + geom_violin()

#HISTOGRAM
cancer |> 
  filter(gene == 'BRCA1') |> 
  ggplot(aes(x = fpkm,fill = tissue_type)) + geom_histogram() + facet_wrap(~tissue_type)

#density plot
cancer |> 
  filter(gene == 'BRCA1') |> 
  ggplot(aes(x = fpkm,fill = tissue_type)) + geom_density() + facet_wrap(~tissue_type)


#scatter plot
#we will take two genes and try to understand the relationship between them
#As BRCA1 and BRCA2 are not column/variables,we have to creat them as variable by spread function
cancer |> 
  filter(gene=='BRCA1'| gene=='BRCA2') |> 
  spread(key = gene,value = fpkm) |> 
  ggplot(aes(x= BRCA1, y = BRCA2,color = tissue_type))+
  geom_point()

#CORRELATION(ADD STATS)

cancer |> 
  filter(gene=='BRCA1'| gene=='BRCA2') |> 
  spread(key = gene,value = fpkm) |> 
  ggplot(aes(x= BRCA1, y = BRCA2,color = tissue_type))+
  geom_point()+
  geom_smooth(method = 'lm',se = F)
#IN NORMAL TISSUE THERE IS A POSITIVE CORRELATION BETWEEN brca1 AND brca2
#IN tumor tissue there is negative correlation between brca1 and brca2

#HEATMAP
gene_of_interest<-c('BRCA1','BRCA2','TP53','ALK','MYCN')
cancer |> 
  filter(gene %in% gene_of_interest ) |>  
  ggplot(aes(x=samples , y= gene,fill = fpkm)) +
  geom_tile()+
  scale_fill_gradient(low = 'white',high = 'red')


