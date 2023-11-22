

rm(list = ls())
library(ggpubr)
source('plotGOFromGeneOntology.R')
df1 <- plotGOFromGeneOntology(file = '../../data/GO_data/GO_roXbound_autosomal_genes.txt',
                              maxNumberOfGOterms = 10,returnDataType = 'dataframe',
                              sortby = 'FDR',showname = 'term')
p <- ggplot(df1,aes(y=GOterms,x=-log10(FDR))) +
  geom_point(aes(size=GeneList)) +
  labs(y=NULL,x=-log[10]~FDR,size='gene number',color=expression(-log[10]~FDR)) +
  theme(plot.background = element_blank(),panel.background = element_rect(fill = 'transparent',color = 'black'),panel.grid = element_blank(),
        legend.position = 'right',legend.background = element_blank(),legend.key = element_blank(),legend.box = 'vertical',
        axis.text = element_text(color = 'black',size = 10),axis.title = element_text(size = 11)) +
  scale_color_gradient()
print(p)

