#plot the chromosome distribution of roX peaks.
rm(list = ls())
library(ggplot2)
library(ggsci)
library(reshape2)
fi <- read.table('../../data/peak_distribution_data/chromosome_distribution_genome_roXpeaks.txt',sep='\t',header=T,check.names = F)
for(i in 2:3){
  #transform into percentage
  fi[,i] <- fi[,i]/sum(fi[,i]) *100
}
dat <- melt(fi,variable.name = 'factor',value.name = 'number')
dat$Chromosome <- factor(dat$Chromosome,levels = unique(dat$Chromosome))
p <- ggplot(dat,aes(x=factor,y=number,fill=Chromosome)) +
  geom_bar(stat = 'identity',position = 'stack') +
  labs(x=NULL,y='Distribution (%)',title = NULL) +
  #scale_fill_aaas() +
  scale_fill_npg() +
  theme(axis.text = element_text(color = 'black'),
        plot.background = element_blank(),legend.background = element_blank(),legend.key = element_blank(),
        panel.background = element_rect(fill = 'transparent',color = 'black'),panel.grid = element_blank(),
        legend.position = 'right',legend.title = element_blank())
print(p)
