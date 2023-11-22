rm(list = ls())
library(ggpubr)
library(reshape2)
fi <- read.table('../../data/gene_expression_data/roXassociatedgenes_expression.txt',sep = '\t',header = T,check.names = F)
fi <- fi[,c('geneid','genename','chrtype','roX-KO-L1_log2FC','msl1kd_log2FC','msl2kd_log2FC','mofkd_log2FC')] #new
names(fi) <- c('geneid','genename','chrtype','roX-KO','MSL1-KD','MSL2-KD','MOF-KD') #new
#roX only
fi_rox <- na.omit(fi[,c("geneid","chrtype","roX-KO")])
fi_rox$FoldChange <- 2^(fi_rox$`roX-KO`)
pv <- ks.test(x=fi_rox$FoldChange[which(fi_rox$chrtype=='autosomes')],y=fi_rox$FoldChange[which(fi_rox$chrtype=='chrX')])$p.value
pv <- ifelse(pv==0,'p < 2.2e-16',paste('p =',pv))
pa <- ggplot(fi_rox,aes(x=FoldChange)) +
  stat_ecdf(aes(color=chrtype)) +
  labs(x='Fold change roX-KO/control',y='Cumulative fraction') +
  scale_x_continuous(limits = c(0,3)) +
  theme(plot.background = element_blank(),panel.background = element_rect(color = 'black',fill = 'transparent',),
        panel.grid = element_blank(),
        legend.position = 'bottom',legend.title = element_blank(),legend.background = element_blank(),legend.key = element_blank(),
        legend.text = element_text(size = 5),legend.margin = margin(0,0,0,0,'pt'),
        axis.text = element_text(color = 'black',size = 5),
        axis.title = element_text(size = 6)) +
  annotate('text',x=1.8,y=.25,size=1.5,label=pv)
print(pa)

#plot violin: roX-KO, MSL1-KD, MSL2-KD, MOF-KD
df <- melt(fi,variable.name = 'genotype',value.name = 'log2FC')
agg <- aggregate(list(pvalue=df$log2FC),list(chrtype=df$chrtype,genotype=df$genotype),FUN = function(x) signif(wilcox.test(x)$p.value,3) ) #https://statisticsglobe.com/set-column-names-within-aggregate-function-in-r
agg$ypos <- c(5,4.5)
pb <- ggplot(df,aes(x=genotype,y=log2FC,fill=chrtype)) +
  geom_violin(size=.2) +
  geom_boxplot(position = position_dodge(width = .9),width=.2,outlier.size = .001,size=.2) +
  geom_hline(yintercept = 0,linetype=2,size=.2) +
  geom_text(data = agg,aes(x=genotype,y=ypos,label=pvalue),size=1.5,position = position_dodge(0.9)) +
  labs(x=NULL,y=expression(log[2]~FoldChange)) +
  scale_y_continuous(limits = c(-5,5)) +
  theme(plot.background = element_blank(),panel.background = element_rect(fill = 'transparent',color = 'black'),panel.grid = element_blank(),axis.title = element_text(size = 6),
        axis.text = element_text(color = 'black',size = 5),legend.position = 'bottom',legend.background = element_blank(),
        legend.key = element_blank(),legend.title = element_blank(),legend.key.size = unit(3,'mm'),legend.text = element_text(size = 5),legend.margin = margin(0,0,0,0,'pt'))
print(pb)
#plot together
pab <- ggarrange(pa,pb,nrow = 1,align = 'h',widths = c(3.5,6.5))
print(pab)