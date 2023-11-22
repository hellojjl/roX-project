#######
##re-plot the binding profiles from deeptools output

rm(list = ls())
source('plotProfileFromDeeptools.R')
library(ggpubr)

# plot the epigenome profile in the vicinity of roX binding sites
df <- plotProfile(infile = '../../data/epigenome_profile_data/epigenome_roX_peakcenter_chrtypes.tab',legend.position = 'bottom',regionlabel = c('roX.peaks.autosomes.bed'='autosomes','roX.peaks.chrX.bed'='chrX'),
                  returnDataType = 'dataframe',ignore.peaktypes = 'chrY',axis.textsize = 8,ytitle = 'Enrichment',nrow = 1,changeFacet = F,xyscales = 'free')

#'H4K16ac','H3K27me3','H3K36me3' only
df2 <- subset(df,df$sample %in% c('H4K16ac','H3K27me3','H3K36me3'))
p <- ggplot(df2,aes(x=binregion,y=enrichment,color=chrtype)) +
  geom_line() +
  labs(x=NULL,y='Enrichment') +
  facet_wrap(~sample,nrow = 1,scales = 'free') +
  theme(plot.background = element_blank(),panel.background = element_rect(color = 'black',fill = 'transparent'),
        panel.grid = element_blank(),strip.background = element_blank(),axis.text = element_text(color = 'black',size = 6),
        axis.title = element_text(size = 6),
        legend.position = 'bottom',legend.background = element_blank(),legend.title = element_blank(),
        legend.box.background = element_blank(),legend.key = element_rect(fill = 'transparent')) +
  scale_x_continuous(breaks = c(0,30,60),labels = c('-3 kb','center','3 kb'))
print(p)

############# plot PRC enrichment on roX binding sites ###########
rm(list = ls())
source('plotProfileFromDeeptools.R')
library(ggpubr)

df <- plotProfile(infile = '../../data/PRC_enrichment_data/PRC_signals_roXpeakcenter_chrtypes_sortedbyroX.tab',replaceSymbol = "()",
                  regionlabel = c('Heatmap1sortedRegions_autosomes.bed'='autosomes','Heatmap1sortedRegions_chrX.bed'='chrX'),
                  returnDataType = 'dataframe')
dat <- subset(df,df$sample %in% c('Pc','Ph','Psc','E(z)','Su(z)12','Pho','Spps'))
p <- ggplot(dat,aes(x=binregion,y=enrichment,color=chrtype)) +
  geom_line() +
  facet_wrap(~sample,nrow = 2,scales = 'free_y') +
  labs(x=NULL) +
  scale_x_continuous(breaks = c(0,30,60),labels = c('-3kb','center','3kb')) +
  theme(plot.background = element_blank(),panel.background = element_rect(color = 'black',fill = 'transparent'),
        panel.grid = element_blank(),strip.background = element_blank(),axis.text = element_text(color = 'black',size = 9),
        axis.title = element_text(size = 11),
        legend.position = c(0.9,0.2),legend.background = element_blank(),legend.title = element_blank(),
        legend.box.background = element_blank(),legend.key = element_rect(fill = 'transparent'))
print(p)


