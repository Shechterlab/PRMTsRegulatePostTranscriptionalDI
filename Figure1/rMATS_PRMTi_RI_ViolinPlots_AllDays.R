library(tidyverse)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\RNAseq')
#load in all the rMATS analyses
d2g_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D2G/RI.MATS.JCEC.txt', header = T, sep = '\t')
d2gm_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D2GM/RI.MATS.JCEC.txt', header = T, sep = '\t')
d2m_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D2M/RI.MATS.JCEC.txt', header = T, sep = '\t')
d4g_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D4G/RI.MATS.JCEC.txt', header = T, sep = '\t')
d4gm_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D4GM/RI.MATS.JCEC.txt', header = T, sep = '\t')
d4m_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D4M/RI.MATS.JCEC.txt', header = T, sep = '\t')
c2_RI <- read.table('rMATS_4.0.2/rmats4.0.2_c1_vs_c2/RI.MATS.JCEC.txt', header = T, sep = '\t')
c3_RI <- read.table('rMATS_4.0.2/rmats4.0.2_c1_vs_c3/RI.MATS.JCEC.txt', header = T, sep = '\t')
#compile a list of RI
RI <- list(d2g_RI = d2g_RI, d2gm_RI = d2gm_RI, d2m_RI = d2m_RI, d4g_RI = d4g_RI, d4gm_RI = d4gm_RI, d4m_RI = d4m_RI, c2_RI = c2_RI, c3_RI = c3_RI)
#ensure all genes being analyzed have FDR < 0.05 and no duplicate genes; multiply ONLY RI by -1
RI <- lapply(RI, function(x) filter(x, FDR < 0.05))
#RI <- lapply(RI, function(x) x[!duplicated(x[c('geneSymbol')]),])
RI <- lapply(RI, function(x) {x$IncLevelDifference <- x$IncLevelDifference * -1;x})
RI <- lapply(RI, function(x) {x$group <- 'RI';x})

df <- bind_rows(RI, .id = 'df')

for (i in unique(df$df)) {
  print(i)
  print(nrow(df[df$df == i,]))
}

order <- c("c3_RI","c2_RI","d4m_RI","d4gm_RI","d4g_RI","d2m_RI","d2gm_RI","d2g_RI")

df$df <- factor(df$df,  order)
#colors <- c('#1eb266','#35afab','#a377b2','#53975e','#037a77','#955da3','#108745','#7d3593')
#for horizontal plot with d2g at the top
colors<-c('#7d3593','#108745','#955da3','#037a77','#53975e','#a377b2','#35afab','#1eb266')


##pdf('../R_figures/rMATS_RI_violinplot_alldays_boxplot_v2.pdf',height=6, width=4)
ggplot(df[df$group =='RI',], aes(x=IncLevelDifference, y=df,  fill=df)) + 
  geom_violin(alpha=0.5, orientation = 'y')+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values=colors)+
  labs(title="",x="", y = "")+
  #geom_jitter(shape=21,width = .2, alpha = .2) +
  geom_vline(xintercept = 0 , linetype = 'dashed')+
  scale_color_manual(values=colors)+
  theme_classic()+
  theme(legend.position = 'none')+
  xlim(-1,1)
##dev.off()

ks.test(df[df$df =='d2g_RI',]$IncLevelDifference, df[df$df =='d2gm_RI',]$IncLevelDifference)
ks.test(df[df$df =='d2g_RI',]$IncLevelDifference, df[df$df =='d2m_RI',]$IncLevelDifference)
ks.test(df[df$df =='d2gm_RI',]$IncLevelDifference, df[df$df =='d2m_RI',]$IncLevelDifference)

ks.test(df[df$df =='d4g_RI',]$IncLevelDifference, df[df$df =='d4gm_RI',]$IncLevelDifference)
ks.test(df[df$df =='d4g_RI',]$IncLevelDifference, df[df$df =='d4m_RI',]$IncLevelDifference)
ks.test(df[df$df =='d4gm_RI',]$IncLevelDifference, df[df$df =='d4m_RI',]$IncLevelDifference)

ks.test(df[df$df =='c2_RI',]$IncLevelDifference, df[df$df =='c3_RI',]$IncLevelDifference)
