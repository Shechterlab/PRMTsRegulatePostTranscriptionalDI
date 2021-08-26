library(tidyverse)
library(reshape2)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\')


#load in dataframes
d2g <- read.table('SKaTER/Sequencing_Analysis/histograms/GSK591rep1_GSK591rep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')
d2m <- read.table('SKaTER/Sequencing_Analysis/histograms/MS023rep1_MS023rep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')
dmso <- read.table('SKaTER/Sequencing_Analysis/histograms/DMSOrep1_DMSOrep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')

#compile a list of dataframes

dataframes <- list(d2g = d2g, d2m = d2m, d4c = dmso)
dataframes <- dataframes %>% reduce(inner_join, by = 'V2')
colnames(dataframes)[5] <- 'd2g'
colnames(dataframes)[9] <- 'd2m'
colnames(dataframes)[13] <- 'dmso'
dataframes <- dataframes[c(1,5,9,13)]

df <- reshape2::melt(dataframes)


colors <- c('#D4D4D4','#008837','#7b3294')

#define order of graph
order <- c( 'dmso','d2g', 'd2m')
#convert dataframe names to factors to allow ordering
df$variable <- factor(df$variable, levels = order)

stat_box_data <- function(y, upper_limit = max(1) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 3), '\n')
    )
  )
}


#pdf('R_figures/skater_violinplot_commongenes_rep_splicingrate.pdf', height = 10, width = 10)
ggplot(df, aes(x=variable, fill=variable, y=log10(value)))+
  scale_fill_manual(values=colors)+
  geom_violin(alpha=0.5, orientation = 'x', width=1, position=position_dodge(width=0.5))+
  #geom_jitter(aes(x=group, y=rate), shape =21)+
  geom_boxplot(width=0.1, fatten=5)+
  geom_hline(yintercept= median(log10(df[df$variable=='dmso',]$value)),size=1, linetype ='dotted') +
  scale_color_manual(values=c('#afadae','#f2a182')) +
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  labs(title="", y = "log10(Splicing Rate)", x="")+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )
#dev.off()

ks.test(dataframes$dmso, dataframes$d2g)
ks.test(dataframes$dmso, dataframes$d2m)
ks.test(dataframes$d2m, dataframes$d2g)
wilcox.test(dataframes$dmso, dataframes$d2g)
wilcox.test(dataframes$dmso, dataframes$d2m)
wilcox.test(dataframes$d2g, dataframes$d2m)