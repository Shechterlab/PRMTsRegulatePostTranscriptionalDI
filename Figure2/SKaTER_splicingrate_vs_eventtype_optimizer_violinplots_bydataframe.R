library(tidyverse)
library(reshape2)
library(RColorBrewer)

setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\')


#load in dataframes
d2g <- read.table('SKaTER/Sequencing_Analysis/results/skater_optimize_v0_5_GSK591rep1_database_sorted_nDB=2_splicingRates_eventType_eventType_boxPlot.txt', header = T, row.names = NULL, sep = '\t', skip =1)
d2m <- read.table('SKaTER/Sequencing_Analysis/results/skater_optimize_v0_5_MS023rep1_database_sorted_nDB=2_splicingRates_eventType_eventType_boxPlot.txt', header = T, row.names = NULL, sep = '\t', skip =1)
dmso <- read.table('SKaTER/Sequencing_Analysis/results/skater_optimize_v0_5_DMSOrep1_database_sorted_nDB=2_splicingRates_eventType_eventType_boxPlot.txt', header = T,row.names = NULL, sep = '\t', skip =1)

#compile a list of dataframes

dataframes <- list(d2g = d2g, d2m = d2m, dmso = dmso)
dataframes <- lapply(dataframes, function(x) {colnames(x) <- c('A5SS', 'constitutive intron', 'casette exon', 'A3SS');x})
dataframes <- bind_rows(dataframes, .id='df')
df <- reshape2::melt(dataframes)
df <- na.omit(df)

colors <- c('#D4D4D4','#008837','#7b3294')
colors <- brewer.pal(4, 'Spectral')
#define order of graph
order <- c( 'dmso','d2g', 'd2m')
#convert dataframe names to factors to allow ordering
df$df <- factor(df$df, levels = order)

stat_box_data <- function(y, upper_limit =max(df$value) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 3), '\n')
    )
  )
}


#pdf('R_figures/skater_violinplot_eventype_splicingrate_optimizer_dmso.pdf', height = 10, width = 10)
ggplot(df[df$df == 'dmso',], aes(x=variable, fill=variable, y=value))+
  scale_fill_manual(values=colors)+
  geom_violin(alpha=0.5, orientation = 'x', width=1, position=position_dodge(width=0.5))+
  #geom_jitter(aes(x=group, y=rate), shape =21)+
  geom_boxplot(width=0.1, fatten=5)+
  geom_hline(yintercept= median(df[df$df=='dmso' & df$variable == 'constitutive intron',]$value),size=1, linetype ='dotted') +
  scale_color_manual(values=colors) +
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  labs(title="", y = "log10(Splicing Rate)", x="")+
  ylim(-10,5)+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )
#dev.off()

ks.test(df[df$df=='dmso' & df$variable == 'A5SS',]$value, df[df$df=='dmso' & df$variable == 'constitutive intron',]$value, alternative = 'greater')
ks.test(df[df$df=='dmso' & df$variable == 'A5SS',]$value, df[df$df=='dmso' & df$variable == 'casette exon',]$value)
ks.test(df[df$df=='dmso' & df$variable == 'A5SS',]$value, df[df$df=='dmso' & df$variable == 'A3SS',]$value)
ks.test(df[df$df=='dmso' & df$variable == 'constitutive intron',]$value, df[df$df=='dmso' & df$variable == 'casette exon',]$value)
ks.test(df[df$df=='dmso' & df$variable == 'constitutive intron',]$value, df[df$df=='dmso' & df$variable == 'A3SS',]$value)
ks.test(df[df$df=='dmso' & df$variable == 'casette exon',]$value, df[df$df=='dmso' & df$variable == 'A3SS',]$value)


t.test(df[df$variable == 'd4c',]$value, df[df$variable =='d2g',]$value)

