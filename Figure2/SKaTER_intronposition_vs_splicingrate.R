library(tidyverse)
library(gtools)
library(reshape2)
library(FSA)

setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\SKaTER\\Sequencing_Analysis\\results')

dmso_1 <- read.table('skater_optimize_v0_5_DMSOrep1_database_sorted_intronPosition_boxPlot.txt', header =T, sep = '\t')
dmso_2 <- read.table('skater_optimize_v0_5_DMSOrep2_database_sorted_intronPosition_boxPlot.txt', header = T, sep = '\t')

gsk_1 <- read.table('skater_optimize_v0_5_gsk591rep1_database_sorted_intronPosition_boxPlot.txt', header =T, sep = '\t')
gsk_2 <- read.table('skater_optimize_v0_5_gsk591rep2_database_sorted_intronPosition_boxPlot.txt', header = T, sep = '\t')

ms023_1 <- read.table('skater_optimize_v0_5_ms023rep1_database_sorted_intronPosition_boxPlot.txt', header =T, sep = '\t')
ms023_2 <- read.table('skater_optimize_v0_5_ms023rep2_database_sorted_intronPosition_boxPlot.txt', header = T, sep = '\t')

dmso <- rbind(dmso_1, dmso_2)
gsk <- rbind(gsk_1, gsk_2)
ms023 <- rbind(ms023_1, ms023_2)

dataframes <- list(dmso = dmso, gsk = gsk, ms023= ms023)
df <- bind_rows(dataframes, .id='df')
df.m <- reshape2::melt(df)
df.m <- na.omit(df.m)
colors <- c('#a3a3a3', '#008837','#7b3294')
#define order of graph
order <- c('dmso', 'gsk', 'ms023')
#convert dataframe names to factors to allow ordering
df.m$df <- factor(df.m$df, levels = order)
#Alyssa's script uses log2 so need to convert back to number before reconverting to log10
df.m$convertedvalue = 2 ^ df.m$value

stat_box_data <- function(y, upper_limit = log10(max(df.m$convertedvalue)) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 3), '\n')
    )
  )
}

#uncomment directly below to save pdf
#pdf('../../../R_figures/intronPosition_vs_SplicingRate.pdf', width = 8, height = 8)
ggplot(df.m, aes(x=variable, y=log10(convertedvalue), fill = df)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = colors) +
  #geom_smooth(method=lm,   # Add linear regression lines
  #            se=TRUE, fullrange = TRUE, color = 'black', fill = 'dimgray') +
  #geom_point(aes(fill = df), size =0.5, shape = 21, position = position_jitterdodge()) +
  #geom_jitter(position = position_dodge2(width = 1), color="black", size=0.1, alpha=0.1) +
  ggtitle("") +
  xlab("") + 
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )+
  ylab('') 
#dev.off()
kruskal.test(value~variable, data = df)
pairwise.wilcox.test(df$value, df$variable,
                     p.adjust.method = "BH")
dunnTest(value~variable, df.m, method='bh')
