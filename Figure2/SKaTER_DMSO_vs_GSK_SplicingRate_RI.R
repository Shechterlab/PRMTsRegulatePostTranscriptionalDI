library(tidyverse)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\')
#load in dataframes
dmso <- read.table('SKaTER/Sequencing_Analysis/histograms/DMSOrep1_DMSOrep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')
d2g <- read.table('SKaTER/Sequencing_Analysis/histograms/GSK591rep1_GSK591rep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')
#d2m <- read.table('SKaTER/Sequencing_Analysis/histograms/MS023rep1_MS023rep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')
#compile a list of dataframes
dataframes <- list(dmso = dmso, d2g = d2g)
dataframes_coordinates <- lapply(dataframes, function(x) str_split_fixed(x$V2, ",", 2))
dataframes_coordinates <- lapply(dataframes_coordinates, function(x) as.data.frame(gsub("[()]", "", x)))
dataframes_coordinates <- lapply(dataframes_coordinates, setNames, c("upstreamEE", "downstreamES"))
df <- mapply(cbind, dataframes, dataframes_coordinates, SIMPLIFY=F)
df <- df %>% reduce(inner_join, by=c("upstreamEE","downstreamES"))
df <- df[c(1,6,7,5,12)]
colnames(df) <- c('gene','upstreamEE','downstreamES','DMSO_rate','GSK_rate')


#df <- lapply(df, function(x) {x$upstreamEE <- gsub("\"", "", x$upstreamEE);x})
df$upstreamEE <- strtoi(df$upstreamEE)
df$downstreamES <- strtoi(df$downstreamES)

d2g <- read.table('RNAseq/rMATS_4.0.2/rmats4.0.2_D4C_vs_D2G/RI.MATS.JCEC.txt', header = T, sep = '\t')
#d2m <- read.table('RNAseq/rMATS_4.0.2/rmats4.0.2_D4C_vs_D2M/RI.MATS.JCEC.txt', header = T, sep = '\t')

d2g <- d2g[d2g$FDR < 0.05,]
#d2g <- d2g[d2g$IncLevelDifference < 0,]

df_RI <- merge(df, d2g, by=c("upstreamEE","downstreamES"))

df_RI <- df_RI[c(3,4,5)]

df <- reshape2::melt(df_RI)

stat_box_data <- function(y, upper_limit = log10(max(df$value)) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 3), '\n')
    )
  )
}


#pdf('R_figures/dmso_vs_d2g_SKaTER_splicingrate_vs_RI.pdf', height = 10, width = 10)
ggplot(df, aes(x=variable, fill=variable, y=log10(value)))+
  scale_fill_manual(values=c('#d4d4d4','#108745'))+
  geom_violin(alpha=0.5, orientation = 'x', width=1, position=position_dodge(width=0.5))+
  #geom_jitter(aes(x=group, y=rate), shape =21)+
  geom_boxplot(width=0.1, fatten=5)+
  geom_hline(yintercept= median(log10(df[df$variable == 'DMSO_rate',]$value)),size=1, linetype ='dotted') +
  scale_color_manual(values=c('#afadae','#108745')) +
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  labs(title="", y = "log10(Splicing Rate)", x="")+
  ylim(-3,1.5)+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )
#dev.off()
ks.test(df[df$variable == 'DMSO_rate',]$value, df[df$variable == 'GSK_rate',]$value)
wilcox.test(df[df$variable == 'DMSO_rate',]$value, df[df$variable == 'GSK_rate',]$value)
