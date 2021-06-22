library(tidyverse)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT')
#load in all the rMATS analyses
d2g <- read.table('SKaTER/Sequencing_Analysis/histograms/fractionCleavageBeforeSplicing/GSK591rep1_GSK591rep2_fractionCleavageBeforeSplicing.txt', header = T, sep = '\t')
d2m <- read.table('SKaTER/Sequencing_Analysis/histograms/fractionCleavageBeforeSplicing/MS023rep1_MS023rep2_fractionCleavageBeforeSplicing.txt', header = T, sep = '\t')
#dmso <- read.table('SKaTER/Sequencing_Analysis/histograms/fractionCleavageBeforeSplicing/DMSOrep1_DMSOrep2_fractionCleavageBeforeSplicing.txt', header = F, sep = '\t')

#compile a list of dataframes
dataframes <- list(d2g = d2g, d2m = d2m)
#dataframes <- lapply(dataframes, function(x) {x <- x[c(1, 6)]})
#colnames <- c('gene', 'rate')
#dataframes <- lapply(dataframes, setNames, colnames)
dataframes_coordinates <- lapply(dataframes, function(x) str_split_fixed(x$intron, ",", 2))
dataframes_coordinates <- lapply(dataframes_coordinates, function(x) as.data.frame(gsub("[()]", "", x)))
dataframes_coordinates <- lapply(dataframes_coordinates, setNames, c("upstreamEE", "downstreamES"))
df <- mapply(cbind, dataframes, dataframes_coordinates, SIMPLIFY=F)

#df <- lapply(df, function(x) {x$upstreamEE <- gsub("\"", "", x$upstreamEE);x})
df <- lapply(df, function(x) {x$upstreamEE <- strtoi(x$upstreamEE);x})
df <- lapply(df, function(x) {x$downstreamES <- strtoi(x$downstreamES);x})

d2g <- read.table('RNAseq/rMATS_4.0.2/rmats4.0.2_D4C_vs_D2G/RI.MATS.JCEC.txt', header = T, sep = '\t')
d2m <- read.table('RNAseq/rMATS_4.0.2/rmats4.0.2_D4C_vs_D2M/RI.MATS.JCEC.txt', header = T, sep = '\t')
df_RI <- list(d2g = d2g, d2m = d2m)
df_RI <- lapply(df_RI, function(x) {x <- x[x$FDR < 0.05,];x})

df_RI <- Map(merge, df, df_RI, by=c("upstreamEE","downstreamES"))
df_RI <- lapply(df_RI, function(x) {x$group <- 'RI';x})
df_RI <- lapply(df_RI, function(x) {x <- x[c(1,8,7,31)]})
df_RI <- lapply(df_RI, setNames, c("upstreamEE","downstreamES","rate","group"))

df <- lapply(df, function(x) {x$group <- 'global';x})
df <- lapply(df, function(x) {x <- x[c(7,8,6,9)]})
df <- lapply(df, setNames, c("upstreamEE","downstreamES","rate","group"))

df <- bind_rows(df, .id='df')
df_RI <- bind_rows(df_RI, .id = 'df')

df <- rbind(df, df_RI)

df <- df[df$df == 'd2g',]

stat_box_data <- function(y, upper_limit = max(df$rate) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 3), '\n')
    )
  )
}


#pdf('R_figures/d2m_SKaTER_probabilityofcleavage_vs_RI_violin.pdf', height = 10, width = 10)
ggplot(df, aes(x=group, fill=group, y=rate))+
  scale_fill_manual(values=c('#afadae','#f2a182'))+
  geom_violin(alpha=0.5, orientation = 'x', width=1, position=position_dodge(width=0.5))+
  #geom_jitter(aes(x=group, y=rate), shape =21)+
  geom_boxplot(width=0.1, fatten=5)+
  geom_hline(yintercept= median(df[df$group == 'global',]$rate),size=1, linetype ='dotted') +
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

ks.test(df[df$group == 'global',]$rate, df[df$group == 'RI',]$rate, alternative = 'greater')
wilcox.test(df[df$group == 'global',]$rate, df[df$group == 'RI',]$rate)

