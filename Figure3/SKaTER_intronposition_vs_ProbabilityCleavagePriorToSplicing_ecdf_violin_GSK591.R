library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(FSA)
#replace PRMTi type to dictate which graphs are made
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT')
#load in all the rMATS analyses
data <- read.table('SKaTER/Sequencing_Analysis/results/GSK591rep1_GSK591rep2_fractionCleavageBeforeSplicing.txt', header = T, sep = '\t')

#compile a list of dataframes
dataframes <- list(data=data)
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

introns <- read.table('SKaTER/Sequencing_Analysis/database/skater_optimize_v0_5_GSK591rep1_database_sorted_skater_optimize_v0_5_GSK591rep2_database_sorted_intronPosition.txt', header = T, sep = '\t')

df_introns <- list(introns=introns)
df_introns <- bind_rows(df_introns, .id= 'df')
df_introns <- melt(df_introns, id=c('df'))
df_introns <- list(data=df_introns)

df_introns_coordinates <- lapply(df_introns, function(x) str_split_fixed(x$value, ",", 2))
df_introns_coordinates <- lapply(df_introns_coordinates, function(x) as.data.frame(gsub("[()]", "", x)))
df_introns_coordinates <- lapply(df_introns_coordinates, setNames, c("upstreamEE", "downstreamES"))
df_introns <- mapply(cbind, df_introns, df_introns_coordinates, SIMPLIFY=F)

#df <- lapply(df, function(x) {x$upstreamEE <- gsub("\"", "", x$upstreamEE);x})
df_introns <- lapply(df_introns, function(x) {x$upstreamEE <- strtoi(x$upstreamEE);x})
df_introns <- lapply(df_introns, function(x) {x$downstreamES <- strtoi(x$downstreamES);x})



df_introns <- Map(merge, df, df_introns, by=c("upstreamEE","downstreamES"))
df_introns <- lapply(df_introns, function(x) {x <- x[c(2,8,1,7,10)]})
df_introns <- lapply(df_introns, setNames, c("gene","upstreamEE","downstreamES","rate","group"))

df <- bind_rows(df_introns, .id = 'df')

stat_box_data <- function(y, upper_limit = max(df$rate) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 3), '\n')
    )
  )
}

#RColorBrewer::display.brewer.all()
colors <- brewer.pal(5, 'Blues')

#define order of graph
order <- c('first', 'second','internal','penultimate','last')
#convert dataframe names to factors to allow ordering
df$group <- factor(df$group, levels = order)

#pdf('R_figures/dmso_SKaTER_probabilityofcleavage_vs_intronPosition_violin.pdf', height = 10, width = 10)
ggplot(df, aes(x=group, fill=group, y=rate))+
  scale_fill_manual(values=colors)+
  geom_violin(alpha=0.5, orientation = 'x', width=1, position=position_dodge(width=0.5))+
  #geom_jitter(aes(x=group, y=rate), shape =21)+
  geom_boxplot(width=0.1, fatten=5)+
  geom_hline(yintercept= median(df$rate),size=1, linetype ='dotted') +
  scale_color_manual(values=colors) +
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  labs(title="", y = "Probability of Cleavage Prior to Splicing", x="")+
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9
  )
#dev.off()

#pdf('R_figures/GSK591_SKaTER_probabilityofcleavage_vs_intronPosition_ecdf.pdf', height = 10, width = 10)
ggplot(df,  aes(x=rate, color=group))+
stat_ecdf(size=1, alpha = 1)+
  scale_color_manual(values=colors) +
  scale_size_manual(values=c(1, 2))+
  geom_hline(data=df, aes(yintercept = 0.5),size=1, linetype ='dotted') +
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  labs(title="Empirical Cumulative \n Density Function", y = "Cumulative Probability", x="Probability of Cleavage Prior to Splicing")
#dev.off()

ks.test(df[df$group == 'first',]$rate, df[df$group == 'last',]$rate, alternative = "less")
wilcox.test(df[df$group == 'global',]$rate, df[df$group == 'RI',]$rate)

kruskal.test(group ~ rate, data=df)
dunnTest(rate ~ group, data=df, method="bh")
