library(tidyverse)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\')
#load in dataframes
d2g <- read.table('SKaTER/Sequencing_Analysis/histograms/GSK591rep1_GSK591rep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')
d2m <- read.table('SKaTER/Sequencing_Analysis/histograms/MS023rep1_MS023rep2_splicing_expCutoff=5_histogram.txt', header = F, sep = '\t')
#compile a list of dataframes
dataframes <- list(d2g = d2g, d2m = d2m)
dataframes_coordinates <- lapply(dataframes, function(x) str_split_fixed(x$V2, ",", 2))
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
df_RI <- lapply(df_RI, function(x) {x <- x[c(1,7,6,30)]})
df_RI <- lapply(df_RI, setNames, c("upstreamEE","downstreamES","rate","group"))

df <- lapply(df, function(x) {x$group <- 'global';x})
df <- lapply(df, function(x) {x <- x[c(6,7,5,8)]})
df <- lapply(df, setNames, c("upstreamEE","downstreamES","rate","group"))

df <- bind_rows(df, .id='df')
df_RI <- bind_rows(df_RI, .id = 'df')

df <- rbind(df, df_RI)

df <- df[df$df == 'd2g',]

stat_box_data <- function(y, upper_limit = log10(max(df$rate)) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 3), '\n')
    )
  )
}


#pdf('R_figures/d2g_SKaTER_splicingrate_vs_RI.pdf', height = 10, width = 10)
ggplot(df, aes(x=group, fill=group, y=log10(rate)))+
  scale_fill_manual(values=c('#afadae','#f2a182'))+
  geom_violin(alpha=0.5, orientation = 'x', width=1, position=position_dodge(width=0.5))+
  #geom_jitter(aes(x=group, y=rate), shape =21)+
  geom_boxplot(width=0.1, fatten=5)+
  geom_hline(yintercept= median(log10(df[df$group == 'global',]$rate)),size=1, linetype ='dotted') +
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

#randomly shuffle and permute data to determine p-value
set.seed(101) ## for reproducibility
nsim <- 1000 ##number of permutations
res <- numeric(nsim)
res2 <-numeric(nsim)## set aside space for results
n = nrow(df[df$group == 'RI',]) ##define number of rows to sample
for (i in 1:nsim) {
  ## standard approach: scramble response value
  random <- sample_n(df[df$group == 'global',], n)
  test <- wilcox.test(df[df$group == 'RI',]$rate, random$rate, alternative = 'less')
  ## compute & store difference in means; store the value
  res[i] <- test$p.value
  res2[i] <- test$statistic
}
## append the observed value to the list of results
median(res)
median(res2)