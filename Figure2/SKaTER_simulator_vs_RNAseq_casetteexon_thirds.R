library(tidyverse)
library(gtools)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\')



#load in dataframes
dmso <- readxl::read_xlsx('SKaTER/Sequencing_Analysis/rMATS_simulator/00_skater_optimize_v0_5_DMSOrep1_DMSOrep2_MERGE_v2_eCutoff=10.0_frac=1_eventUnion_psi.xlsx')
d2g <- readxl::read_xlsx('SKaTER/Sequencing_Analysis/rMATS_simulator/00_skater_optimize_v0_5_GSK591rep1_GSK591rep2_MERGE_v2_eCutoff=10.0_frac=1_eventUnion_psi.xlsx')
d2m <- readxl::read_xlsx('SKaTER/Sequencing_Analysis/rMATS_simulator/00_skater_optimize_v0_5_MS023rep1_MS023rep2_MERGE_v2_eCutoff=10.0_frac=1_eventUnion_psi.xlsx')


#compile a list of dataframes
dataframes <- list(d2m = d2m, d2g=d2g, dmso= dmso)

dataframes <- lapply(dataframes, function(x) {x$sampleavg <- (x$sample1avg + x$sample2avg)/2;x})

dataframes <- lapply(dataframes, function(x) {x$thirds <- quantcut(x$avg_simPSI, q=3, na.rm=TRUE);x})

df <- bind_rows(dataframes, .id='df')

df <- as.data.frame(df)

df$thirds <- as.character(df$thirds)

df$group <- NA

for (i in 1:nrow(df)) {
  if (df[i,11] == "[0.0344,0.75]" | df[i,11] == "[0.0946,0.679]" | df[i,11] == '[0.0276,0.712]') {
    df[i,18] <- 'bottom'
  } 
  else if (df[i,11] == "(0.75,0.979]" | df[i,11] == '(0.679,0.982]' | df[i,11] == "(0.712,0.975]") {
    df[i,18] <- "middle"
  } 
  else {
    df[i,18] <- "top"
  }
}



#colors: 
d2g_col <- c('#108745')
d2m_col <- c('#7d3593')
d4c_col <- c('#d4d4d4')

colors <- c('#D4D4D4','#008837','#7b3294')
#define order of graph
order <- c( 'dmso','d2g', 'd2m')
#convert dataframe names to factors to allow ordering
df$df <- factor(df$df, levels = order)
#uncomment directly below to save pdf
#pdf('R_figures/SKaTER_SE_simulator_vs_RNAseq_dmso.pdf', width = 8, height = 8)
#ggplot(dataframes$d4c, aes(x=thirds, y=sampleavg, fill = d4c_col)) +
#  geom_boxplot(alpha = 0.8) +
#  scale_fill_manual(values = d4c_col) +
#  geom_smooth(method=lm,   # Add linear regression lines
#              se=TRUE, fullrange = TRUE, color = 'black', fill = 'dimgray') +
  #geom_point(aes(fill = df), size =0.5, shape = 21, position = position_jitterdodge()) +
  #geom_jitter(position = position_dodge2(), color="black", size=0.1, alpha=0.1) +
#  ggtitle("") +
#  xlab("SKaTER PSI") + 
#  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
#  theme(panel.background = element_blank()) +
#  theme(panel.grid.minor = element_blank()) +
#  ylab('RNA-seq PSI')
#dev.off()

#pdf('R_figures/SKaTER_SE_simulator_vs_RNAseq_all.pdf', width = 8, height = 8)
ggplot(df, aes(x=group, y=sampleavg, fill = df)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = colors) +
  geom_smooth(method=lm,   # Add linear regression lines
              se=TRUE, fullrange = TRUE, color = 'black', fill = 'dimgray') +
  #geom_point(aes(fill = df), size =0.5, shape = 21, position = position_jitterdodge()) +
  #geom_jitter(position = position_dodge2(), color="black", size=0.1, alpha=0.1) +
  ggtitle("") +
  xlab("SKaTER PSI") + 
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylab('RNA-seq PSI')
#dev.off()

cor.test(dataframes$dmso$sampleavg, dataframes$dmso$avg_simPSI, method = 'spearman')
cor.test(dataframes$d2g$sampleavg, dataframes$d2g$avg_simPSI, method = 'spearman')
cor.test(dataframes$d2m$sampleavg, dataframes$d2m$avg_simPSI, method = 'spearman')
