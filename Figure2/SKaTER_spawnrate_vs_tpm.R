library(tidyverse)
library(gtools)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\')

D2G_1 <- read.table('RNAseq\\Kallisto\\D2G1/abundance.tsv', sep='\t', header =T)
D2G_2 <- read.table('RNAseq\\Kallisto/D2G2/abundance.tsv', sep='\t', header =T)
D2G_3 <- read.table('RNAseq\\Kallisto/D2G3/abundance.tsv', sep='\t', header =T)
D2G <- merge(D2G_1, D2G_2, by = 'target_id')
D2G <- merge(D2G, D2G_3, by = 'target_id')

D2M_1 <- read.table('RNAseq\\Kallisto\\D2M1/abundance.tsv', sep='\t', header =T)
D2M_2 <- read.table('RNAseq\\Kallisto/D2M2/abundance.tsv', sep='\t', header =T)
D2M_3 <- read.table('RNAseq\\Kallisto/D2M3/abundance.tsv', sep='\t', header =T)
D2M <- merge(D2M_1, D2M_2, by = 'target_id')
D2M <- merge(D2M, D2M_3, by = 'target_id')


D4C_1 <- read.table('RNAseq\\Kallisto\\D4C1/abundance.tsv', sep='\t', header =T)
D4C_2 <- read.table('RNAseq\\Kallisto/D4C2/abundance.tsv', sep='\t', header =T)
D4C_3 <- read.table('RNAseq\\Kallisto/D4C3/abundance.tsv', sep='\t', header =T)
D4C <- merge(D4C_1, D4C_2, by = 'target_id')
D4C <- merge(D4C, D4C_3, by = 'target_id')

tpm <- list(d2g = D2G, d2m = D2M, d4c = D4C)
tpm <- lapply(tpm, function(x) {x$avgtpm <- rowMeans(x[c('tpm.x', 'tpm.y', 'tpm')], na.rm=TRUE);x})
tpm <- lapply(tpm, function(x) {x <- x[c(1, 14)]})


tpm_names <- lapply(tpm, function(x) data.frame(do.call('rbind', strsplit(as.character(x$target_id), '|', fixed = TRUE))))
tpm_names <- lapply(tpm_names, function(x) {x <- x[c(6)]})

tpm <- Map(cbind, tpm, tpm_names)
tpm <- lapply(tpm, function(x) {x <- x[c(3, 2)]})
colnames <- c('gene', 'tpm')
tpm <- lapply(tpm, setNames, colnames)
tpm <- lapply(tpm, function(x) aggregate(tpm ~ gene, data=x, FUN=sum))
tpm <- lapply(tpm, function(x) {x <- x[x$tpm > 1,];x})

#load in all the rMATS analyses
d2g <- read.table('SKaTER/Sequencing_Analysis/histograms/GSK591rep1_GSK591rep2_spawn_expCutoff=5_histogram.txt', header = F, sep = '\t')
d2m <- read.table('SKaTER/Sequencing_Analysis/histograms/MS023rep1_MS023rep2_spawn_expCutoff=5_histogram.txt', header = F, sep = '\t')
dmso <- read.table('SKaTER/Sequencing_Analysis/histograms/DMSOrep1_DMSOrep2_spawn_expCutoff=5_histogram.txt', header = F, sep = '\t')

#compile a list of dataframes
dataframes <- list(d2g = d2g, d2m = d2m, d4c= dmso)
dataframes <- lapply(dataframes, function(x) {x <- x[c(1, 5)]})
colnames <- c('gene', 'rate')
dataframes <- lapply(dataframes, setNames, colnames)

dataframes_merged<-mapply(merge, dataframes, tpm, SIMPLIFY = FALSE)

dataframes_merged <- lapply(dataframes_merged, function(x) {x$thirds <- quantcut(x$rate, q=3, na.rm=TRUE);x})


#combine all the dataframes while keeping the info for original condition analyzed
df <- bind_rows(dataframes_merged, .id = 'df')
df$thirds <- as.character(df$thirds)

df$group <- NA

for (i in 1:nrow(df)) {
  if (df[i,5] == "[0.002,1.03]" | df[i,5] == "[0.0044,1.15]" | df[i,5] == '[0.00207,1.23]') {
    df[i,6] <- 'bottom'
  } 
  else if (df[i,5] == "(1.03,2.71]" | df[i,5] == '(1.15,2.69]' | df[i,5] == "(1.23,3.12]") {
    df[i,6] <- "middle"
  } 
  else {
    df[i,6] <- "top"
  }
}



#df.m <- melt(combined)
#define colors for graph
colors <- c('#D4D4D4','#008837','#7b3294')
#define order of graph
order <- c( 'd4c','d2g', 'd2m')
#convert dataframe names to factors to allow ordering
df$df <- factor(df$df, levels = order)
#uncomment directly below to save pdf
#pdf('R_figures/SKaTER_SpawnRate_vs_TPM_rep_092120.pdf', width = 8, height = 8)
ggplot(df, aes(x=group, y=log(tpm), fill=df)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = colors) +
  #geom_point(aes(fill = df), size =0.5, shape = 21, position = position_jitterdodge()) +
  #geom_jitter(position = position_dodge2(), color="black", size=0.1, alpha=0.1) +
  ggtitle("Spawn Rate vs TPM") +
  xlab("Spawn Rate") + 
  theme(legend.position = "right", legend.key = element_rect(colour = NA, fill = NA)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ylim(0,10)+
  ylab('ln(TPM)')
#dev.off()

cor.test(df[df$df == 'd4c',]$rate, df[df$df == 'd4c',]$tpm, method = 'spearman')
cor.test(df[df$df == 'd2g',]$rate, df[df$df == 'd2g',]$tpm, method = 'spearman')
cor.test(df[df$df == 'd2m',]$rate, df[df$df == 'd2m',]$tpm, method = 'spearman')



#KS test
summary(df[df$df == 'd4c',][,5])
#1st and 3rd tertile; compare tpm distribution via KS test
t.test(df[df$df == 'd4c' & df$group == 'bottom',][,4] , df[df$df == 'd4c' & df$group == 'middle',][,4])
ks.test(df[df$df == 'd4c' & df$thirds == '(2.17e+03,3.37e+03]',][,4] , df[df$df == 'd4c' & df$thirds == '(3.37e+03,1.14e+04]',][,4])
ks.test(df[df$df == 'd4c' & df$thirds == '[149,2.17e+03]',][,4] , df[df$df == 'd4c' & df$thirds == '(2.17e+03,3.37e+03]',][,4])



#ANOVA
a1 <- aov(tpm ~ rate + df, data = df)
summary(a1)
TukeyHSD(a1, which = 'df')
#get aggregate stats
group_by(df, tpm, rate) %>% 
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  )
