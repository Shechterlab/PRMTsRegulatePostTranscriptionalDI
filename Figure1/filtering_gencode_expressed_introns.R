
#prep TPM
library(tidyverse)
library(gtools)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\')
D2G_1 <- read.table('RNAseq\\Kallisto\\D2G1/abundance.tsv', sep='\t', header =T)
D2G_2 <- read.table('RNAseq\\Kallisto/D2G2/abundance.tsv', sep='\t', header =T)
D2G_3 <- read.table('RNAseq\\Kallisto/D2G3/abundance.tsv', sep='\t', header =T)
D2M_1 <- read.table('RNAseq\\Kallisto\\D2M1/abundance.tsv', sep='\t', header =T)
D2M_2 <- read.table('RNAseq\\Kallisto/D2M2/abundance.tsv', sep='\t', header =T)
D2M_3 <- read.table('RNAseq\\Kallisto/D2M3/abundance.tsv', sep='\t', header =T)
D4C_1 <- read.table('RNAseq\\Kallisto\\D4C1/abundance.tsv', sep='\t', header =T)
D4C_2 <- read.table('RNAseq\\Kallisto/D4C2/abundance.tsv', sep='\t', header =T)
D4C_3 <- read.table('RNAseq\\Kallisto/D4C3/abundance.tsv', sep='\t', header =T)
tpm <- list(D2G_1 = D2G_1, D2G_2 = D2G_2, D2G_3 = D2G_3, D2M_1 = D2M_1, D2M_2 = D2M_2, D2M_3 = D2M_3, D4C_1 = D4C_1, D4C_2 = D4C_2, D4C_3 = D4C_3)
tpm <- tpm %>% reduce(inner_join, by='target_id')
tpm$avgtpm <- rowMeans(tpm[c(grep("^tpm",colnames(tpm)))], na.rm = TRUE)
tpm <- tpm[c(1,38)]
tpm_names <-data.frame(do.call('rbind', strsplit(as.character(tpm$target_id), '|', fixed = TRUE)))
tpm_names <-tpm_names[c(2)]
tpm <- cbind(tpm, tpm_names)
tpm <- tpm[c(3,2)]
colnames(tpm) <- c('gene_id', 'tpm')
tpm <- aggregate(tpm ~ gene_id, data=tpm, FUN=sum)
tpm <- tpm[tpm$tpm > 10,]
tpm$gene_id <- gsub("\\..*","",tpm$gene_id)

##load in all introns gtf
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
hg19 <-rtracklayer::import('RNAseq/gencode.v19.annotation.gtf')
txdb <- makeTxDbFromGRanges(hg19)
all.introns <- intronicParts(txdb)
all.introns <- as.data.frame(all.introns)
#export(all.introns, "gencode.v19.allintrons.gtf")
#introns <- as.data.frame(rtracklayer::import('RNAseq/gencode.v19.allintrons.gtf'))
introns$gene_id <- gsub("\\..*","",introns$gene_id)
##build expressed intron dataframe
expressed_introns <- merge(tpm, introns, by= 'gene_id')
#export(expressed_introns, "gencode.v19.A549expressed.introns.gtf")
