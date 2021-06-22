#the goal of this analysis is to compare the gene length fo RI genes to adjacent genes vs. the A549 expressed hg19 universe
library(tidyverse)
library(reshape2)
library(rtracklayer)
#set your working directory
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\RNAseq')
#load in bed files
hg19 <-  as.data.frame(rtracklayer::import('gencode.v19.annotation.gtf'))
hg19$gene_id <- gsub("\\..*","",hg19$gene_id)
keep <- c('lincRNA', 'snRNA', 'snoRNA', 'protein_coding','miRNA','misc_RNA')
hg19_filtered <- hg19[hg19$gene_type %in% keep,]
expressed.introns <-as.data.frame(rtracklayer::import('gencode.v19.A549expressed.introns.gtf'))
#filter out hg19 genes that are not in A549 D4C RNA-seq
hg19_expressed <- hg19_filtered[(hg19_filtered$gene_id %in% expressed.introns$gene_id),]
hg19_expressed <- hg19_expressed[hg19_expressed$type == 'gene',]
hg19_expressed <- merge(hg19_expressed, expressed.introns, by ='gene_id')

for (i in 1:nrow(hg19_expressed)) {
  if (hg19_expressed[i,grep("strand.x",colnames(hg19_expressed))] == '-') {
    hg19_expressed[i,40] <- hg19_expressed[i,grep("start.y",colnames(hg19_expressed))] - hg19_expressed[i,grep("start.x",colnames(hg19_expressed))]
    } 
  else {
    hg19_expressed[i,40] <- hg19_expressed[i,grep("end.x",colnames(hg19_expressed))] - hg19_expressed[i,grep("end.y",colnames(hg19_expressed))]
  }
}

df <- hg19_expressed[c(1,2,3,4,6,5,14,27,28,29,30,36,38,39,40)]

colnames(df) <- c('gene_id','chr','start','stop','strand','gene_length','gene_name','intron_chr','intron_start','intron_stop','intron_length','tpm','tx_id','tx_name','distance_to_TES')

#write.table(df,'expressed_intron_distance_to_TES.bed',sep='\t', row.names = F, quote = F)
