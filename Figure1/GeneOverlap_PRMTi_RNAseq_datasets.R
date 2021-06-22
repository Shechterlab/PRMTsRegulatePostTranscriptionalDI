library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(GeneOverlap)

setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\RNAseq')
d2g_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D2G/RI.MATS.JCEC.txt', header = T, sep = '\t')
d2gm_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D2GM/RI.MATS.JCEC.txt', header = T, sep = '\t')
d2m_RI <- read.table('rMATS_4.0.2/rmats4.0.2_D4C_vs_D2M/RI.MATS.JCEC.txt', header = T, sep = '\t')

d2g_RI <- d2g_RI[d2g_RI$FDR < 0.05,]
d2m_RI <- d2m_RI[d2m_RI$FDR < 0.05,]
d2gm_RI <- d2gm_RI[d2gm_RI$FDR < 0.05,]


#d2g_RI <- d2g_RI[d2g_RI$IncLevelDifference < 0,]
#d2m_RI <- d2m_RI[d2m_RI$IncLevelDifference > 0,]
#d2gm_RI <- d2gm_RI[d2gm_RI$IncLevelDifference < 0,]


d2g_RI <- d2g_RI[!duplicated(c(d2g_RI$upstreamEE, d2g_RI$downstreamES)),]
d2m_RI <- d2m_RI[!duplicated(c(d2m_RI$upstreamEE, d2m_RI$downstreamES)),]
d2gm_RI <- d2gm_RI[!duplicated(c(d2gm_RI$upstreamEE, d2gm_RI$downstreamES)),]

d2g_RI <- d2g_RI[!is.na(d2g_RI$riExonStart_0base),]
d2m_RI <- d2m_RI[!is.na(d2m_RI$riExonStart_0base),]
d2gm_RI <- d2gm_RI[!is.na(d2gm_RI$riExonStart_0base),]


d2g_RI$region <- paste(d2g_RI$upstreamEE, d2g_RI$downstreamES, sep = '')
d2m_RI$region <- paste(d2m_RI$upstreamEE, d2m_RI$downstreamES, sep = '')
d2gm_RI$region <- paste(d2gm_RI$upstreamEE, d2gm_RI$downstreamES, sep = '')


BoutzCancerCell <- read.table('rMATS_4.0.2/PRJNA362466_BoutzCancerCell_EPZ015666/RI.MATS.JCEC.txt', header = T, sep = '\t')
BoutzCancerCell <- BoutzCancerCell[BoutzCancerCell$FDR < 0.05,]
#BoutzCancerCell <- BoutzCancerCell[BoutzCancerCell$IncLevelDifference < 0,]
BoutzCancerCell <- BoutzCancerCell[!duplicated(c(BoutzCancerCell$upstreamEE, BoutzCancerCell$downstreamES)),]
BoutzCancerCell <- BoutzCancerCell[!is.na(BoutzCancerCell$riExonStart_0base),]
BoutzCancerCell$region <- paste(BoutzCancerCell$upstreamEE, BoutzCancerCell$downstreamES, sep = '')

FedoriwCancerCell_GSK591 <- read.table('rMATS_4.0.2/PRJNA522801_Fedoriw_CancerCell_MTAPwt_GSK591/RI.MATS.JCEC.txt', header = T, sep = '\t')
FedoriwCancerCell_GSK591 <- FedoriwCancerCell_GSK591[FedoriwCancerCell_GSK591$FDR < 0.05,]
#FedoriwCancerCell_GSK591 <- FedoriwCancerCell_GSK591[FedoriwCancerCell_GSK591$IncLevelDifference < 0,]
FedoriwCancerCell_GSK591 <- FedoriwCancerCell_GSK591[!duplicated(c(FedoriwCancerCell_GSK591$upstreamEE, FedoriwCancerCell_GSK591$downstreamES)),]
FedoriwCancerCell_GSK591 <- FedoriwCancerCell_GSK591[!is.na(FedoriwCancerCell_GSK591$riExonStart_0base),]
FedoriwCancerCell_GSK591$region <- paste(FedoriwCancerCell_GSK591$upstreamEE, FedoriwCancerCell_GSK591$downstreamES, sep = '')

FedoriwCancerCell_GSK712 <- read.table('rMATS_4.0.2/PRJNA522801_Fedoriw_CancerCell_MTAPwt_GSK712/RI.MATS.JCEC.txt', header = T, sep = '\t')
FedoriwCancerCell_GSK712 <- FedoriwCancerCell_GSK712[FedoriwCancerCell_GSK712$FDR < 0.05,]
#FedoriwCancerCell_GSK712 <- FedoriwCancerCell_GSK712[FedoriwCancerCell_GSK712$IncLevelDifference < 0,]
FedoriwCancerCell_GSK712 <- FedoriwCancerCell_GSK712[!duplicated(c(FedoriwCancerCell_GSK712$upstreamEE, FedoriwCancerCell_GSK712$downstreamES)),]
FedoriwCancerCell_GSK712 <- FedoriwCancerCell_GSK712[!is.na(FedoriwCancerCell_GSK712$riExonStart_0base),]
FedoriwCancerCell_GSK712$region <- paste(FedoriwCancerCell_GSK712$upstreamEE, FedoriwCancerCell_GSK712$downstreamES, sep = '')

FedoriwCancerCell_combination <- read.table('rMATS_4.0.2/PRJNA522801_Fedoriw_CancerCell_MTAPwt_combination/RI.MATS.JCEC.txt', header = T, sep = '\t')
FedoriwCancerCell_combination <- FedoriwCancerCell_combination[FedoriwCancerCell_combination$FDR < 0.05,]
#FedoriwCancerCell_combination <- FedoriwCancerCell_combination[FedoriwCancerCell_combination$IncLevelDifference < 0,]
FedoriwCancerCell_combination <- FedoriwCancerCell_combination[!duplicated(c(FedoriwCancerCell_combination$upstreamEE, FedoriwCancerCell_combination$downstreamES)),]
FedoriwCancerCell_combination <- FedoriwCancerCell_combination[!is.na(FedoriwCancerCell_combination$riExonStart_0base),]
FedoriwCancerCell_combination$region <- paste(FedoriwCancerCell_combination$upstreamEE, FedoriwCancerCell_combination$downstreamES, sep = '')

Guccione_GSK591 <- read.table('rMATS_4.0.2/PRJNA509783_GuccioneCancerCell_K562wt_GSK591/RI.MATS.JCEC.txt', header = T, sep = '\t')
Guccione_GSK591 <- Guccione_GSK591[Guccione_GSK591$FDR < 0.05,]
#Guccione_GSK591 <- Guccione_GSK591[Guccione_GSK591$IncLevelDifference < 0,]
Guccione_GSK591 <- Guccione_GSK591[!duplicated(c(Guccione_GSK591$upstreamEE, Guccione_GSK591$downstreamES)),]
Guccione_GSK591 <- Guccione_GSK591[!is.na(Guccione_GSK591$riExonStart_0base),]
Guccione_GSK591$region <- paste(Guccione_GSK591$upstreamEE, Guccione_GSK591$downstreamES, sep = '')

Guccione_cotreat <- read.table('rMATS_4.0.2/PRJNA509783_GuccioneCancerCell_K562wt_GSK591+MS023/RI.MATS.JCEC.txt', header = T, sep = '\t')
Guccione_cotreat <- Guccione_cotreat[Guccione_cotreat$FDR < 0.05,]
#Guccione_cotreat <- Guccione_cotreat[Guccione_cotreat$IncLevelDifference < 0,]
Guccione_cotreat <- Guccione_cotreat[!duplicated(c(Guccione_cotreat$upstreamEE, Guccione_cotreat$downstreamES)),]
Guccione_cotreat <- Guccione_cotreat[!is.na(Guccione_cotreat$riExonStart_0base),]
Guccione_cotreat$region <- paste(Guccione_cotreat$upstreamEE, Guccione_cotreat$downstreamES, sep = '')

Guccione_MS023 <- read.table('rMATS_4.0.2/PRJNA509783_GuccioneCancerCell_K562wt_MS023/RI.MATS.JCEC.txt', header = T, sep = '\t')
Guccione_MS023 <- Guccione_MS023[Guccione_MS023$FDR < 0.05,]
#Guccione_MS023 <- Guccione_MS023[Guccione_MS023$IncLevelDifference < 0,]
Guccione_MS023 <- Guccione_MS023[!duplicated(c(Guccione_MS023$upstreamEE, Guccione_MS023$downstreamES)),]
Guccione_MS023 <- Guccione_MS023[!is.na(Guccione_MS023$riExonStart_0base),]
Guccione_MS023$region <- paste(Guccione_MS023$upstreamEE, Guccione_MS023$downstreamES, sep = '')

Helin_sgPRMT5 <- read.table('rMATS_4.0.2/rMATS_4.0.2_sgNeg_vs_sgPRMT5_gencodev19/RI.MATS.JCEC.txt', header = T, sep = '\t')
Helin_sgPRMT5 <- Helin_sgPRMT5[Helin_sgPRMT5$FDR < 0.05,]
#Helin_sgPRMT5 <- Helin_sgPRMT5[Helin_sgPRMT5$IncLevelDifference < 0,]
Helin_sgPRMT5 <- Helin_sgPRMT5[!duplicated(c(Helin_sgPRMT5$upstreamEE, Helin_sgPRMT5$downstreamES)),]
Helin_sgPRMT5 <- Helin_sgPRMT5[!is.na(Helin_sgPRMT5$riExonStart_0base),]
Helin_sgPRMT5$region <- paste(Helin_sgPRMT5$upstreamEE, Helin_sgPRMT5$downstreamES, sep = '')


library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\RNAseq')
hg19 <-rtracklayer::import('gencode.v19.annotation.gtf')
hg19 <- hg19[!hg19@elementMetadata@listData$gene_type == 'pseudogene',]
hg19 <- hg19[hg19@elementMetadata@listData$gene_status == 'KNOWN',]
hg19 <- hg19[!hg19@seqnames@values == 'chrM',]
txdb <- makeTxDbFromGRanges(hg19)
all.introns <- intronicParts(txdb)
all.introns <- as.data.frame(all.introns)

#gsk591 <- list(A549_GSK591 = d2g_RI$region, U87_EPZ015666 = BoutzCancerCell$region, K562_GSK591 = Guccione_GSK591$region, THP1_sgPRMT5 = Helin_sgPRMT5$region)
#gom.self_GSK591 <- newGOM(gsk591, genome.size=nrow(all.introns))

#pdf('../R_figures/oddsratio_overlap_RI_from_multipleRNAseq_d2gGSK_odds.pdf', height=8, width=8)
#drawHeatmap(gom.self_GSK591, what='Jaccard', grid.col='Greens')
#dev.off()
#pdf('../R_figures/oddsratio_overlap_RI_from_multipleRNAseq_d2gGSK_JaccardIndex.pdf', height=8, width=8)
#drawHeatmap(gom.self_GSK591, what='Jaccard', grid.col='Greens')
#dev.off()

#ms023 <- list(A549_MS023 = d2m_RI$region, K562_MS023 = Guccione_MS023$region)
#gom.self_MS023 <- newGOM(ms023, genome.size=nrow(all.introns))

#pdf('../R_figures/oddsratio_overlap_RI_from_multipleRNAseq_d2gMS023.pdf', height=8, width=8)
#drawHeatmap(gom.self_MS023)
#dev.off()

#cotreat <- list(A549_GSK591_MS023 = d2gm_RI$region, K562_GSK591_MS023 = Guccione_cotreat$region)
#gom.self_cotreat <- newGOM(cotreat, genome.size=nrow(all.introns))

#pdf('../R_figures/oddsratio_overlap_RI_from_multipleRNAseq_d2gMS023.pdf', height=8, width=8)
#drawHeatmap(gom.self_cotreat)
#dev.off()

#gom.obj <- newGOM(list(d2m_RI$region), list(Guccione_MS023$region), genome.size = 364334)
#drawHeatmap(gom.obj)

global <- list(A549_GSK591 = d2g_RI$region, U87_EPZ015666 = BoutzCancerCell$region, K562_GSK591 = Guccione_GSK591$region, THP1_sgPRMT5 = Helin_sgPRMT5$region, A549_MS023 = d2m_RI$region, K562_MS023 = Guccione_MS023$region, A549_GSK591_MS023 = d2gm_RI$region, K562_GSK591_MS023 = Guccione_cotreat$region, PANC03.27_GSK591 = FedoriwCancerCell_GSK591$region, PANC03.27_GSK712 = FedoriwCancerCell_GSK712$region, PANC03.27_GSK591_GSK712 = FedoriwCancerCell_combination$region)

#273902 comes from intersecting intergenic and exonic portions of gtf file used for rMATS; alternative is to use nrow(all.introns)
gom.self_global <- newGOM(global, genome.size=273902)

pdf('../R_figures/102320_RI_overlap_multipleRNAseq_globalRI_log2odds_adjp_oranges.pdf', height= 16, width= 16)
drawHeatmap(gom.self_global, what='odds.ratio', grid.col='Oranges', log.scale=T, adj.p = T, note.col='black', , ncolused = 9)
dev.off()

pdf('../R_figures/102320_RI_overlap_multipleRNAseq_globalRI_JaccardIndex_adjp_oranges.pdf', height= 16, width= 16)
drawHeatmap(gom.self_global, what='Jaccard', grid.col='Oranges', log.scale=F, adj.p = T, note.col='black', ncolused = 9)
dev.off()


