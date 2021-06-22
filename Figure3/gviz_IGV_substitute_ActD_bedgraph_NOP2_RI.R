library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
setwd('F:\\ActinomycinD\\bedgraph/')
dmso_minus <- 'A1_Aligned.sortedByCoord.RPKM.bedgraph'
dmso_plus <- 'A4_Aligned.sortedByCoord.RPKM.bedgraph'
d2g_minus <- 'A2_Aligned.sortedByCoord.RPKM.bedgraph'
d2g_plus <- 'A5_Aligned.sortedByCoord.RPKM.bedgraph'
d2m_minus <- 'A3_Aligned.sortedByCoord.RPKM.bedgraph'
d2m_plus <- 'A6_Aligned.sortedByCoord.RPKM.bedgraph'

dTrack1 <- DataTrack(range = dmso_minus, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bg",ylim= c(0,5000))
dTrack2 <- DataTrack(range = dmso_plus, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bg", ylim= c(0,5000))
dTrack3 <- DataTrack(range = d2g_minus, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bg", ylim= c(0,5000))
dTrack4 <- DataTrack(range = d2g_plus, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bg", ylim= c(0,5000))
dTrack5 <- DataTrack(range = d2m_minus, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bg", ylim= c(0,5000))
dTrack6 <- DataTrack(range = d2m_plus, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bg", ylim= c(0,5000))

#bm <- useMart(host = "grch37.ensembl.org", 
#              biomart = "ENSEMBL_MART_ENSEMBL", 
#              dataset = "hsapiens_gene_ensembl")
#bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr7",
#                              start = 39655462, end = 39755127,
#                              name = "ENSEMBL", biomart = bm, filter = list(with_refseq_mrna = FALSE))
txTr <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, chromosome = 'chr12',
                        start = 6669263, end = 6669738)
#plotTracks(dTrack1, chromosome = "1", from = 173833623, to = 173833840)

#plotTracks(c(txTr, dTrack1), chromsome = "12", from = 6669265, to = 6669732)
#plotTracks(bmt)
#plotTracks(txTr)

pdf('NOP2_RI.pdf', height = 20, width = 30)

plotTracks(c(txTr, dTrack1, dTrack2, dTrack3, dTrack4, dTrack5, dTrack6), chromsome = "chr12", from = 6669263, to = 6669738)
dev.off()
