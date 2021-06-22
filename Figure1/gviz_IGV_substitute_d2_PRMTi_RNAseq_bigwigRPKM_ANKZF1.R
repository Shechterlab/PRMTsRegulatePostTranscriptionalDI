library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\RNAseq\\bigwig\\RPKM_bw')
bwFile1 <- 'D4C1_Aligned.sortedByCoord.out.RPKM.bw'
bwFile2<- 'D2G1_Aligned.sortedByCoord.out.RPKM.bw'
bwFile3 <- 'D2M1_Aligned.sortedByCoord.out.RPKM.bw'



dTrack1 <- DataTrack(range = bwFile1, genome = "hg19", type = "l", 
                     chromosome = "chr2", name = "bigWig", ylim= c(0,5000))
dTrack2 <- DataTrack(range = bwFile2, genome = "hg19", type = "l", 
                     chromosome = "chr2", name = "bigWig", ylim= c(0,5000))
dTrack3 <- DataTrack(range = bwFile3, genome = "hg19", type = "l", 
                     chromosome = "chr2", name = "bigWig", ylim= c(0,5000))

#bm <- useMart(host = "grch37.ensembl.org", 
#              biomart = "ENSEMBL_MART_ENSEMBL", 
#              dataset = "hsapiens_gene_ensembl")
#bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr1",
#                              start = 173832227, end = 173837439,
#                              name = "ENSEMBL", biomart = bm, filter = list(with_refseq_mrna = FALSE))
txTr <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, chromosome = 'chr2', start = 220094451, end = 220101505)

#plotTracks(dTrack2, chromsome = "chr5", from = 179032266, to = 179052362)
#plotTracks(bmt)
#plotTracks(txTr)

pdf('2d_PRMTi_ANKZF1.RPKM.pdf', height = 10, width = 30)
plotTracks(c(txTr, dTrack1, dTrack2, dTrack3), chromsome = "chr2", from = 220094451, to = 220101505)
dev.off()