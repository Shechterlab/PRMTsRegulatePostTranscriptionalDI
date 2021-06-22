library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\SKaTER\\Sequencing_Analysis\\skater_2\\bigwig\\')
bwFile1 <- 'D1_CKDL200156383-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'
bwFile2<- 'D2_CKDL200156384-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'
bwFile3 <- 'D4_CKDL200156386-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'
bwFile4<- 'D5_CKDL200156387-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'
bwFile5 <- 'D6_CKDL200156388-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'
bwFile6 <- 'D7_CKDL200156389-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'
bwFile7<- 'D8_CKDL200156390-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'
bwFile8 <- 'D9_CKDL200156391-1a_H5MJMDSXY_L4_Aligned.sortedByCoord.PEonly_marked.bw'






dTrack1 <- DataTrack(range = bwFile1, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim= c(0,200))
dTrack2 <- DataTrack(range = bwFile2, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim= c(0,200))
dTrack3 <- DataTrack(range = bwFile3, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim= c(0,200))
dTrack4 <- DataTrack(range = bwFile4, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim= c(0,200))
dTrack5 <- DataTrack(range = bwFile5, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim= c(0,200))
dTrack6 <- DataTrack(range = bwFile6, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim= c(0,200))
dTrack7 <- DataTrack(range = bwFile7, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim= c(0,200))
dTrack8 <- DataTrack(range = bwFile8, genome = "hg19", type = "l", 
                     chromosome = "chr7", name = "bigWig", ylim = c(0,200))
#bm <- useMart(host = "grch37.ensembl.org", 
#              biomart = "ENSEMBL_MART_ENSEMBL", 
#              dataset = "hsapiens_gene_ensembl")
#bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr7",
#                              start = 39655462, end = 39755127,
#                              name = "ENSEMBL", biomart = bm, filter = list(with_refseq_mrna = FALSE))
txTr <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, chromosome = 'chr7',
                        start = 123319997, end = 123391057)

#plotTracks(dTrack2, chromsome = "chr7", from = 179032266, to = 179052362)
#plotTracks(bmt)
#plotTracks(txTr)

pdf('WASL_SKaTER.pdf', height = 10, width = 30)
plotTracks(c(txTr, dTrack1, dTrack2, dTrack3, dTrack4, dTrack5, dTrack6, dTrack7, dTrack8), chromsome = "chr7", from = 123319997, to = 123391057)
dev.off()
