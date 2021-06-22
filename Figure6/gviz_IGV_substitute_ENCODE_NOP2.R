library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\RNAseq\\ENCODE\\')
bwFile1 <- 'ENCLB555AWP_Aligned.sortedByCoord.PEonly_marked.RPKM.bedgraph'
bwFile2<- 'ENCLB555ASA_Aligned.sortedByCoord.PEonly_marked.RPKM.bedgraph'



dTrack1 <- DataTrack(range = bwFile1, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bigWig", ylim= c(0,15000))
dTrack2 <- DataTrack(range = bwFile2, genome = "hg19", type = "l", 
                     chromosome = "chr12", name = "bigWig", ylim= c(0,15000))

#bm <- useMart(host = "grch37.ensembl.org", 
#              biomart = "ENSEMBL_MART_ENSEMBL", 
#              dataset = "hsapiens_gene_ensembl")
#bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr1",
#                              start = 173832227, end = 173837439,
#                              name = "ENSEMBL", biomart = bm, filter = list(with_refseq_mrna = FALSE))
txTr <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, chromosome = 'chr12', start = 6665835, end = 6677629)

#plotTracks(dTrack2, chromsome = "chr5", from = 179032266, to = 179052362)
#plotTracks(bmt)
#plotTracks(txTr)

pdf('ENCODE_Cyto_Nuc_NOP2.RPKM.pdf', height = 10, width = 30)
plotTracks(c(txTr, dTrack1, dTrack2), chromsome = "chr12", from = 6665835, to = 6677629)
dev.off()
