library(seqinr)
library(ggplot2)
library(ggseqlogo)
setwd('C:\\Users\\maxim\\Dropbox (EinsteinMed)\\CloudStation\\PRMT\\RNAseq')
#use bedtools getfasta function to extract nucleotide sequence with strand consideration from hg19 fasta file using bed file with coordinates
downstream_exon <- read.fasta('rMATS_4.0.2/commonRIsequence_downstreamexon.fasta', as.string = TRUE)
intron <- read.fasta('rMATS_4.0.2/commonRIsequence.fasta', as.string = TRUE)


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

intron_end <- substrRight(intron, 30)

downstream <- substr(downstream_exon, 1, 5)

downstream_splice <- as.data.frame(cbind(intron_end, downstream))
downstream_splice$sequence <- paste(downstream_splice$intron_end, downstream_splice$downstream, sep = '')

downstream_splice <- data.frame(lapply(downstream_splice, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

seqs <- list(downstream_splice$sequence)


ggseqlogo(seqs, method = 'prob', seq_type = 'dna')

#ggplot() + 
#  annotate('rect', xmin = 0.5, xmax = 3.5, ymin = -0.05, ymax = 1.9, alpha = .1, col='black', fill='yellow') +
#  geom_logo(seqs, stack_width = 0.90) + 
#  annotate('segment', x = 4, xend=8, y=1.2, yend=1.2, size=2) + 
#  annotate('text', x=6, y=1.3, label='Text annotation') + 
#  theme_logo()
