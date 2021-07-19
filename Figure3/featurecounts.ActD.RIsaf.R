setwd('/home/mmaron2/maxim/featurecounts/ActD/')
library('Rsubread')

args = commandArgs(trailingOnly = TRUE)
if (length(args)<=1) {
  stop("Please provide the metadata (csv), the SAF, and the output file name", call.=FALSE)
}

sampleTable <- read.csv(args[1], row.names=1)
filenames <- paste0(sampleTable$run, '.bam')
saffile <- args[2]
fc <- featureCounts(files = filenames, annot.ext=saffile, isGTFAnnotationFile=FALSE, isPairedEnd=TRUE, nthreads = 35)
colnames(fc$counts) <- sampleTable$run
countdata <- fc$counts
write.csv(countdata, paste0(args[3], '.ActD.countdata.csv'))