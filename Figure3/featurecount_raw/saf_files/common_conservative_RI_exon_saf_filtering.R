setwd('C:/Users/maxim/Dropbox (EinsteinMed)/CloudStation/PRMT/RNAseq/rMATS_4.0.2/')
gsk <- read.table('D4C_vs_D2G_RIconservative.saf')
ms023 <- read.table('D4C_vs_D2M_RIconservative.saf')

common <- merge(gsk, ms023, by=c('V2','V3','V4'))
common <- common[c(4,1:3,5,10,7)]

#write.table(common, 'D2G_D2M_commonRI_conservative.saf', col.names = F, row.names = F, quote = F, sep = '\t')


rmats <- read.table('rmats4.0.2_D4C_vs_D2G/RI.MATS.JCEC.txt', header = T, sep = '\t')
rmats$V1.x <- paste(rmats$geneSymbol, rmats$ID, sep='')

data <- merge(common, rmats, by = 'V1.x')

data <- data[c(8:ncol(data))]

#write.table(data, 'D2G_D2M_commonRI.MATS.JCEC.txt', col.names = F, row.names = F, quote = F, sep = '\t')


df <- read.table('D2G_D2M_commonRI_flankingExons.saf')
df.sort <- df[order(df$V3, df$V4),]

df.sort[,8] <- NA
for (i in 2:nrow(df.sort)) {
    if (df.sort[i,3] == df.sort[i-1,3] & df.sort[i,4] > df.sort[i,3]) {
      df.sort[i,8] <- TRUE
    }
}


df.sort <- df.sort[!grepl('TRUE',df.sort$V8),]

df.sort <- df.sort[df.sort$V1 %in% names(which(table(df.sort$V1) > 1)),]

df.sort <- df.sort[c(1:7)]

#write.table(df.sort, 'D2G_D2M_commonRI_flankingExons.filtered.saf', col.names = F, row.names = F, quote = F, sep = '\t')
