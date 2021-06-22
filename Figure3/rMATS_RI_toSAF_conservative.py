#!/usr/bin/env python
import os
import argparse

#os.chdir('/mnt/c/Users/maxim/Dropbox (EinsteinMed)/CloudStation/PRMT/RNAseq/')

parser=argparse.ArgumentParser()
parser.add_argument("-r","--rMATS",help="name of rMATS file")
#parser.add_argument("-s","--strand", help="which strand to output")
parser.add_argument("-FDR","--FDR",help="Value of FDR cutoff")
parser.add_argument("-f","--fileName",help="Name of output file .gtf")
args = parser.parse_args()

#rmats = '/mnt/c/Users/maxim/Dropbox (EinsteinMed)/CloudStation/PRMT/RNAseq/rMATS_4.0.2/rmats4.0.2_D4C_vs_D2G/RI.MATS.JCEC.txt'
rmats = args.rMATS
genes = []
outF = open(args.fileName, 'w')
with open(rmats, 'r') as f:
    next(f)
    for line in f:
        ID, GeneID, geneSymbol, chr, strand, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, ID, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, SkipFormLen, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference = line.strip().split('\t')
        if float(FDR) < float(args.FDR):
            upstreamEE = int(upstreamEE) + 10
            downstreamES = int(downstreamES) - 10
            outF.write('\t'.join([geneSymbol.strip('"')+ID, chr, str(upstreamEE), str(downstreamES), strand, geneSymbol.strip('"'), GeneID.strip('"')])+'\n')