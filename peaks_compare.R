################## 1.多个peak的比较 ##################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
BiocManager::install("ChIPseeker")
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")

getwd()
setwd("C:/Users/99039/Desktop/atac/")
library("GenomicFeatures")
library("ChIPseeker")
library("GenomicRanges")
library("rtracklayer")
library("dplyr")
library("ggplot2")

GFF_file<-"C:/Users/99039/Desktop/atac/Arabidopsis_thaliana.TAIR10.gff3"
txdb <- makeTxDbFromGFF(GFF_file)
peakSRR5874657 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874657.single_peaks.narrowPeak.bed")
peakSRR5874658 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874658.single_peaks.narrowPeak.bed")
peakSRR5874661 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874661.single_peaks.narrowPeak.bed")
peakSRR5874662 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874662.single_peaks.narrowPeak.bed")

# 多个peak set注释时，先构建list,然后用lapply.list(name1=bed_file1,name2=bed_file2)??RYBP的数据有问题，这里加上去，会一直报错。
peaks <- list(peakSRR5874657=peakSRR5874657,peakSRR5874658=peakSRR5874658,peakSRR5874661=peakSRR5874661,peakSRR5874662=peakSRR5874662)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)

pdf('plotAvgProf-1.pdf')
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

pdf('plotAvgProf-2.pdf')
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()

pdf('tagHeatmap.pdf')
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()
