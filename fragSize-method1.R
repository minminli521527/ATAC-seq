if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ATACseqQC")

getwd()
setwd("C:/Users/99039/Desktop/atac")
## load the library
library(ATACseqQC)
## input the bamFile
bamFiles = "GL1.bam"      #需提供bam和bam.bai两个文件
bamFiles.labels <- sub(".bam", "", basename(bamFiles))
## generate fragement size distribution
pdf("ATAC-GL1.pdf")
fragSizeDist(bamFiles, bamFiles.labels)
dev.off()
