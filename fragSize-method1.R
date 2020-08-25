if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ATACseqQC")

getwd()
setwd("C:/Users/99039/Desktop/atac")
## load the library
library(ATACseqQC)
## input the bamFile
bamFiles = "GL1.bam"      #���ṩbam��bam.bai�����ļ�
bamFiles.labels <- sub(".bam", "", basename(bamFiles))
## generate fragement size distribution
pdf("ATAC-GL1.pdf")
fragSizeDist(bamFiles, bamFiles.labels)
dev.off()