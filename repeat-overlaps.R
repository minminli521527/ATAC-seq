options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c( "ChIPseeker", "ChIPpeakAnno"))




################# Method /两个样本 #################
library(ChIPseeker)
library(ChIPpeakAnno)
getwd()
setwd("C:/Users/99039/Desktop/atac")

#highPeak <- readPeakFile( 'SRR5874657.single_peaks.narrowPeak')
#uniquePeak <- readPeakFile( 'SRR5874658.single_peaks.narrowPeak')
#ol <- findOverlapsOfPeaks(highPeak, uniquePeak)
SRR5874657 <- readPeakFile( 'SRR5874657.single_peaks.narrowPeak')
SRR5874658 <- readPeakFile( 'SRR5874658.single_peaks.narrowPeak')
ol <- findOverlapsOfPeaks(SRR5874657, SRR5874658)
################# 查看/导出overlaps #################
# SRR5874657和SRR5874658都有的peaks数量
length(ol$peaklist$"SRR5874657///SRR5874658")
# 仅SRR5874657有的peaks数量
length(ol$peaklist$"SRR5874657")
# 仅SRR5874658有的peaks数量
length(ol$peaklist$"SRR5874658")
# 输出重叠的peaks信息
overlaps <- ol$peaklist$"SRR5874657///SRR5874658"
write.csv (overlaps, file ="overlaps-2.csv")
################# 韦恩图 #################
pdf('overlapVenn-2.pdf')
makeVennDiagram(ol, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2")) # label color, keep same as circle border colors
dev.off()




########### Method /多个样本 ###########
library(ChIPseeker)
library(ChIPpeakAnno)
#highPeak <- readPeakFile( 'SRR5874657.single_peaks.narrowPeak')
#uniquePeak <- readPeakFile( 'SRR5874658.single_peaks.narrowPeak')
#ol <- findOverlapsOfPeaks(highPeak, uniquePeak)
SRR5874657 <- readPeakFile( 'SRR5874657.single_peaks.narrowPeak')
SRR5874658 <- readPeakFile( 'SRR5874658.single_peaks.narrowPeak')
SRR5874661 <- readPeakFile( 'SRR5874661.single_peaks.narrowPeak')
SRR5874662 <- readPeakFile( 'SRR5874662.single_peaks.narrowPeak')
ol <- findOverlapsOfPeaks(SRR5874657, SRR5874658,SRR5874661, SRR5874662)
################# 查看/导出overlaps #################
# SRR5874657和SRR5874658、SRR5874661、SRR5874662都有的peaks数量
length(ol$peaklist$"SRR5874657///SRR5874658///SRR5874661///SRR5874662")
# 仅SRR5874657有的peaks数量
length(ol$peaklist$"SRR5874657")
# 输出重叠的peaks信息
overlaps <- ol$peaklist$"SRR5874657///SRR5874658///SRR5874661///SRR5874662"
write.csv (overlaps, file ="overlaps-4.csv")
################# 韦恩图 #################
pdf('overlapVenn-4.pdf')
makeVennDiagram(ol, fill=c("#009E73", "#F0E442", "#D55E00", "#0072B2"), # circle fill color
                col=c("#009E73", "#F0E442", "#D55E00", "#0072B2"), #circle border color
                cat.col=c("#009E73", "#F0E442", "#D55E00", "#0072B2")) # label color, keep same as circle border colors
dev.off()



