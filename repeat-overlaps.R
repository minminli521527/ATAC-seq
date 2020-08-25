options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c( "ChIPseeker", "ChIPpeakAnno"))




################# Method /�������� #################
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
################# �鿴/����overlaps #################
# SRR5874657��SRR5874658���е�peaks����
length(ol$peaklist$"SRR5874657///SRR5874658")
# ��SRR5874657�е�peaks����
length(ol$peaklist$"SRR5874657")
# ��SRR5874658�е�peaks����
length(ol$peaklist$"SRR5874658")
# ����ص���peaks��Ϣ
overlaps <- ol$peaklist$"SRR5874657///SRR5874658"
write.csv (overlaps, file ="overlaps-2.csv")
################# Τ��ͼ #################
pdf('overlapVenn-2.pdf')
makeVennDiagram(ol, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2")) # label color, keep same as circle border colors
dev.off()




########### Method /������� ###########
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
################# �鿴/����overlaps #################
# SRR5874657��SRR5874658��SRR5874661��SRR5874662���е�peaks����
length(ol$peaklist$"SRR5874657///SRR5874658///SRR5874661///SRR5874662")
# ��SRR5874657�е�peaks����
length(ol$peaklist$"SRR5874657")
# ����ص���peaks��Ϣ
overlaps <- ol$peaklist$"SRR5874657///SRR5874658///SRR5874661///SRR5874662"
write.csv (overlaps, file ="overlaps-4.csv")
################# Τ��ͼ #################
pdf('overlapVenn-4.pdf')
makeVennDiagram(ol, fill=c("#009E73", "#F0E442", "#D55E00", "#0072B2"), # circle fill color
                col=c("#009E73", "#F0E442", "#D55E00", "#0072B2"), #circle border color
                cat.col=c("#009E73", "#F0E442", "#D55E00", "#0072B2")) # label color, keep same as circle border colors
dev.off()


