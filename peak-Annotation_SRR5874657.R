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


######################## 1.输入数据 ###########################
# 获取本地下载的参考基因组的gff3文件，自己生成txdb
# 载入gff3文件
GFF_file<-"C:/Users/99039/Desktop/atac/Arabidopsis_thaliana.TAIR10.gff3"
# 生成txdb
txdb <- makeTxDbFromGFF(GFF_file)
# 读入bed文件，这个地方是macs2产生的peak文件，但最好是改名加bed后缀，不然在R中转化为数据框的时候，第一列的抬头有错位，目前还不知道是什么原因。
peakSRR5874657 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874657.single_peaks.narrowPeak.bed")



######################## 2.peak信息 ###########################
peak_SRR5874657<-annotatePeak(SRR5874657,tssRegion = c(-2500,2500),TxDb = txdb,addFlankGeneInfo = T,flankDistance = 5000)
# 在R里打输出的对象peak_SRR5874657，下面Console会显示ChIPseq的位点落在基因组上什么样的区域，分布情况如何。
peak_SRR5874657
as.GRanges(peak_SRR5874657)
# 输出as.GRanges(peak_SRR5874657)，peak信息
write.csv (as.GRanges(peak_SRR5874657), file ="peak_SRR5874657.csv")



################## 3.观察peaks在染色体的分布 #####################
# covplot函数可以直接接受macs2-callpeak的输出bed文件进行画图。
pdf('peakdistribution_SRR5874657.pdf')
covplot(peakSRR5874657, weightCol="V5") 
dev.off()



################## 4.固定窗口的peaks分布 ##################
#选择promoter区作为窗口，使用转录起始位点，然后指定上下游，使用getPromoters函数，通过提取数据库里的数据，找到promoters区域上下游3000bp范围，准备好窗口
promoter <- getPromoters(TxDb=txdb, 
                         upstream=3000, downstream=3000)
#使用getTagMatrix函数，把peak比对到这个窗口，并生成矩阵，供分析、可视化
tagMatrix <- getTagMatrix(peakSRR5874657,windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), 
           color="red")
#peakHeatmap函数只要给文件名，就可以直接出peaks的热图，是针对promoter区上下游的窗口进行绘图
pdf('peakHeatmap_SRR5874657.pdf')
peakHeatmap(peakSRR5874657,
            weightCol = NULL, 
            TxDb = txdb,
            upstream = 3000, downstream = 3000, 
            xlab = "", ylab = "", title = NULL, 
            color = "red",
            verbose = TRUE)
dev.off()



################## 5.谱图的形式来展示结合的强度 ##################
#选择promoter区作为窗口，使用转录起始位点，然后指定上下游，使用getPromoters函数，通过提取数据库里的数据，找到promoters区域上下游3000bp范围，准备好窗口
promoter <- getPromoters(TxDb=txdb, 
                         upstream=3000, downstream=3000)
#使用getTagMatrix函数，把peak比对到这个窗口，并生成矩阵，供分析、可视化
tagMatrix <- getTagMatrix(peakSRR5874657,windows=promoter)
pdf('plotAvgProf_SRR5874657.pdf')
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency")
#添加置信区间,峰高证明结合的紧密
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), 
            conf = 0.95, resample = 1000)
dev.off()



################## 6.注释信息的可视化 ##################
# 主要就是注释一下，peaks峰落在的区域的分布，看看都落在了哪里
peakAnno <-annotatePeak(peakSRR5874657, 
                        tssRegion = c(-3000, 3000), 
                        TxDb = txdb, 
                        level = "transcript", 
                        assignGenomicAnnotation = TRUE, 
                        genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                        annoDb = NULL, 
                        addFlankGeneInfo = FALSE, 
                        flankDistance = 5000, 
                        sameStrand = FALSE, 
                        ignoreOverlap = FALSE, 
                        ignoreUpstream = FALSE, 
                        ignoreDownstream = FALSE, 
                        overlap = "TSS", 
                        verbose = TRUE)
# 6.1 饼图
pdf('plotAnnoPie_SRR5874657.pdf')
plotAnnoPie(peakAnno)
dev.off()
# 6.2 柱状分布图
pdf('plotAnnoBar_SRR5874657.pdf')
plotAnnoBar(peakAnno)
dev.off()
# 6.2 vennpie
pdf('vennpie_SRR5874657.pdf')
vennpie(peakAnno)
dev.off()
# 6.3 UpSetR
pdf('upsetplot_SRR5874657.pdf')
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)#融合vennpie
dev.off()



################## 7.可视化TSS区域的TF binding loci ##################
pdf('plotDistToTSS_SRR5874657.pdf')
plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()




