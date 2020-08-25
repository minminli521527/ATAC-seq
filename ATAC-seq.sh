### 1. Download genome data, build index ###
# 1.1 Download software
#Generally for DNA data, comparison software is generally used: bowtie2, bwa; for RNA data, generally used for comparison: hisat2, star
$ conda create -n atac python=2 bwa
$ conda install bedtools sra-tools deeptools homer bowtie2 trim-galore multiqc sambamba
$ conda install -c bioconda gatk4=4.1.1.0=0
# 1.2 Build an index
$ mkdir index
$ cd index
$ bowtie2-build Arabidopsis_thaliana.TAIR10.dna.toplevel.fa arabiopsis_tair10



### 2.fastac ###
#fastqc
$ mkdir clean align peaks motif qc
$ cd qc
$ find ../data -name "*.gz" |xargs fastqc -t 3 -o ./
#or
$ fastqc -t 5 ../data/*gz -o ./
#multiqc
$ multiqc -n ../qc ./ #qc is the folder name



### 3. Data Filtering###
# 3.1 Trim-galore filters double-ended *gz files
$ trim_galore -q 15 --phred33 --length 26 -e 0.1 --stringency 3 --paired -o clean/ data/SRR5874657.1_1.fastq.gz data/SRR5874657.1_2.fastq.gz
# Trim_galore on double-ended *gz file filtering/batch
************************************************** *********************************
$ cd /clean
$ ls ../data/*_1.fastq.gz> ./1
$ ls ../data/*_2.fastq.gz> ./2
$ paste ./1 ./2> ./config
#New loop script qc.sh
$ cat> ./qc.sh
#Write the loop in the qc.sh script:
source activate atac
dir='/home/minminli/Bioinformatics/atac/clean'
cat /home/minminli/Bioinformatics/atac/clean/config |while read id
do
           arr=($id)
           fq1=${arr[0]}
           fq2=${arr[1]}
trim_galore -q 15 --phred33 --length 26 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
done
source deactivate atac
#Run script
$ bash /home/minminli/Bioinformatics/atac/clean/qc.sh config
************************************************** *********************************

# 3.2 Use fastqc+multiqc software to check data quality
$ cd /clean
$ mkdir fastqc
$ cd fastqc
$ fastqc -t 5 ../*gz -o ./
$ multiqc -n ../fastqc ./



### 4. Comparison ###
# 4.1 Sample processing/single
# (4.1.1) Mapping
$ cd align
#$ bowtie2 -x /home/minminli/Bioinformatics/atac/data/index/arabiopsis_tair10 -1 ../clean/SRR5874657.1_1_val_1.fq.gz -2 ../clean/SRR5874657.1_2_val_2.fq.gz |samtools sort -@ 5 -O bam -o ./SRR5874657.bam-
# Among them, the bowtie2 comparison adds the -X 2000 parameter, which is the largest insert, and a wide insert range (10-1000bp)
$ bowtie2 -p 5 --very-sensitive -X 2000 -x /home/minminli/Bioinformatics/atac/data/index/arabiopsis_tair10 -1 ../clean/SRR5874657.1_1_val_1.fq.gz -2 ../clean/ SRR5874657.1_2_val_2.fq.gz |samtools sort -O bam -@ 5 -o-> SRR5874657.bam
# Count the ratio of mitochondria and chloroplast reads
$ ls *.bam | xargs -i samtools index {}
$ samtools idxstats SRR5874657.bam | cut -f 1,3> SRR5874657.txt #Statistics of Mt/Pt based on the .txt text file
# Count the ratio of mitochondria and chloroplast reads/batch processing
$ ls *.bam | while read id ;do (samtools idxstats $id | cut -f 1,3> ${id%%.*}.txt);done

# 4.2 Sample processing/batch: see tutorial: https://www.jianshu.com/p/09e05bcd6981



### 5. Removal of PCR duplicates + keep multiple aligned reads + high quality + removal of mitochondrial chloroplasts
### (5.1) Remove PCR duplication
# (a) Use GATK picard to remove PCR repeats
$ gatk MarkDuplicates -I SRR5874657.bam -O SRR5874657.picard.rmdup.bam --REMOVE_DUPLICATES true -M test.log
### Use GATK picard to PCR repeat/batch processing
$ ls SRR58746??.bam | while read id ;do (gatk MarkDuplicates -I $id -O ${id%%.*}.picard.rmdup.bam --REMOVE_DUPLICATES true -M ${id%%.*} .log);done
# (b) Or use sambamba to PCR repeat
$ sambamba markdup -r SRR5874657.bam SRR5874657.sambamba.rmdup.bam
# samtools View and compare bam information before and after repeating
$ samtools flagstat SRR5874657.picard.rmdup.bam
$ samtools flagstat SRR5874657.bam

# (5.2) Next, only keep two reads to compare to the same chromosome (Proper paired) (-f 2), and high-quality comparison results (Mapping quality>=30) (-q 30), by the way Filter mitochondrial reads (grep -v chrM)
# View bam
$ samtools view SRR5874657.picard.rmdup.bam |wc
# View bam after -f2 -q30
$ samtools view -f 2 -q 30 SRR5874657.picard.rmdup.bam |wc
# View -f2 -q30 bam after removing mitochondrial Mt
$ samtools view -f 2 -q 30 SRR5874657.picard.rmdup.bam |grep -v Mt|wc
# View -f2 -q30 bam after removing chloroplast Pt
$ samtools view -f 2 -q 30 SRR5874657.picard.rmdup.bam |grep -v Pt|wc
# Output -f2 -q30 bam after removing mitochondrial Mt and chloroplast Pt sort
$ samtools view -h -f 2 -q 30 SRR5874657.picard.rmdup.bam |grep -v Mt|grep -v Pt| samtools sort -O bam -@ 5 -o-> SRR5874657.picard.rmdup.last.bam
# Output -f2 -q30 bam/batch processing after removing mitochondrial Mt and chloroplast Pt sort
$ ls *.picard.rmdup.bam | while read id ;do (samtools view -h -f 2 -q 30 $id |grep -v Mt|grep -v Pt| samtools sort -O bam -@ 5 -o- > ${id%%.*}.picard.rmdup.last.bam );done

# (5.3) Build index for *.picard.rmdup.bam file
$ ls *.picard.rmdup.last.bam | xargs -i samtools index {}
# igvVisualization


### 6. Peak shift ###
#For the reads on the positive strand to be shifted by 4bp to the right, the starting position of the alignment is increased by 4; for the reads on the negative strand, the starting position is shifted to the left by 5bp, and the starting position of the alignment is reduced by 5bp. For the original bam file, use bedtools to convert it into a bed file, and then compare the position offset.
# (6.1) bam to bed-take the intersection of two reads
# $ bedtools bamtobed ​​-bedpe -i SRR5874657.picard.rmdup.last.bam> SRR5874657.picard.rmdup.last.paired.bed
# bam转bed- Intersection of two reads/batch
# $ ls *.picard.rmdup.last.bam | while read id ;do (bedtools bamtobed ​​-bedpe -i $id> ${id%%.*}.picard.rmdup.last.paired.bed);done

# (6.1) bam to bed-single reads
$ bedtools bamtobed ​​-i SRR5874657.picard.rmdup.last.bam> SRR5874657.picard.rmdup.last.single.bed
# bam to bed-single reads/batch
$ ls *.picard.rmdup.last.bam | while read id ;do (bedtools bamtobed ​​-i $id> ${id%%.*}.picard.rmdup.last.single.bed);done

# (6.2) shift
# $ cat SRR5874657.picard.rmdup.last.bed | awk -F'\t''BEGIN {OFS=FS}{if ($6=="+"){$2=$2+4&&$3=$3+4}else if ($6=="-"){$3=$3-5&&$2=$2-5}print $0'}> SRR5874657.picard.rmdup.last.shift.bed
$ cat SRR5874657.picard.rmdup.last.single.bed | awk -F'\t''BEGIN {OFS=FS}{if ($6=="+"){$2=$2+4}else if ($6= ="-"){$3=$3-5}print $0'}> SRR5874657.picard.rmdup.last.single.shift.bed
# shift shift/batch
$ ls *.picard.rmdup.last.single.bed | while read id ;do (cat $id | awk -F'\t''BEGIN {OFS=FS}{if ($6=="+"){$2 =$2+4}else if ($6=="-"){$3=$3-5}print $0'}> ${id%%.*}.picard.rmdup.last.single.shift.bed);done

# (6.3) bed to bam
# Download the genome file Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
$ gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
$ samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
$ cut -f1,2 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai> Arabidopsis_thaliana.TAIR10.genome
$ bedToBam -i SRR5874657.picard.rmdup.last.single.shift.bed -g Arabidopsis_thaliana.TAIR10.genome> SRR5874657.picard.rmdup.last.single.shift.bam
# /Batch
$ ls *.picard.rmdup.last.single.shift.bed | while read id ;do (bedToBam -i $id -g Arabidopsis_thaliana.TAIR10.genome> ${id%%.*}.picard.rmdup.last. single.shift.bam);done

# (6.4) Bam file build index
# Build an index for *.picard.rmdup.last.single.shift.sort.bam files
$ samtools sort SRR5874657.picard.rmdup.last.single.shift.bam -o SRR5874657.picard.rmdup.last.single.shift.sort.bam
$ samtools index SRR5874657.picard.rmdup.last.single.shift.sort.bam
# Build index/batch for *.picard.rmdup.last.single.shift.sort.bam file
************************************************** *********************************
$ ls *.picard.rmdup.last.single.shift.bam | while read id ;do (samtools sort $id -o ${id%%.*}.picard.rmdup.last.single.shift.sort.bam );done
$ ls *.picard.rmdup.last.single.shift.sort.bam | xargs -i samtools index {}
************************************************** *********************************
# igv Visualization



### 7. Calculate the insert length ###
### 7.1 Method1 ###
# (7.1.1) Generate *_length.txt file
# (a) Generate *_length.txt file
$ samtools view SRR5874657.picard.rmdup.last.single.shift.bam |awk'{print $9}'> SRR5874657_length.txt
# (b) Generate *_length.txt file/batch processing
************************************************** *********************************
#Create a config.last_bam file, which contains the name of the bam file
$ ls *.last.single.shift.bam> 1
$ paste 1 1> config.last_bam
#Create the sh file indel_length.sh that extracts the ninth column indel insertion length information of the bam file, and the content is as follows:
cat config.last_bam |while read id;
do
arr=($id)
sample=${arr[0]}
sample_name=${arr[1]}
samtools view $sample |awk'{print $9}'> ${sample_name%%.*}_length.txt
done
#Run indel_length.sh
$ bash ./indel_length.sh
************************************************** *********************************

# (7.1.2) Run on RStudio
************************************************** *********************************
getwd()
setwd("C:/Users/99039/Desktop/atac")
a <-read.table("SRR5874657_length.txt")
a=abs(as.numeric(a[,1]))
pdf('SRR5874657_hist.pdf')
hist(a,
     main="Insertion Size distribution",
     ylab="Read Count",xlab="Insert Size",
     xaxt="n",
     breaks=seq(0,max(a),by=10),col = "lightblue",border="lightblue"
)
axis(side=1,
     at=seq(0,max(a),by=100),
     labels=seq(0,max(a),by=100)
)
dev.off()
************************************************** *********************************

### 7.2 Method2 ###
#Run on RStudio
************************************************** *********************************
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ATACseqQC")
#View working path
getwd()
setwd("C:/Users/99039/Desktop/atac")
## load the library
library(ATACseqQC)
## input the bamFile
bamFiles = "SRR5874657.picard.rmdup.last.single.shift.bam" #Need to have two files bam and bam.bai
bamFiles.labels <- sub(".bam", "", basename(bamFiles))
## generate fragement size distribution
pdf("SRR5874657.picard.rmdup.last.single.shift.pdf")
fragSizeDist(bamFiles, bamFiles.labels)
dev.off()
************************************************** *********************************



### 8. Use macs2 to find peaks ###
### 8.1 Parameters ###
# -t: treat group; -c: control or mock (non-specific antibody, such as IgG) group; -n: name the MACS2 output file; –outdir: MACS2 result file save path; -f: MACS2 read file
#Format,"ELAND", "BED", "ELANDMULTI", "ELANDEXPORT", "ELANDMULTIPET" (for pair-end tags), "SAM", "BAM", "BOWTIE", "BAMPE" or "BEDPE", >The input file format is automatically detected by default, so files in different formats can be used; -g: effective genome size (comparable genome size), the software sets several default values, hs:2.7e9, mm:1.87e9, ce:9e7, dm:1.2e8. -Nomodel: MACS does not build a model.
#ATAC-seq cares about where to cut, the breakpoint is the center of the peak, so use the shift model, --shift -75 or -100
# peak calling
$ macs2 callpeak -t SRR5874657.picard.rmdup.last.bed -g dm --nomodel --shift -100 --extsize 200 -n SRR5874657 --outdir ../peaks/
# peak calling/Batch processing
$ ls *.bed | while read id ;do (macs2 callpeak -t $id -g dm --nomodel --shift -100 --extsize 200 -n ${id%%.*} --outdir ../peaks /) ;done #$(id%%.*) take the prefix
### 8.2 MACS2 output file interpretation ###
# (1) NAME_peaks.xls
#File storing peak information: chromosome name; peak start position; peak end position; peak region length; peak summit position; peak summit position accumulation signal; -log10(pvalue); fold enrichment for this peak summit against random Poisson distribution with local lambda; -log10(qvalue) at peak summit; peak name
# (2) NAME_peaks.narrowPeak
# *narrowPeak files are in BED6+4 format and can be uploaded to UCSC for viewing. The information in columns 1-10 of the output file contains: 1: chromosome number; 2: peak start position; 3: end position; 4: peak name; 5: int(-10*log10qvalue); 6: positive and negative Chain; 7: fold change; 8: -log10pvalue; 9: -log10qvalue; 10: relative summit position to peak start
# (3) NAME_summits.bed
# BED format, including the peak summits (peak highest point) position; if you want to find the motifs of the binding site, it is recommended to use this file.
# (4) NAME_model.r
# R program, after running, generate model pictures based on input data
# $ Rscript NAME_model.r




### 9. Calculate FRiP value ###
# 9.1 Calculate FRiP value/single sample
# Calculate the number of reads in the file in peaks
# -a total files; -b peaks files
$ bedtools intersect -a ../align/SRR5874657.picard.rmdup.last.single.shift.bed -b ./SRR5874657.single_peaks.narrowPeak |wc -l
#16034374
# Calculate the number of reads in the total bed file
$ wc -l ../align/SRR5874657.picard.rmdup.last.single.shift.bed
#28307698
# Therefore the FRiP of 2-cell-1 is 16034374/28307698=0.566431576

# 9.2 Calculate FRiP value/batch
# Copy the original .bed file to this folder
$ cp ../align/*.single.shift.bed ./
# Run
************************************************** *********************************
ls *narrowPeak|while read id;
do
echo $id
# Modify according to the actual situation, $id ".single_peaks.narrowPeak" is only the prefix
bed=$(basename $id ".single_peaks.narrowPeak").picard.rmdup.last.single.shift.bed
Reads=$(bedtools intersect -a $bed -b $id |wc -l|awk'{print $1}')
totalReads=$(wc -l $bed|awk'{print $1}')
echo $Reads $totalReads
echo'==> FRiP value:' $(bc <<< "scale=2;100*$Reads/$totalReads")'%'
done
************************************************** *********************************
# The following results appear:
#16034374 28307698
#==> FRiP value: 56.64%
#SRR5874658.single_peaks.narrowPeak
#21411243 37626506
#==> FRiP value: 56.90%
#SRR5874661.single_peaks.narrowPeak
#11061749 18766996
#==> FRiP value: 58.94%
#SRR5874662.single_peaks.narrowPeak
#7372967 13157406
#==> FRiP value: 56.03%



### 10. Calculate duplication ###
# 10.1 IDR software calculation
# (10.1.1) Calculation
$ conda create -n IDR -y python=3 idr
$ conda activate IDR
$ idr -h
$ idr --samples ../SRR5874657.single_peaks.narrowPeak ../SRR5874658.single_peaks.narrowPeak --plot --log-output-file SRR5874657-58.idr.log
# (10.1.2) Interpretation of results
# (1) idrValues.txt
# In the output file, save all peak results. The format is similar to the input file format, but with a few more columns of information. The first 10 columns are the standard narrowPeak format files, which contain peaks information after repeated sample integration.
$ wc -l idrValues.txt
# score int: such as min(int(log2(-125IDR), 1000), then IDR=0, the scaled IDR is 1000; IDR=0.05, int(-125log2(0.05)) = 540; IDR=1.0, int( -125log2(1.0) = 0.
# If you want to see IDR<0.05, you can filter by the information in the fifth column:
$ awk'{if($5 >= 540) print $0}' idrValues.txt | wc -l
#
# (2) SRR5874657-8.idr.log
The #log file will give the ratio of peaks passing IDR <0.05.
#
# (3) idrValues.txt.png
# The png file includes 4 pictures: upper left: Rep1 peak ranks vs Rep2 peak ranks, peaks that have not passed a specific IDR threshold are displayed in red. Upper right: Rep1 log10 peak scores vs Rep2 log10 peak scores. Peaks that did not pass a specific IDR threshold are displayed in red. The following two figures: Peak rank vs IDR scores. The box plot shows the distribution of IDR values. By default, the IDR value threshold is -1E-6.


# 10.2 Import the narrowPeak file to the local and use local R for visualization
************************************************** *********************************
# Installation package
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c( "ChIPseeker", "ChIPpeakAnno"))
# Load package
library(ChIPseeker)
library(ChIPpeakAnno)
#highPeak <- readPeakFile('SRR5874657.single_peaks.narrowPeak')
#uniquePeak <- readPeakFile('SRR5874658.single_peaks.narrowPeak')
#ol <- findOverlapsOfPeaks(highPeak, uniquePeak)
SRR5874657 <- readPeakFile('SRR5874657.single_peaks.narrowPeak')
SRR5874658 <- readPeakFile('SRR5874658.single_peaks.narrowPeak')
SRR5874661 <- readPeakFile('SRR5874661.single_peaks.narrowPeak')
SRR5874662 <- readPeakFile('SRR5874662.single_peaks.narrowPeak')
ol <- findOverlapsOfPeaks(SRR5874657, SRR5874658,SRR5874661, SRR5874662)
################# Output overlap data ################
# SRR5874657/58/61/62 common overlaps quantity
length(ol$peaklist$"SRR5874657///SRR5874658///SRR5874661///SRR5874662")
# Only SRR5874657 quantity
length(ol$peaklist$"SRR5874657")
# Output SRR5874657/58/61/62 common overlaps data
overlaps <- ol$peaklist$"SRR5874657///SRR5874658///SRR5874661///SRR5874662"
write.csv (overlaps, file ="overlaps-4.csv")
################# Wayne Diagram ################
pdf('overlapVenn-4.pdf')
makeVennDiagram(ol, fill=c("#009E73", "#F0E442", "#D55E00", "#0072B2"), # circle fill color
                col=c("#009E73", "#F0E442", "#D55E00", "#0072B2"), #circle border color
                cat.col=c("#009E73", "#F0E442", "#D55E00", "#0072B2")) # label color, keep same as circle border colors
dev.off()
************************************************** *********************************

# 10.3 Use Bedtools to perform simple overlap and merge duplicate samples
$ bedtools intersect -a ../SRR5874657.single_summits.bed -b ../SRR5874658.single_summits.bed> SRR5874657-58.bed



### 11. deeptools visualization ###
# 11.1 igv visualization
# (11.1.1) bam to bw (bigwig)
# --normalizeUsing RPKM normalization
$ bamCoverage --normalizeUsing RPKM -b SRR5874657.picard.rmdup.last.single.shift.sort.bam -o SRR5874657.picard.rmdup.last.single.shift.bw
# (11.1.1) bam to bw (bigwig)/batch
$ ls *picard.rmdup.last.single.shift.sort.bam | while read id ;do (bamCoverage --normalizeUsing RPKM -b $id -o ${id%%.*}.picard.rmdup.last.single .shift.bw) ;done
#Load the bw/bam/*submits.bed file in igv and observe the peak

# 11.2 Check the signal strength near TSS
# (11.2.1) Build the genome bed
# Install bedops
$ conda activate atac
$ conda install bedops
# First download the corresponding version of Arabidopsis_thaliana.TAIR10.46.gff3 --bed
$ convert2bed --input=gff [--output=bed] <Arabidopsis_thaliana.TAIR10.46.gff3> Arabidopsis_thaliana.TAIR10.46.gff3.bed
# computeMatrix Only the first three columns of the bed file identified by computeMatrix, just take the first three columns of the above bed file
$ cut -f 1-3 Arabidopsis_thaliana.TAIR10.46.gff3.bed> Arabidopsis_thaliana.TAIR10_GFF3.bed

# (11.2.2) computeMatrix command
# reference-point mode
$ computeMatrix reference-point --referencePoint TSS -p 15 -b 2500 -a 2500 -R ./Arabidopsis_thaliana.TAIR10_GFF3.bed -S ./SRR5874657.picard.rmdup.last.single.shift.bw --skipZeros -out SRR5874657_TSS_GFF3 .gz --outFileSortedRegions SRR5874657_TSS_GFF3.bed
# Result visualization
$ plotHeatmap -m SRR5874657_TSS_GFF3.gz -out SRR5874657_TSS_Heatmap.png
# $ plotHeatmap -m SRR5874657_TSS_GFF3.gz -out SRR5874657_TSS_Heatmap.pdf --plotFileFormat pdf --dpi 720
$ plotProfile -m SRR5874657_TSS_GFF3.gz -out SRR5874657_TSS_Profile.png
# $ plotProfile -m SRR5874657_TSS_GFF3.gz -out SRR5874657_TSS_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720
# reference-point mode + result visualization / batch
**********************************tss.sh************* *****************************
bed=/home/minminli/Bioinformatics/atac/peaks/tss/Arabidopsis_thaliana.TAIR10_GFF3.bed
for id in /home/minminli/Bioinformatics/atac/peaks/tss/*.picard.rmdup.last.single.shift.bw
do
echo $id
file=$(basename $id)
sample=${file%%.*}
echo $sample

computeMatrix reference-point --referencePoint TSS -p 15 -b 2500 -a 2500 -R $bed -S $id --skipZeros -out ${sample}_TSS_GFF3.gz --outFileSortedRegions ${sample}_TSS_GFF3.bed

plotHeatmap -m ${sample}_TSS_GFF3.gz -out ${sample}_TSS_Heatmap.png
# plotHeatmap -m ${sample}_TSS_GFF3.gz -out ${sample}_TSS_Heatmap.pdf --plotFileFormat pdf --dpi 720
plotProfile -m ${sample}_TSS_GFF3.g -out ${sample}_TSS_Profile.png
# plotProfile -m ${sample}_TSS_GFF3.gz -out ${sample}_TSS_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720

done
************************************************** *********************************
$ nohup bash tss.sh 1>tss.log &



### 12. peak-Annotation ###
#Use local R for visualization
***********************************************************************************
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
######################## 1. Enter the data ##################### #####
# Get the gff3 file of the reference genome downloaded locally and generate txdb by yourself
# Load gff3 file
GFF_file<-"C:/Users/99039/Desktop/atac/Arabidopsis_thaliana.TAIR10.gff3"
# Generate txdb
txdb <- makeTxDbFromGFF(GFF_file)
# Read in the bed file. This place is the peak file generated by macs2, but it is best to change the name and add the bed suffix. Otherwise, when it is converted into a data frame in R, the header of the first column is misplaced. I don't know why. .
peakSRR5874657 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874657.single_peaks.narrowPeak.bed")
######################## 2.peak Information##################### #####
peak_SRR5874657<-annotatePeak(SRR5874657,tssRegion = c(-2500,2500),TxDb = txdb,addFlankGeneInfo = T,flankDistance = 5000)
# Type the output object peak_SRR5874657 in R, and the Console below will show what region of the ChIPseq site falls on the genome and how it is distributed.
peak_SRR5874657
as.GRanges(peak_SRR5874657)
# Output as.GRanges(peak_SRR5874657), peak information
write.csv (as.GRanges(peak_SRR5874657), file ="peak_SRR5874657.csv")
################## 3. Observe the distribution of peaks in chromosomes #####################
# covplot function can directly accept the output bed file of macs2-callpeak for drawing.
pdf('peakdistribution_SRR5874657.pdf')
covplot(peakSRR5874657, weightCol="V5")
dev.off()
################## 4. Peaks distribution of fixed window ##################
#Select the promoter area as the window, use the transcription start site, then specify the upstream and downstream, use the getPromoters function, by extracting the data in the database, find the upstream and downstream 3000bp range of the promoters area, and prepare the window
promoter <- getPromoters(TxDb=txdb,
                         upstream=3000, downstream=3000)
#Use the getTagMatrix function to compare the peak to this window and generate a matrix for analysis and visualization
tagMatrix <- getTagMatrix(peakSRR5874657,windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000),
           color="red")
#peakHeatmap function Just give the file name, you can directly output the heat map of peaks, which is for the windows upstream and downstream of the promoter area.
pdf('peakHeatmap_SRR5874657.pdf')
peakHeatmap(peakSRR5874657,
            weightCol = NULL, 
            TxDb = txdb,
            upstream = 3000, downstream = 3000, 
            xlab = "", ylab = "", title = NULL, 
            color = "red",
            verbose = TRUE)
dev.off()
################## 5. Show the strength of the combination in the form of a spectrum ##################
#Select the promoter area as the window, use the transcription start site, and then specify the upstream and downstream, use the getPromoters function, by extracting the data in the database, find the upstream and downstream 3000bp range of the promoters area, and prepare the window
promoter <- getPromoters(TxDb=txdb,
                         upstream=3000, downstream=3000)
#Use the getTagMatrix function to compare the peak to this window and generate a matrix for analysis and visualization
tagMatrix <- getTagMatrix(peakSRR5874657,windows=promoter)
pdf('plotAvgProf_SRR5874657.pdf')
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')",
            ylab = "Read Count Frequency")
#Add confidence interval, peak height proves that the combination is tight
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            conf = 0.95, resample = 1000)
dev.off()
################## 6. Visualization of comment information ##################
# Mainly just to note the distribution of the area where the peaks fall, and see where they all fall
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
# 6.1 Pie Chart
pdf('plotAnnoPie_SRR5874657.pdf')
plotAnnoPie(peakAnno)
dev.off()
# 6.2 Histogram
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
upsetplot(peakAnno, vennpie=TRUE)#fusion vennpie
dev.off()
################## 7. Visualizing TF binding loci ##################
pdf('plotDistToTSS_SRR5874657.pdf')
plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()
################## 8. Comparison of multiple peaks ##################
GFF_file<-"C:/Users/99039/Desktop/atac/Arabidopsis_thaliana.TAIR10.gff3"
txdb <- makeTxDbFromGFF(GFF_file)
peakSRR5874657 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874657.single_peaks.narrowPeak.bed")
peakSRR5874658 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874658.single_peaks.narrowPeak.bed")
peakSRR5874661 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874661.single_peaks.narrowPeak.bed")
peakSRR5874662 <- readPeakFile("C:/Users/99039/Desktop/atac/SRR5874662.single_peaks.narrowPeak.bed")
# For multiple peak set comments, first construct a list, and then use lapply.list(name1=bed_file1,name2=bed_file2)?? There is a problem with the data of RYBP, if you add it here, you will always get an error.
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
***********************************************************************************



### 13. motif analysis ###
# (13.1) Install software
#Search homer in bioconda
$ conda activate atac
$ conda search homer -c bioconda
#Install homeer
$ conda install -c bioconda homer
# (13.2) Run
# (13.2.1) Prepare genome file
# the chromosome names followed by white space, i.e. ">chr blahblahblah", not ">chr1-blahblahblah", or prefereably only ">chr1"
# Download the genome file Arabidopsis_thaliana.TAIR10.dna.toplevel.fa, replace> with >chr
$ find -name'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa' | xargs perl -pi -e's|>|>chr|g'
# (13.2.2) Prepare bed file
# Use awk to process the macs2 output file SRR5874657.single_summits.bed, print out the 4, 1, 2, and 3 columns of the original file in turn, each column is separated by the tab key \t, and chr is added before the second column, the processed file Named as SRR5874657.single_summits_homer.bed
# $ awk'{print $4""$1""$2""$3"+"}' SRR5874657.single_summits.bed> SRR5874657.single_summits_homer.bed
$ awk'{print $4"\t""chr"$1"\t"$2"\t"$3"\t+"}' SRR5874657.single_summits.bed> SRR5874657.single_summits_homer.bed
# (13.2.3) Run
# -size # or -size given, default: 200; returning motifs enriched with a p-value less than 0.05; By default, HOMER will search for motifs of len 8, 10, and 12 bp (change using -len <#, #,#> with no spaces between the numbers, ie "-len 6,10,15,20");-p <#>", default 1;-S <#>", default 25;
$ findMotifsGenome.pl -p 15 ../SRR5874657.single_summits_homer.bed ./Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ./ -size 200 -len 8,10,12
# Two peaks comparison: Since HOMER uses a differential motif discovery algorithm, different types of background sequences can be chosen to produce different results. For example, you may want to compare the ChIP-Seq peaks specific in one cell type versus the peaks that are specific to another. To do this, create a second peak/BED file and use it with the argument "-bg <peak/BED file>".
$ findMotifsGenome.pl ../SRR5874657.single_summits_homer.bed ./Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ./ -size 200 -len 8,10,12 -bg ../SRR5874661.single_summits_homer.bed

