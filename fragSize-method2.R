getwd()
setwd("C:/Users/99039/Desktop/atac")
a <-read.table("SRR5874657_length.txt")    #统计bam文件第九列，输出*_length.txt
a=abs(as.numeric(a[,1]))
pdf('SRR5874657_hist.pdf')
hist(a,
     main="Insertion Size distribution",
     ylab="Read Count",xlab="Insert Size",
     xaxt="n",
     breaks=seq(0,max(a),by=10),col = "lightblue",border="lightblue"
)                         #去掉col = "lightblue",border="lightblue"就是黑色边框
axis(side=1,
     at=seq(0,max(a),by=100),
     labels=seq(0,max(a),by=100)
)
dev.off() 

