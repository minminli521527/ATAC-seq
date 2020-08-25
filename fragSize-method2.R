getwd()
setwd("C:/Users/99039/Desktop/atac")
a <-read.table("SRR5874657_length.txt")    #ͳ��bam�ļ��ھ��У����*_length.txt
a=abs(as.numeric(a[,1]))
pdf('SRR5874657_hist.pdf')
hist(a,
     main="Insertion Size distribution",
     ylab="Read Count",xlab="Insert Size",
     xaxt="n",
     breaks=seq(0,max(a),by=10),col = "lightblue",border="lightblue"
)                         #ȥ��col = "lightblue",border="lightblue"���Ǻ�ɫ�߿�
axis(side=1,
     at=seq(0,max(a),by=100),
     labels=seq(0,max(a),by=100)
)
dev.off() 
