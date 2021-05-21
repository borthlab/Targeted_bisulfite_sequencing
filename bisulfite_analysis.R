require("ggplot2")
require("reshape2")
require("scatterpie")
require("bsseq")

setwd("results/")
infile<-list.files(recursive = TRUE, pattern="pe.bismark.cov.gz$")
bismarkBSseq <- bsseq::read.bismark(files = infile, 
                                 rmZeroCov = TRUE,
                                 strandCollapse = TRUE,
                  verbose = TRUE)
bismarkBSseq
cells<-c(rep("BFP",3), rep("F8", 3))
type<-c("ori", "neg", "pos", "ori", "neg", "pos")
pData(bismarkBSseq)$cells<-cells
pData(bismarkBSseq)$type<-type
bs<-bismarkBSseq
sampleNames(bs)<-c("BFP.ori","BFP.neg", "BFP.pos", "F8.ori", "F8.neg", "F8.pos")
pData(bs)
bs
chr<-as.character(seqnames(bs))
start<-start(bs)
stop<-end(bs)
cvg<-getCoverage(bs, type="Cov")
colnames(cvg)<-c("BFP.ori","BFP.neg", "BFP.pos", "F8.ori", "F8.neg", "F8.pos")
bs.df<-as.data.frame(cbind(chr,start,stop, cvg))

### Extract methylation levels for CMV
infile<-97036:97038
for(i in infile){
  assign(paste("sample",i, sep="_"), read.table(paste(i, "_1_val_1_bismark_bt2_pe.bismark.cov.gz", sep="")))
}
cmv.df<-as.data.frame(cbind(sample_97036[,1:4],sample_97037[,4],sample_97038[,4]))
colnames(cmv.df)<-c("chr", "start", "stop","BFP.ori","BFP.neg", "BFP.pos")
write.table(cmv.df,file ="../bs.df.cmv.txt",sep="\t",quote=F,row.names=F,col.names=T)

### Extract methylation levels for Fut8
infile<-97039:97041
for(i in infile){
  assign(paste("sample",i, sep="_"), read.table(paste(i, "_1_val_1_bismark_bt2_pe.bismark.cov.gz", sep="")))
}
sample_97039<-sample_97039[which(sample_97039$V5+sample_97039$V6>10),]
sample_97040<-sample_97040[which(sample_97040$V5+sample_97040$V6>10),]
sample_97041<-sample_97041[which(sample_97041$V5+sample_97041$V6>10),]
f8.df<-as.data.frame(cbind(sample_97039[,1:4],sample_97040[,4],sample_97041[,4]))
colnames(f8.df)<-c("chr", "start", "stop","Fut8.ori","Fut8.neg", "Fut8.pos")
write.table(f8.df,file ="../bs.df.fut8.135.txt",sep="\t",quote=F,row.names=F,col.names=T)

### Plot methylation levels for CMV
cmv<-read.table("bs.df.cmv.txt", header=T)
fut8<-read.table("bs.df.fut8.txt", header=T)

cmv.m<-melt(cmv[,c(2,4:6)], id="start", value=c("BFP.ori", "BFP.neg", "BFP.pos"))
cmv.m$unmeth.value<-100-cmv.m$value
cmv.m$id=0
cmv.m$Pos<-seq(from = 1, by = 5, length.out = 14)
cmv.m[which(cmv.m$variable=="BFP.neg"),]$id=3
cmv.m[which(cmv.m$variable=="BFP.pos"),]$id=6
cmv.m[which(cmv.m$variable=="BFP.ori"),]$id=9
colnames(cmv.m)<-c("Start", "Methylation.status", "Methylated","Unmethylated","Sample", "Position")

cmv.p<-ggplot() + 
geom_scatterpie(aes(x=Position, y=Sample, group=Methylation.status), 
				data=cmv.m, 
				cols=c("Methylated", "Unmethylated"), 
				color="#000000", lwd=0.1) + 
	scale_fill_manual(values=c("black", "white" )) +
	scale_x_continuous(breaks = unique(cmv.m$Position), labels = unique(cmv.m$Start)) + 
	scale_y_continuous(breaks = unique(cmv.m$Sample), labels = unique(cmv.m$Methylation.status)) + 
	coord_equal() +
	labs(title="CMV", fill = "State (%)") +
	theme_bw() +
	theme(legend.position = "none")

	
### Plot methylation levels for Fut8
fut8.m<-melt(fut8[,c(2,4:6)], id="start", value=c("Fut8.ori", "Fut8.neg", "Fut8.pos"))
fut8.m$unmeth.value<-100-fut8.m$value
fut8.m$id=0
fut8.m$Pos<-seq(from = 1, by = 5, length.out = 25)
fut8.m[which(fut8.m$variable=="Fut8.neg"),]$id=5
fut8.m[which(fut8.m$variable=="Fut8.pos"),]$id=10
fut8.m[which(fut8.m$variable=="Fut8.ori"),]$id=15
colnames(fut8.m)<-c("Start", "Methylation.status", "Methylated","Unmethylated","Sample", "Position")
fut8.m$Position<-as.numeric(fut8.m$Position)

fut8.p<-ggplot() + 
	geom_scatterpie(aes(x=Position, y=Sample, group=Methylation.status, r=2), 
				data=fut8.m, 
				cols=c("Methylated", "Unmethylated"), 
				color="#000000", lwd=0.1) + 
	scale_x_continuous(breaks = unique(fut8.m$Position), labels = unique(fut8.m$Start)) + 
	scale_y_continuous(breaks = unique(fut8.m$Sample), labels = unique(fut8.m$Methylation.status)) + 
	scale_fill_manual(values=c("black", "white")) +
	coord_equal() +
	labs(title="Fut8", fill = "State (%)") +
	theme_bw() +
	theme(legend.position = "bottom") +
	theme(axis.text.x = element_text(angle = 45)) 

png(file="methylated_pie.png")
plot_grid(cmv.p, fut8.p, nrow = 2, ncol = 1)
dev.off()

pdf(file="methylated_pie.pdf")
plot_grid(cmv.p, fut8.p, nrow = 2, ncol = 1)
dev.off()

