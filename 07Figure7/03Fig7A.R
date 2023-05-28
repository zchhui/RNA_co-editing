rm(list=ls())
pat <- "D:/HFU/RNA_editing_combine_regulation/res/deg/";
options(stringsAsFactors = F)
setwd(pat);
da1 <- read.csv("TCGA_mRNA_DE_group_Syn_vs_Nosyn.csv",header = T,sep=",",as.is=T) 
da2 <- read.csv("TCGA_mRNA_DE_group_Nosyn_vs_normal.csv",header = T,sep=",",as.is=T) 
da3 <- read.csv("TCGA_mRNA_DE_group_Syner_vs_normal.csv",header = T,sep=",",as.is=T) 
da4 <- read.csv("TCGA_mRNA_DE_group_tumor_vs_normal.csv",header=T,as.is=T,sep=",")
Syn.Nosyn.up <- subset(da1,padj < 0.05 & log2FoldChange > 1 )### FDR<0.05,FC>2
Syn.Nosyn.up <- Syn.Nosyn.up[,1]
Syn.Nosyn.dn <- subset(da1,padj < 0.05 &  log2FoldChange < -1)
Syn.Nosyn.dn <- Syn.Nosyn.dn[,1]

NoSyn.normal.up <- subset(da2,padj < 0.05 & log2FoldChange > 1 )### FDR<0.05,FC>2
NoSyn.normal.up <- NoSyn.normal.up[,1]
NoSyn.normal.dn <- subset(da2,padj < 0.05 &  log2FoldChange < -1)
NoSyn.normal.dn <- NoSyn.normal.dn[,1]

Syn.normal.up <- subset(da3,padj < 0.05 & log2FoldChange > 1 )### FDR<0.05,FC>2
Syn.normal.up <- Syn.normal.up[,1]
Syn.normal.dn <- subset(da3,padj < 0.05 &  log2FoldChange < -1)
Syn.normal.dn <- Syn.normal.dn[,1]

Tum.normal.up <- subset(da4,padj < 0.05 & log2FoldChange > 1 )### FDR<0.05,FC>2
Tum.normal.up <- Tum.normal.up[,1]
Tum.normal.dn <- subset(da4,padj < 0.05 &  log2FoldChange < -1)
Tum.normal.dn <- Tum.normal.dn[,1]


library(gplots)
par(mfrow=c(1,2))
data1 <- list(Tum2Nor.up=Tum.normal.up,Syn2Nor.up=Syn.normal.up,Nosyn2Nor=NoSyn.normal.up,Syn2Nosyn=Syn.Nosyn.up);
venn(data1)
data2 <- list(Tum2Nor.dn=Tum.normal.dn,Syn2Nor.dn=Syn.normal.dn,Nosyn2Nor=NoSyn.normal.dn,Syn2Nosyn=Syn.Nosyn.dn);
venn(data2)
com.up <- intersect(intersect(intersect(Tum.normal.up,Syn.normal.up),NoSyn.normal.up),Syn.Nosyn.up)
com.dn <- intersect(intersect(intersect(Tum.normal.dn,Syn.normal.dn),NoSyn.normal.dn),Syn.Nosyn.dn)
up.dn <- rbind(cbind(com.up,"up"),cbind(com.dn,"dn"));
write.table(up.dn,"../deg_up_dn.txt",sep="\t",col.names=F,row.names=F,quote=F)

