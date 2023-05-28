rm(list=ls())
###Deseq2
exp <- read.csv("D:/HFU/RNA_editing_combine_regulation/data/exp/TCGA-LIHC_count_90per_symbole.csv",row.names=1,header=T,as.is=T);
colnames(exp) <- substr(colnames(exp),1,15);
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
anno <- read.table("HCC_sample_group_total.txt",header=T,as.is=T)
colnames(anno)<-c("sample","group")
anno[which(anno[,2]=="Tum.S1" |anno[,2]=="Tum.S2" |anno[,2]=="Tum.S3" |anno[,2]=="Tum.S4"),2] <- "High.risk";
anno[which(anno[,2]=="Tum.S0"),2] <- "Low.risk";

library(DESeq2)
colData <- anno

exp <- round(exp)
exp <- exp[,anno[,1]]
dds <- DESeqDataSetFromMatrix(countData = exp,colData = anno,design = ~ group)
dds <- DESeq(dds)##normalization

res1 <- results(dds,contrast=c("group","High.risk","Low.risk")) #High risk vs Low.risk
write.csv(res1,"./deg/TCGA_mRNA_DE_group_Syn_vs_Nosyn.csv",row.names = T)

res2 <- results(dds,contrast=c("group","Low.risk","Normal")) # Low.risk vs Normal
write.csv(res2,file= "./deg/TCGA_mRNA_DE_group_Nosyn_vs_normal.csv",row.names = T)

res3 <- results(dds,contrast=c("group","High.risk","Normal")) #High.risk vs Normal
write.csv(res3,file= "./deg/TCGA_mRNA_DE_group_Syner_vs_normal.csv",row.names = T)
