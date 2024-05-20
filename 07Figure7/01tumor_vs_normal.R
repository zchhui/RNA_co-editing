rm(list=ls())
###Deseq2
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
exp1 <- read.csv("D:/HFU/RNA_editing_combine_regulation/data/exp/TCGA-LIHC_count_90per_symbole.csv",row.names=1,header=T,as.is=T);
sam <- substr(colnames(exp1),14,15);
tum <- colnames(exp1)[sam=="01" | sam=="02"];
nor <- colnames(exp1)[sam=="11"];

library(DESeq2)
colData <- data.frame(sample=c(tum,nor),group=c(rep("Tumor",length(tum)),rep("Normal",length(nor))));
exp <- exp1[,colData[,1]] 

exp <- round(exp)
dds <- DESeqDataSetFromMatrix(countData = exp,colData = colData,design = ~ group)
dds <- DESeq(dds)##normalization

res <- results(dds,contrast=c("group","Tumor","Normal")) 
res <- res[order(res$padj),] # resOrdered=as.data.frame(res)
diff_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))#FDR<0.05 & FC>=2
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
write.csv(res,file= "./deg/TCGA_mRNA_DE_group_tumor_vs_normal.csv",row.names = T)# 得到csv格式的差异表达分析结果

