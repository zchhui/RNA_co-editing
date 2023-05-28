rm(list = ls())
options(stringsAsFactors = F)
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
setwd(pat);
exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/exp/LIHC.htseq_fpkm_symbol_90per_unique.txt",row.names=1,header=T,as.is=T);
de <- read.table("deg_up_dn.txt",header=F,as.is=T);###shared dysregulated genes in Figure 7A
sam.gr <- read.table("HCC_sample_group_total.txt",header=T,as.is=T)
high.risk <- sam.gr[which(sam.gr[,2]!="Normal" & sam.gr[,2]!="Tum.S0" ),1]
low.risk <- sam.gr[which(sam.gr[,2]=="Tum.S0" ),1]
gene <- intersect(de[,1],rownames(exp));

#exp
exp.deg <- exp[gene,c(high.risk,low.risk)]
exp.degs <- apply(exp.deg,1,function(x)log2(x+0.05));##log2(FPKM+0.05))
exp.deg <- t(exp.degs);

cli <- read.csv("prog/clinical.csv",header=T,stringsAsFactors=F);
rownames(cli) <- cli[,1];
cli <- cli[-which(cli[,"stage.T"]=="TX"),]
cli[which(cli[,"stage.T"] =="T2a" | cli[,"stage.T"] =="T2b"),"stage.T"] <- "T2";
cli[which(cli[,"stage.T"] =="T3a" | cli[,"stage.T"] =="T3b"),"stage.T"] <- "T3";


data <- data.frame(sample=c(high.risk,low.risk),group=c(rep("high.risk",length(high.risk)),rep("Low.risk",length(low.risk))));
rownames(data) <- c(high.risk,low.risk);
sam <- which(substr(rownames(data),1,12)%in%rownames(cli));
data <- data[sam,]; #TCGA.DD.AACA  duplicated

data$gender <- cli[substring(rownames(data),1,12),"gender"]
data$grade <- cli[substring(rownames(data),1,12),"grade"]
data$stage.T <- cli[substring(rownames(data),1,12),"stage.T"]
data <- data[rownames(data)%in%colnames(exp.deg),]


timer <- read.table("D:/HFU/RNA_editing_combine_regulation/data/download/TIMER-TCGA-immuneEstimation.txt",header=T,as.is=T,sep="\t");
timer[,1] <- gsub("-",".", timer[,1])
colnames(timer)[1] <- "sample"
da.time <- merge(data,timer,by="sample")
write.table(da.time,"annotation_Sample_TIMER.txt",sep="\t",row.names=F,quote=F)
rownames(da.time) <- da.time[,1]

high.risk <- subset(da.time,group=="high.risk",sample)[,1];
low.risk <- subset(da.time,group=="Low.risk",sample)[,1];
da.time <- da.time[,-1]

library(pheatmap);
cols <- colorRampPalette(c("darkred", "red","white", "blue","darkblue"))(100);
ann_colors = list(
    gender = c(MALE="lightblue", FEMALE="pink"),
    group = c(high.risk="#619CFF",Low.risk="#F8766D"),
    stage.T = c(T1 = "grey90", T2 = "lightyellow",T3="yellow",T4="orange"),
    grade = c(G1 = "grey", G2 = "pink",G3="red",G4="darkred")
    )
library(ComplexHeatmap)
pheatmap(exp.deg[,c(high.risk,low.risk)],cluster_rows=T,cluster_cols=F,
     color = cols,annotation_col=da.time,scale="row",show_colnames=F,
     show_rownames=F,annotation_colors = ann_colors)


