###Figure 5B PCA
##part 1 group 
rm(list=ls())
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
mt <- read.table("co_edit_sample_matrix_cancer.txt",header=T,sep="\t",stringsAsFactors=F); 
result <- read.table("CO_prog_fdr05_5per_new_adjust_single_ed.txt",sep="\t",header=T,as.is=T);
###fdr<0.05
a <- ncol(result); ##fdr<0.05
ind <- which(result[,a]<0.05);
res <- result[ind,1:2];
colnames(mt)[1:2] <- colnames(res);
mt.r <- merge(mt,res,by=colnames(res));
mt.r1 <- mt.r[,-c(1,2)];
num <- apply(mt.r1,2,sum);
sam.gr <- cbind(names(num),paste0("Tum.S",num))

exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/exp/LIHC.htseq_fpkm_symbol_90per_unique.txt",row.names=1,header=T,as.is=T);
colnames(exp) <- substr(colnames(exp),1,15);
sam <- substr(colnames(exp),14,15);
sam.n <- colnames(exp)[sam=="11"];
sam.n <- cbind(sam.n,"Normal");
colnames(sam.n) <- colnames(sam.gr);
sam.gr <- rbind(sam.gr,sam.n)
write.table(sam.gr,"HCC_sample_group_total.txt",sep="\t",col.names=T,row.names=F,quote=F)


##PCA  Fig5B
sam.gr <- read.table("HCC_sample_group_total.txt",header=T,as.is=T)
library(ggfortify)
###sample classify
anno <- sam.gr
colnames(anno)<-c("sample","group1")
exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/exp/LIHC.htseq_fpkm_symbol_90per_unique.txt",row.names=1,header=T,as.is=T);
colnames(exp) <- substr(colnames(exp),1,15);
exp.T <- apply(exp,1,function(x)log2(x+0.05))
exp.T <- exp.T[anno[,1],]
anno <- data.frame(sam.gr,group="risk")
colnames(anno)<-c("sample","group1","group")
anno[which(anno[,2]!="Normal" & anno[,2]!="Tum.S0" ),3] <- "High-risk";
anno[which(anno[,2]=="Tum.S0"),3] <- "Low-risk";
anno[which(anno[,2]=="Normal"),3] <- "Normal";
autoplot(prcomp(exp.T,scale=T),size=1,data=anno,colour ='group')#
write.csv(anno,"High_low_risk_subgroup.csv",row.names=F,quote=F)

