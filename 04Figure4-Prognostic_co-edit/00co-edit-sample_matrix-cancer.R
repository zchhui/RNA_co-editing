###################################
rm(list=ls())
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
coco <- read.table("RADAR_CO_network.txt",as.is=T,header=T,sep="\t"); 
edi <- read.table("D:/HFU/RNA_editing_combine_regulation/data/RADAR_edit.txt",sep="\t",header=T);

can <- which(substr(colnames(edi),14,15)=="01" | substr(colnames(edi),14,15)=="02");
edi.c <- edi[,can];#373

mt <- matrix(nrow=nrow(coco),ncol=ncol(edi.c));

for(i in 1:nrow(coco)){
 for(j in 1:ncol(edi.c)){
mt[i,j] <- as.numeric(!is.na(edi.c[coco[i,1],j]) & !is.na(edi.c[coco[i,2],j]));
 }
}
colnames(mt) <- colnames(edi.c);
mt.r <- cbind(coco[,1:2],mt);
write.table(mt.r,"co_edit_sample_matrix_cancer.txt",sep="\t",row.names=F,quote=F);
