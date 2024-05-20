rm(list=ls())
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
de <- read.table("deg_up_dn.txt",header=F,as.is=T);
de <- de[,1];

tar <- read.table("RADAR_targets_fdr05_fc2.txt",header=F,sep="\t",as.is=T);
result <- read.table("CO_prog_fdr05_5per_new_adjust_single_ed.txt",sep="\t",header=T,as.is=T);
a <- ncol(result); ##fdr<0.05
ind <- which(result[,a]<0.05);
result <- result[ind,1:2]; #single site adjust co-editing PAIRS

tar.ed <- list();
com.tar <- c();
for(i in 1:nrow(result)){
 tar1 <- tar[tar[,1]==result[i,1],2];
 tar2 <- tar[tar[,1]==result[i,2],2];
 co.tar <- intersect(tar1,tar2);
tar.ed[[i]] <- co.tar;
com.tar <- c(com.tar,co.tar);
 }
tar3 <- names(table(com.tar))[which( table(com.tar)>=3)];
gen <- intersect(tar3,de)
write.table(gen,"co_tar_deg_gene.txt",sep="\t",row.names=F,quote=F);

library(gplots)
venn(list(co.target=tar3,DEG=de));  #plot Fig7C
