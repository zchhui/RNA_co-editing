##plot Fig6B
rm(list=ls());
gc();
pa.res <- "D:/HFU/RNA_editing_combine_regulation/res/";
setwd(pa.res);
exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/exp/LIHC.htseq_fpkm_symbol_90per_unique.txt",row.names=1,header=T,as.is=T);
sam.gr <- read.table("HCC_sample_group_total.txt",header=T,as.is=T)
gene <- c("HLA-A","HLA-B","HLA-C","B2M","CD274");#MHC related genes

exp.g <- exp[gene,];
low.r <- sam.gr[sam.gr[,2]=="Tum.S0",1]
high.r <- sam.gr[(sam.gr[,2]!="Tum.S0" & sam.gr[,2]!="Normal") ,1]
nor <- sam.gr[ sam.gr[,2]=="Normal" ,1]
p <- c()
for(i in 1:5){
 p <- c(p,wilcox.test(as.numeric(exp.g[i,high.r]),as.numeric(exp.g[i,low.r]))$p.value);
 }
 
 
par(mfrow=c(2,3))
pp1 <- c();
pp2 <- c();
pp3 <- c();
for(i in 1:5){
da1 <- list(high.risk=log2(as.numeric(exp.g[i,high.r])+0.05),
            low.risk=log2(as.numeric(exp.g[i,low.r])+0.05),
            Normal=log2(as.numeric(exp.g[i,nor])+0.05)
            )
boxplot(da1,col=c("#619CFF","#F8766D","#00BA38"),main=rownames(exp.g)[i],ylab="Expression level");
points(c(rep(1,length(high.r)),rep(2,length(low.r)),rep(3,length(nor))),c(log2(as.numeric(exp.g[i,high.r])+0.05),log2(as.numeric(exp.g[i,low.r])+0.05),log2(as.numeric(exp.g[i,nor])+0.05)),col=c(rep("#619CFF",length(high.r)),rep("#F8766D",length(low.r)),rep("#00BA38",length(nor)))); 
p1 <- wilcox.test(log2(as.numeric(exp.g[i,high.r])+0.05),log2(as.numeric(exp.g[i,low.r])+0.05))$p.value;
p2 <- wilcox.test(log2(as.numeric(exp.g[i,low.r])+0.05),log2(as.numeric(exp.g[i,nor])+0.05))$p.value;
p3 <- wilcox.test(log2(as.numeric(exp.g[i,high.r])+0.05),log2(as.numeric(exp.g[i,nor])+0.05))$p.value;
pp1 <- c(pp1,round(p1,4));
pp2 <- c(pp2,round(p2,4));
pp3 <- c(pp3,round(p3,4));
 } 

pp1
[1] 0.8563 0.0447 0.0498 0.0153 0.5545
pp2
[1] 0.0000 0.0001 0.0000 0.0616 0.0000
pp3
[1] 0.0000 0.0616 0.0007 0.5713 0.0048

