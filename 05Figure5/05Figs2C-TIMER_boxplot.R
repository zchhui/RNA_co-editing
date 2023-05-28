rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/HFU/RNA_editing_combine_regulation/res/")
timer <- read.table("D:/HFU/RNA_editing_combine_regulation/data/download/TIMER-TCGA-immuneEstimation.txt",header=T,as.is=T,sep="\t");
timer[,1] <- gsub("-",".", timer[,1])
colnames(timer)[1] <- "sample"
rownames(timer) <- timer[,1]
timer <- timer[,-1]

ann <- read.table("HCC_sample_group_total.txt",header=T,as.is=T)
timer <- timer[ann[,1],];
timer <- t(timer)

l.r <- ann[ann[,2]=="Tum.S0",1];
h.r <- ann[!(ann[,2]=="Tum.S0" |ann[,2]=="Normal"),1];
nor <- ann[ann[,2]=="Normal",1];

par(mfrow=c(2,3))
pp1 <- c();
pp2 <- c();
pp3 <- c();
for(i in 1:5){
da1 <- list(high.risk=as.numeric(timer[i,h.r]),
            low.risk=as.numeric(timer[i,l.r]),
            Normal=as.numeric(timer[i,nor])
            )
boxplot(da1,col=c("#619CFF","#F8766D","#00BA38"),main=rownames(timer)[i],
            ylab="Inflitration level",ylim=c(0,0.3+round(median(as.numeric(timer[i,h.r])),2)));
library(PMCMR)
shapiro.test(as.numeric(timer[1,]))$p.value #p<0.05, wilcox test

p1 <- wilcox.test(as.numeric(timer[i,h.r]),as.numeric(timer[i,l.r]))$p.value;
p2 <- wilcox.test(as.numeric(timer[i,l.r]),as.numeric(timer[i,nor]))$p.value;
p3 <- wilcox.test(as.numeric(timer[i,h.r]),as.numeric(timer[i,nor]))$p.value;
pp1 <- c(pp1,round(p1,4));
pp2 <- c(pp2,round(p2,4));
pp3 <- c(pp3,round(p3,4));
 }
i=6;
da1 <- list(high.risk=as.numeric(timer[i,h.r]),
            low.risk=as.numeric(timer[i,l.r]),
            Normal=as.numeric(timer[i,nor])
            )
boxplot(da1,col=c("#619CFF","#F8766D","#00BA38"),main=rownames(timer)[i],
            ylab="Inflitration level",ylim=c(0,0.6+round(median(as.numeric(timer[i,h.r])),2)));
p1 <- wilcox.test(as.numeric(timer[i,h.r]),as.numeric(timer[i,l.r]))$p.value;
p2 <- wilcox.test(as.numeric(timer[i,l.r]),as.numeric(timer[i,nor]))$p.value;
p3 <- wilcox.test(as.numeric(timer[i,h.r]),as.numeric(timer[i,nor]))$p.value;
pp1 <- c(pp1,round(p1,4));
pp2 <- c(pp2,round(p2,4));
pp3 <- c(pp3,round(p3,4));
