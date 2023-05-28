rm(list=ls())
library(survival);

pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
cli <- read.csv("./prog/clinical.csv",header=T,stringsAsFactors=F);
rownames(cli) <- cli[,1];
mt <- read.table("co_edit_sample_matrix_cancer.txt",header=T,sep="\t",stringsAsFactors=F); 
mtr <- mt[,substr(colnames(mt),1,12)%in%rownames(cli)];
mt <- cbind(mt[,1:2],mtr);
num <- apply(mtr,1,function(x){sum(x==1)});
cuf <- round(0.05*ncol(mtr));##at least 5% co-edit in HCC 
index <- which(num>=cuf);
mtr <- mtr[index,];
mt <- mt[index,];
cli <- cli[intersect(substr(colnames(mt),1,12),rownames(cli)),];

result <- read.table("CO_prog_fdr05_5per_new_adjust_single_ed.txt",sep="\t",header=T,as.is=T);
###fdr<0.05
a <- ncol(result); ##fdr<0.05
ind <- which(result[,a]<0.05);
res <- result[ind,1:2];
colnames(mt)[1:2] <- colnames(res);
mt.r <- merge(mt,res,by=colnames(res));
mt.r1 <- mt.r[,-c(1,2)];
num <- apply(mt.r1,2,sum);
 
p.val <- c();
par(mfrow=c(1,1));
clic <- cbind(cli,cli[,1]);
clic[substr(names(num)[num==0],1,12),ncol(clic)] <- 0;  ##co.ed=0
clic[substr(names(num)[num==1],1,12),ncol(clic)] <- 1;  ###co.ed=1
clic[substr(names(num)[num==2],1,12),ncol(clic)] <- 2;  ###co.ed=2
clic[substr(names(num)[num==3],1,12),ncol(clic)] <- 3;  ###co.ed=3
clic[substr(names(num)[num==4],1,12),ncol(clic)] <- 4;  ###co.ed=4
colnames(clic)[ncol(clic)] <- "co_ed";

prog.dif <- survdiff(Surv(time,status)~co_ed,data=clic);
prog.fit <- survfit(Surv(time,status)~co_ed,data=clic);
pp <- prog.dif$pvalue
colorss <- c("#7F847F","#CCCCFF","#9999FF","#9966FF","#6633FF");
plot(prog.fit,col=colorss,lwd=2,xlab="Survival time(days)",ylab="Survival(%)")
legend(2000,1,c("No_Co-Edit","Risk Co-Edit=1","Risk Co-Edit=2","Risk Co-Edit=3","Risk Co-Edit=4"),col=colorss,lty=1,text.col=colorss,lwd=2);
text(2500,0.31,paste0("No_Edit: n=",length(which(num==0))),col="#7F847F");
text(2500,0.26,paste0("1Risk Co-Edit: n=",length(which(num==1))),col="#CCCCFF");
text(2500,0.21,paste0("2Risk Co-Edit: n=",length(which(num==2))),col="#9999FF");
text(2500,0.16,paste0("3Risk Co-Edit: n=",length(which(num==3))),col="#9966FF");
text(2500,0.11,paste0("4Risk Co-Edit: n=",length(which(num==4))),col="#6633FF");
text(500,0.1,paste0("logrank p-value=",signif(pp,3)));

