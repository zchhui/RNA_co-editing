rm(list=ls())
library(survival);

pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
cli <- read.csv("prog/clinical.csv",header=T,stringsAsFactors=F);
rownames(cli) <- cli[,1];
mt <- read.table("co_edit_sample_matrix_cancer.txt",header=T,sep="\t",stringsAsFactors=F); 
mtr <- mt[,substr(colnames(mt),1,12)%in%rownames(cli)];
mt <- cbind(mt[,1:2],mtr);
num <- apply(mtr,1,function(x){sum(x==1)});
cuf <- round(0.05*ncol(mtr));##co-edit in at least 5% HCC samples
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
clic[substr(names(num)[num==0],1,12),ncol(clic)] <- "Low risk";  ##co.ed=0
clic[substr(names(num)[num>=1],1,12),ncol(clic)] <- "High risk";  ###have at least 1 risk co-edit pairs
colnames(clic)[ncol(clic)] <- "co_ed";

prog.dif <- survdiff(Surv(time,status)~co_ed,data=clic);
pp <- signif(prog.dif$pvalue,3)
prog.fit <- survfit(Surv(time,status)~co_ed,data=clic);
colorss <- c("#619CFF","#F8766D");
plot(prog.fit,col=colorss,lwd=2,xlab="Survival time(days)",ylab="Survival(%)")
legend(2000,1,c("High risk","Low risk"),col=colorss,lty=1,text.col=colorss,lwd=2);
text(500,0.1,paste0("logrank p-value=",pp));
