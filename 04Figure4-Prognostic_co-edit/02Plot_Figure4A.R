###Plot Figure 4A
rm(list=ls());
gc();
library(survival);
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
cli <- read.csv("prog/clinical.csv",header=T,stringsAsFactors=F);
rownames(cli) <- cli[,1];
mt <- read.table("co_edit_sample_matrix_cancer.txt",header=T,sep="\t",stringsAsFactors=F); 
mtr <- mt[,substr(colnames(mt),1,12)%in%rownames(cli)]; ####shared HCC samples in follow-up data

mt <- cbind(mt[,1:2],mtr);
num <- apply(mtr,1,function(x){sum(x==1)});
cuf <- round(0.05*ncol(mtr));## co-edited in at least 5% percent HCC samples 
index <- which(num>=cuf);
mtr <- mtr[index,];
mt <- mt[index,];
result <- read.table("CO_prog_fdr05_5per_new.txt",sep="\t",header=T,as.is=T);
edi <- read.table("D:/HFU/RNA_editing_combine_regulation/data/RADAR_edit.txt",header=T,as.is=T);
sam <- intersect(colnames(edi),colnames(mt)); 
edi <- edi[,sam];
cli <- cli[intersect(substr(colnames(mt),1,12),rownames(cli)),];


p.val <- c();
par(mfrow=c(2,3));
for(m in 1:nrow(result)){
 index <- which(mt[,1]==result[m,1] & mt[,2]==result[m,2]);
 i <- index;
 ed.s <- colnames(mtr)[mtr[i,]==1];
 noed.s <- colnames(mtr)[mtr[i,]==0];
 clic <- cbind(cli,cli[,1]);
 clic[substr(ed.s,1,12),ncol(clic)] <- "Co-edit";  ##co.ed
 clic[substr(noed.s,1,12),ncol(clic)] <- "No_co-edit";  ###noco.ed
 colnames(clic)[ncol(clic)] <- "co_ed";
 ###single edit sites site1
 clic1 <- cbind(clic,clic[,1]);
 ed.s1 <- colnames(edi)[!is.na(edi[mt[i,1],])];
 noed.s1 <- colnames(edi)[is.na(edi[mt[i,1],])];
 clic1[substr(ed.s1,1,12),ncol(clic1)] <- 1;  ##ed1=1
 clic1[substr(noed.s1,1,12),ncol(clic1)] <- 0;  ###noed1=0
 colnames(clic1)[ncol(clic1)] <- "ed1";
 ######single editing sites：site2
 clic2 <- cbind(clic1,clic1[,1]);
 ed.s2 <- colnames(edi)[!is.na(edi[mt[i,2],])];
 noed.s2 <- colnames(edi)[is.na(edi[mt[i,2],])];
 clic2[substr(ed.s2,1,12),ncol(clic2)] <- 1;  ##ed2=1
 clic2[substr(noed.s2,1,12),ncol(clic2)] <- 0;  ###noed2=0
 colnames(clic2)[ncol(clic2)] <- "ed2";

 bancanshu <- coxph(Surv(time,status)~co_ed+ed1+ed2,data=clic2);
 coef <- summary(bancanshu)$ coef
 hr <- coef[1,2];
 p <- coef[1,5];##co-edit
 p.val <- c(p.val,p);
}
fdr <- p.adjust(p.val,"BH");
res <- data.frame(result,ed.adj.p.val=p.val,ed.adj.fdr=fdr);
write.table(res,"CO_prog_fdr05_5per_new_adjust_single_ed.txt",sep="\t",col.names=T,row.names=F,quote=F);



###plot 4A
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
res <- read.table("CO_prog_fdr05_5per_new_adjust_single_ed.txt",header=T,as.is=T,sep="\t")
par(mfrow=c(2,3));
mm <- which(res[,"ed.adj.fdr"]<0.05)
result <- res[mm,]
for(m in 1:5){
index <- which(mt[,1]==result[m,1] & mt[,2]==result[m,2]);
i <- index;
ed.s <- colnames(mtr)[mtr[i,]==1];
noed.s <- colnames(mtr)[mtr[i,]==0];
clic <- cbind(cli,cli[,1]);
clic[substr(ed.s,1,12),ncol(clic)] <- 1;  ##co.ed=1
clic[substr(noed.s,1,12),ncol(clic)] <- 0;  ###noco.ed=0
colnames(clic)[ncol(clic)] <- paste0("co_ed",m);
###single editing sites：site1
clic1 <- cbind(clic,clic[,1]);
ed.s1 <- colnames(edi)[!is.na(edi[mt[i,1],])];
noed.s1 <- colnames(edi)[is.na(edi[mt[i,1],])];
clic1[substr(ed.s1,1,12),ncol(clic1)] <- 1;  ##ed1=1
clic1[substr(noed.s1,1,12),ncol(clic1)] <- 0;  ###noed1=0
colnames(clic1)[ncol(clic1)] <- "ed1";
###single editing sites：site2
clic2 <- cbind(clic1,clic1[,1]);
ed.s2 <- colnames(edi)[!is.na(edi[mt[i,2],])];
noed.s2 <- colnames(edi)[is.na(edi[mt[i,2],])];
clic2[substr(ed.s2,1,12),ncol(clic2)] <- 1;  ##ed2=1
clic2[substr(noed.s2,1,12),ncol(clic2)] <- 0;  ###noed2=0
colnames(clic2)[ncol(clic2)] <- "ed2";

prog.dif <- survdiff(Surv(time,status)~co_ed+ed1+ed2,data=clic2);
prog.fit <- survfit(Surv(time,status)~co_ed+ed1+ed2,data=clic2);
plot(prog.fit,col=c("grey","orange","chocolate","darkred"),lwd=2,xlab="Survival time(days)",ylab="Survival(%)",main=paste0(paste(result[m,1:2],collapse="-"),"\n",paste(result[m,c(7,9)],collapse="-")))
legend(2000,1,c("No_Edit","Edit2","Edit1","co-edit"),col=c("grey","orange","chocolate","darkred"),lty=1,text.col=c("grey","orange","chocolate","darkred"),lwd=2)
text(2500,0.24,paste0("No_Edit: n=",length(which(clic2[,ncol(clic)]==0 & clic2[,ncol(clic1)]==0 & clic2[,ncol(clic2)]==0 ))),col="grey")
text(2500,0.1,paste0("Edit2: n=",length(which(clic2[,ncol(clic)]==0 & clic2[,ncol(clic1)]==0 & clic2[,ncol(clic2)]==1 ))),col="orange")
text(2500,0.17,paste0("Edit1: n=",length(which(clic2[,ncol(clic)]==0 & clic2[,ncol(clic1)]==1 & clic2[,ncol(clic2)]==0 ))),col="chocolate")
text(2500,0.03,paste0("co-edit: n=",length(which(clic2[,ncol(clic)]==1 & clic2[,ncol(clic1)]==1 & clic2[,ncol(clic2)]==1 ))),col="darkred")
pp <- signif(prog.dif$pvalue,3)
text(500,0.1,paste0("logrank p-value=",pp));
}

