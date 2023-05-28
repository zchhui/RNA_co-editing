###Table S2: The results of multivariate survival analysis of 5 prognostic RNA co-editing pairs.
rm(list=ls());
gc();
library(survival);
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
cli <- read.csv("prog/clinical.csv",header=T,stringsAsFactors=F);
rownames(cli) <- cli[,1];
mt <- read.table("co_edit_sample_matrix_cancer.txt",header=T,sep="\t",stringsAsFactors=F); 
mtr <- mt[,substr(colnames(mt),1,12)%in%rownames(cli)]; ##临床有关的mt
mt <- cbind(mt[,1:2],mtr);
num <- apply(mtr,1,function(x){sum(x==1)});
cuf <- round(0.05*ncol(mtr));##至少5%样本出现协同编辑，不是罕见现象
index <- which(num>=cuf);
mtr <- mtr[index,];
mt <- mt[index,];
result <- read.table("CO_prog_fdr05_5per_new_adjust_single_ed.txt",sep="\t",header=T,as.is=T);
result <- result[result[,"ed.adj.fdr"]<0.05,]
edi <- read.table("D:/HFU/RNA_editing_combine_regulation/data/RADAR_edit.txt",header=T,as.is=T);
sam <- intersect(colnames(edi),colnames(mt)); ###
edi <- edi[,sam];
cli <- cli[intersect(substr(colnames(mt),1,12),rownames(cli)),];
cli$stage.T[cli$stage.T=="T2a" | cli$stage.T=="T2b"] <- "T2";
cli$stage.T[cli$stage.T=="T3a" | cli$stage.T=="T3b"] <- "T3";
ind1 <- which(cli$stage.T=="TX");##去掉无法评估的样本；
cli <- cli[-ind1,];

p.val <- c();
par(mfrow=c(2,3));
for(m in 1:nrow(result)){
 index <- which(mt[,1]==result[m,1] & mt[,2]==result[m,2]);
 i <- index;
 ed.s <- colnames(mtr)[mtr[i,]==1];
 noed.s <- colnames(mtr)[mtr[i,]==0];
 clic <- cbind(cli,cli[,1]);
 clic[substr(ed.s,1,12),ncol(clic)] <- "Yes";  ##co.ed=1
 clic[substr(noed.s,1,12),ncol(clic)] <- "No";  ###noco.ed=0
 colnames(clic)[ncol(clic)] <- "co_ed";
 
 bancanshu <- coxph(Surv(time,status)~co_ed+stage.T+grade+ages+gender+fetoprotein_outcome_value,data=clic);
 summary(bancanshu)
 coef <- summary(bancanshu)$ coef
 hr <- coef[1,2];
 p <- coef[1,5];
 p.val <- c(p.val,p);
}
fdr <- p.adjust(p.val,"BH");
res <- data.frame(result,ed.adj.p.val=p.val,ed.adj.fdr=fdr);
write.table(res,"CO_prog_fdr05_5per_new_adjust_single_ed_stage_grade.txt",sep="\t",col.names=T,row.names=F,quote=F);




