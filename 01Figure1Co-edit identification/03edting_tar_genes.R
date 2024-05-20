########################################03 target genes regulated by RNA editing
rm(list=ls())
options(stringsAsFactors=F);
exp <- read.table("/NFSdata/CJ/CO_edit/data/LIHC.htseq_fpkm_symbol_example.txt",header=T);
exp.s <- colnames(exp);
exp.s <- substr(exp.s,1,15)
colnames(exp) <- exp.s;


#
ed <- read.table("/NFSdata/CJ/CO_edit/data/RADAR_edit.txt",header=T,sep="\t");
co.si <- read.table("/NFSdata/CJ/CO_edit/res/RADAR_edit_CO_fdr05_per1000_fdr05.txt",sep="\t",header=T);
si <- union(co.si[,1],co.si[,2]);
ed <- ed[si,];


### HCC tumor samples
can <- which(substr(colnames(ed),14,15)=="01" | substr(colnames(ed),14,15)=="02");
can.s.ed <- colnames(ed)[can];
can.s <- intersect(colnames(exp),can.s.ed)

ed <- ed[,can.s]
exp <- exp[,can.s]


###
for(i in 1:nrow(ed)){
res <- matrix(ncol=6,nrow=0); #c("ed","gene","p","fc","ed_num","noed_num")
ed.sam <- colnames(ed)[!is.na(ed[i,])];
noed.sam <- colnames(ed)[is.na(ed[i,])];
syms <- rownames(exp);
for(j in 1:nrow(exp)){
sym <- syms[j];
 exp.ed <- exp[sym,ed.sam];
 exp.noed <- exp[sym,noed.sam];
 if(length(ed.sam)>=3 & length(noed.sam)>=3){
 fc <- mean(as.numeric(exp.ed))/mean(as.numeric(exp.noed));
 p <- wilcox.test(as.numeric(exp.ed),as.numeric(exp.noed))$p.value;
 ed.num <- length(ed.sam);
 noed.num <- length(noed.sam);
 res <- rbind(res,c(rownames(ed)[i],sym,p,fc,ed.num,noed.num));
 }
 }
print(i) ;
write.table(res,"/NFSdata/CJ/CO_edit/res/edit_targets.txt",sep="\t",quote=F,row.names=F,col.names=F,append=T);
}

#FDR <0.05, FC>=2 filter target genes
ed.t <- read.table("/NFSdata/CJ/CO_edit/res/edit_targets.txt",header=F)
fdr <- p.adjust(ed.t[,3],"BH");
fc <- ed.t[,4]
index <- which(fdr<0.05 &  abs(log2(fc))>=1);####FC>=2 &fdr<0.05
ed.tar <- cbind(ed.t[index,],fdr[index]);
write.table(ed.tar,"/NFSdata/CJ/CO_edit/res/RADAR_targets_fdr05_fc2.txt",sep="\t",quote=F,row.names=F,col.names=F);




########################################editing-associated genes adjust by age, gender, BMI,tumor grades and stages.
rm(list=ls())
options(stringsAsFactors=F);
exp <- read.table("/NFSdata/CJ/CO_edit/data/LIHC.htseq_fpkm_symbol_example.txt",header=T);
exp.s <- colnames(exp);
exp.s <- substr(exp.s,1,15)
colnames(exp) <- exp.s;

ed <- read.table("/NFSdata/CJ/CO_edit/data/RADAR_edit.txt",header=T,sep="\t");
ed2tar <- read.table("/NFSdata/CJ/CO_edit/res/RADAR_targets_fdr05_fc2.txt",sep="\t",header=F);

### HCC tumor samples
can <- which(substr(colnames(ed),14,15)=="01" | substr(colnames(ed),14,15)=="02");
can.s.ed <- colnames(ed)[can];
can.s <- intersect(colnames(exp),can.s.ed)

ed <- ed[,can.s]
exp <- exp[,can.s]

cli <- read.csv("/NFSdata/CJ/CO_edit/res/clinical.csv",header=T,stringsAsFactors=F);
rownames(cli) <- cli[,1];
colnames(exp) <- substr(colnames(exp),1,12)
cli$stage.T[cli$stage.T=="T2a" | cli$stage.T=="T2b"] <- "T2";
cli$stage.T[cli$stage.T=="T3a" | cli$stage.T=="T3b"] <- "T3";
ind1 <- which(cli$stage.T=="TX");##去掉无法评估的样本；
cli <- cli[-ind1,];
sam <- intersect(colnames(exp),rownames(cli))
exp <- exp[,sam]
cli <- cli[sam,c(1,6:9,12)]
colnames(ed) <- substr(colnames(ed),1,12)
ed <- ed[,sam]

###
res <- data.frame(ed2tar,cof=NA,p.glm=NA)
for(i in 1:nrow(ed2tar)){
si <- ed2tar[i,1]
ed.sam <- colnames(ed)[!is.na(ed[si,])];
noed.sam <- colnames(ed)[is.na(ed[si,])];
sym <- ed2tar[i,2];
clic <- cli[c(ed.sam,noed.sam),-1]
clic$gene <- as.numeric(exp[sym,c(ed.sam,noed.sam)]);
group <- c(rep("ed",length(ed.sam)),rep("noed",length(noed.sam)))
clic$group <- group;
clic$group <- as.factor(clic$group)
log1 <- glm(group ~.,family=binomial(link='logit'),data=clic)
res[i,8] <-  summary(log1)$coefficients["gene",1]
res[i,9] <-  summary(log1)$coefficients["gene",4]
ed.num <- length(ed.sam);
noed.num <- length(noed.sam);
}
fdr <- p.adjust(res[,9],"BH")
resu <- data.frame(res,fdr) 
index <- which(fdr<0.05)
sig <- resu[index,]
write.table(sig,"/NFSdata/CJ/CO_edit/res/RADAR_targets_fdr05_fc2_adjust_fdr05.txt",sep="\t",quote=F,row.names=F,col.names=F);

