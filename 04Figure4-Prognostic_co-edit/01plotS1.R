####identify prognostic related co-editing pairs 
rm(list=ls());
gc();
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
cli <- read.csv("./prog/clinical.csv",header=T,stringsAsFactors=F);
rownames(cli) <- cli[,1];
mt <- read.table("co_edit_sample_matrix_cancer.txt",header=T,sep="\t",stringsAsFactors=F); 
mtr <- mt[,substr(colnames(mt),1,12)%in%rownames(cli)]; ##shared HCC samples in follow-up data

mt <- cbind(mt[,1:2],mtr);
num <- apply(mtr,1,function(x){sum(x==1)});
cuf <- round(0.05*ncol(mtr));        ## co-edited in at least 5% percent HCC samples 
index <- which(num>=cuf);
mtr <- mtr[index,];
mt <- mt[index,];
cli <- cli[intersect(substr(colnames(mt),1,12),rownames(cli)),];

res <- c();
HR <- c();
for(i in 1:nrow(mt)){
 ed.s <- colnames(mtr)[mtr[i,]==1];
 noed.s <- colnames(mtr)[mtr[i,]==0];
clic <- cbind(cli,cli[,1]);
clic[substr(ed.s,1,12),ncol(clic)] <- "Co-edit";  ##co.ed=1
clic[substr(noed.s,1,12),ncol(clic)] <- "No_co-edit";  ###noco.ed=0
colnames(clic)[ncol(clic)] <- "co_ed";
bancanshu <- coxph(Surv(time,status)~co_ed,data=clic);
coef <- summary(bancanshu)$ coef
beta <- coef[1,1];
hz <- coef[1,2];
p <- coef[1,5];
res <- c(res,p);
HR <- c(HR,hz);
  }
fdr <- p.adjust(res,"BH");
resu <- cbind(mt[,1:2],HR,res,fdr);
colnames(resu) <- c("co1","co2","hazard ratio","surv.p","surv.fdr")
#write.table(resu,"CO_prog_fdr_new.txt",sep="\t",col.names=T,row.names=F,quote=F);
result <- resu[fdr<0.05,];
typ <- read.table("D:/HFU/RNA_editing_combine_regulation/data/RADAR_edit_type.txt",header=F,as.is=T);
typ <- typ[,c(3,5,6)];
colnames(typ) <- c("co1","ge.ty.co1","ge.co1");
re <- merge(result,typ,by="co1");
colnames(typ) <- c("co2","ge.ty.co2","ge.co2");
re2 <- merge(re,typ,by="co2");
re2 <- re2[,c(2,1,3:ncol(re2))];
write.table(re2,"CO_prog_fdr05_5per_new.txt",sep="\t",col.names=T,row.names=F,quote=F);


################## plot Figure S1
re2 <- read.table("D:/HFU/RNA_editing_combine_regulation/res/CO_prog_fdr05_5per_new.txt",header=T,as.is=T,sep="\t")
par(mfrow=c(3,4))
for(m in 1:nrow(re2)){
index <- which(mt[,1]==re2[m,1] & mt[,2]==re2[m,2]);
i <- index;
ed.s <- colnames(mtr)[mtr[i,]==1];
 noed.s <- colnames(mtr)[mtr[i,]==0];
clic <- cbind(cli,cli[,1]);
clic[substr(ed.s,1,12),ncol(clic)] <- "Co-edit";  ##co.ed=1
clic[substr(noed.s,1,12),ncol(clic)] <- "No_co-edit";  ###noco.ed=0
colnames(clic)[ncol(clic)] <- "co_ed";
prog.dif <- survdiff(Surv(time,status)~co_ed,data=clic); #logrank
prog.fit <- survfit(Surv(time,status)~co_ed,data=clic);
p <- prog.dif$pvalue

plot(prog.fit,col=c("darkblue","darkred"),lwd=2,xlab="Survival time(days)",ylab="Survival(%)",main=paste0(paste(re2[m,1:2],collapse="-"),"\n",paste(re2[m,c(7,9)],collapse="-")))
legend(2000,0.9,c("No co-edit","Co-edit"),col=c("darkblue","darkred"),lty=1,text.col=c("darkblue","darkred"),lwd=2)
text(2500,0.08,paste0("Noco-edit: n=",length(which(clic[,ncol(clic)]==0))),col="darkblue")
text(2500,0.03,paste0("Co-edit: n=",length(which(clic[,ncol(clic)]==1))),col="darkred")
text(2500,1,paste0("logrank p = ",signif(p,3)),lwd=2);
}

