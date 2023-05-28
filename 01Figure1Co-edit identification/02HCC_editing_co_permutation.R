#02Random perturbations are used to identify candidate RNA co-editing pair
rm(list=ls());
gc();
library("vegan")
pa1 <- "/NFSdata/CJ/CO_edit/data/";
pa2 <- "/NFSdata/CJ/CO_edit/res/";
options(stringsAsFactors = F)

IN1 <- paste0(pa1,"RADAR_edit.txt");
IN2 <- paste0(pa2,"RADAR_CO_fdr05.txt");
OUT1 <- paste0(pa2,"permutation1000.RData");
OUT2 <- paste0(pa2,"RADAR_edit_CO_fdr05_per1000_fdr05.txt")

edi <- read.table(IN1,sep="\t",header=T);
co <- read.table(IN2,sep="\t",header=T);

#0-1 matrix
ed <- apply(edi,1,function(x){
 x[!is.na(x)] <- 1;
 x[is.na(x)] <- 0;
 return(x);
 })
ed <- t(ed);
out <- permatswap(ed,fixedmar="both", shuffle="both", times = 1000, mtype = "prab");
save(out,OUT)

###site1 edit  site2 edit
edi.com <- function(x,y){
 ind1 <- which(x==1);
 ind2 <- which(y==1);
 ed.com <- length(intersect(ind1,ind2));
 return(ed.com);
 }

re <- matrix(0,ncol=1000,nrow=nrow(co));
resu <- cbind(co[,c(1:3,7)],re); #site1 ,site2, N_co-edit, phyper_fdr

for(p in 1:1000){
ed.per <- out$perm[[p]];
##co-occur
for(i in 1:nrow(co))
 { 
res1 <- edi.com(ed.per[co[i,1],],ed.per[co[i,2],])
resu[i,p+4] <- res1 
  }
}

res.f <- resu[,1:5];
for(i in 1:nrow(resu)){
  num.h <- which(resu[i,5:1004]>=resu[i,3]);
  p.val <- length(num.h)/1000;
  res.f[i,5] <- p.val;
}
fdr <- p.adjust(res.f[,5],"BH")
res.f <- cbind(res.f,fdr);
colnames(res.f) <- c("site1","site2","N_co-edit","BH-p.val-phyper","p.val-per","BH-P.val-perm");
res.sig <- res.f[fdr<0.05,];#permutation 1000 times,FDR<0.05
write.table(res.sig,OUT2,sep="\t",col.names=T,row.names=F,quote=F);

