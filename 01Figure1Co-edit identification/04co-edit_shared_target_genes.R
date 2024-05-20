#Co-edit network was constructed
rm(list=ls());
ed.t <- "/NFSdata/CJ/CO_edit/res/RADAR_targets_fdr05_fc2_adjust_fdr05.txt";
co <- "/NFSdata/CJ/CO_edit/res/RADAR_edit_CO_fdr05_per1000_fdr05.txt";
ed.t <- read.table(ed.t,header=F,as.is=T);
co <- read.table(co,header=T,as.is=T);
M <- length(unique(ed.t[,2]));
hyp.p <- matrix(NA,ncol=4,nrow=nrow(co));

for(i in 1:nrow(co)){ 
tar1 <- ed.t[ed.t[,1]==co[i,1],2];
tar2 <- ed.t[ed.t[,1]==co[i,2],2];
q <- length(intersect(tar1,tar2));
if(q >= 3){
ta1 <- length(tar1);
nota1 <- M - ta1;
ta2 <- length(tar2);
p.va1 <- phyper(q-1, ta1, nota1, ta2, lower.tail = F);
res <- c(q,ta1,ta2,p.va1);
hyp.p[i,] <- res; 
}
 }
co.r <- data.frame(co,hyp.p);
co.res <- co.r[!is.na(hyp.p[,1]),];
#write.table(co.res,"/NFSdata/CJ/CO_edit/res/co_tar_pval.txt",sep="\t",quote=F,row.names=F,col.names=F);

fdr <- p.adjust(co.res[,ncol(co.res)],"BH");
co.re <- data.frame(co.res,fdr);
co.r.sig <- co.re[fdr<0.05,];
colnames(co.r.sig) <- c("site1","site2","N_co-edit","BH-p.val-phyper","p.val-per","BH-P.val-perm","N.shared.tar","N.tar1","N.tar2","p.val-target","BH-P.val-target");
write.table(co.r.sig,"/NFSdata/CJ/CO_edit/res/RADAR_co_tar_fdr05.txt",sep="\t",quote=F,row.names=F,col.names=T);
co <- co.r.sig[,1:3]
write.table(co,"/NFSdata/CJ/CO_edit/res/RADAR_CO_network.txt",sep="\t",quote=F,row.names=F,col.names=T);



