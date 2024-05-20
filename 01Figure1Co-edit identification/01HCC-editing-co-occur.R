##01 Identification of potential RNA co-editing pairs by hypergeometric distribution
rm(list=ls());
gc();
pat <- "/NFSdata/CJ/CO_edit/data/";
pat2 <- "/NFSdata/CJ/CO_edit/res/";

options(stringsAsFactors = F)
setwd(pat);

edit <- read.table(paste0(pat,"RADAR_edit.txt"),sep="\t",header=T);

##hcc
can <- which(substr(colnames(edit),14,15)=="01" | substr(colnames(edit),14,15)=="02");
edi <- edit[,can];

#co-edited in at least 3 HCC samples
edit.num <- apply(edi,1,function(x){
  sum(as.numeric(!is.na(x)))}
  )
edi <- edi[edit.num>=3,] 
  
###site1 edit  site2 edit
edi.com.p <- function(x,y){
 ind1 <- which(!is.na(x));
 ind2 <- which(!is.na(y));
 ed.com <- length(intersect(ind1,ind2));
 ed1 <- length(ind1);
 noed1 <- length(x)-ed1;
 ed2 <- length(ind2);
 p.val <- phyper(ed.com-1, ed1, noed1, ed2, lower.tail = F); ##phyper test
 res <- c(ed.com,ed1,ed2,p.val);
 return(res);
 }

result <- matrix(nrow=0,ncol=6);
for (i in 1:(nrow(edi)-1)){
  for(j in (i+1):nrow(edi)){
site1 <- rownames(edi)[i]
site2 <- rownames(edi)[j]
res1 <- edi.com.p(edi[i,],edi[j,])
res <- c(site1,site2,res1)
result <- rbind(result,res)
  }
}
colnames(result) <- c("site1","site2","N_co-edit","N_ed1","N_ed2","phyper_p")
fdr <- p.adjust(result[,6],"BH")
result <- data.frame(result,phyper_fdr=fdr)
result <- result[fdr<0.05,]
write.table(result,paste0(pat2,"RADAR_CO_fdr05.txt"),sep="\t",col.names=T,row.names=F,quote=F);

