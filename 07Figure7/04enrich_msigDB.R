rm(list=ls());
gc();
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
setwd(pat);
de <- read.table("deg_up_dn.txt",header=F,as.is=T);
##all gtf genes as backgroud
gtf <- read.table("D:/HFU/RNA_editing/data/GENCODE/gene_type_sort.bed",stringsAsFactors=F,header=T,sep="\t");
gens <- unique(gtf[gtf[,5]=="protein_coding",6]); #19776

for(i in c("up","dn")){
gene <- de[de[,2]==i,1];
gene <- intersect(gene,gens);#

path <- "D:/HFU/RNA_editing_combine_regulation/data/download/msigDB/";
msig.go <- readLines(paste0(path,"c2.cp.kegg.v2022.1.Hs.symbols.gmt"));
msig.go <- strsplit(msig.go,"\t")
go.name <- sapply(msig.go,function(x){x[[1]][1]})
go.name <- tolower(go.name)
go.name <- strsplit(go.name,"_")
go.name <- sapply(go.name,function(x){paste(x,collapse=" ")})##

res <- c()
pp <- c()
qq <- c()
enri.ge.h <- c();

for(g in 1:length(msig.go)){
kegg <- msig.go[[g]][-c(1:2)] 
q <- length(intersect(kegg,gene)) 
if(q >= 3 ){
m <- length(kegg)
n <- length(gens) - m
k <- length(gene) 
p.value <- phyper(q-1,m,n,k,lower.tail=F)  #P（X>q） 
res <- c(res,go.name[g])
pp <- c(pp,p.value)
qq <- c(qq,q)
ge.h <- intersect(kegg,gene);
go.ge.h <- paste(ge.h,collapse=";");
enri.ge.h <- c(enri.ge.h,go.ge.h);
}
}
res.pp <- data.frame(res,pp,p.adjust(pp,"BH"),qq,enri.ge.h);
if(nrow(res.pp)!=0 ){
if(length(which(p.adjust(pp,"BH")<0.05))>0){ #fdr0.05
ind <- which(p.adjust(pp,"BH")<0.05)
res.pp <- res.pp[ind,]#fdr<0.1
com.go <- res.pp[order(res.pp[,3]),];
write.table(com.go,paste0("./deg/kegg_enrich_",i,".txt"),sep="\t",row.names=F,quote=F);
}
}
}
