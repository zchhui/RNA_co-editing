###K3 mcode
setwd("D:/HFU/RNA_editing_combine_regulation/res/");
ed.g <- read.table("../data/RADAR_edit_type.txt",sep="\t",header=F,stringsAsFactors=F)
mod <- read.table("./mcode_net/mcode_K3.txt",header=F,as.is=T,sep="\t");
##all gtf genes as backgroud
gtf <- read.table("D:/HFU/RNA_editing/data/GENCODE/gene_type_sort.bed",stringsAsFactors=F,header=T,sep="\t");
bg <- unique(gtf[gtf[,5]=="protein_coding",6]); #19776

for(i in 1:nrow(mod)){
module <- paste(mod[i,5],collapse=", ") #module i
modu <- strsplit(module,", ")[[1]];
ed.mod <- ed.g[ed.g[,1]%in%modu,]
gene <- unique(ed.mod[,6]);  ##gene name of editing sites located on 
#ed.g first line
#chr1;100023798 chr1;100023798;A;G chr1;100023798;A;G  + protein_coding SLC35A3
gene <- intersect(gene,bg);

path <- "D:/HFU/RNA_editing_combine_regulation/data/download/msigDB/";
msig.go <- readLines(paste0(path,"c2.cp.kegg.v2022.1.Hs.symbols.gmt"));
msig.go <- strsplit(msig.go,"\t")
go.name <- sapply(msig.go,function(x){x[[1]][1]})
go.name <- tolower(go.name)
go.name <- strsplit(go.name,"_")
go.name <- sapply(go.name,function(x){paste(x,collapse=" ")})

res <- c()
pp <- c()
qq <- c()
enri.si.ge <- c();

for(g in 1:length(msig.go)){
msig.gen <- msig.go[[g]][-c(1:2)] 
q <- length(intersect(msig.gen,gene))  #intersect genes
if(q >= 3 ){
m <- length(msig.gen)
n <- length(bg) - m
k <- length(gene) 
p.value <- phyper(q-1,m,n,k,lower.tail=F)  #P（X>q） 
res <- c(res,go.name[g])
pp <- c(pp,p.value)
qq <- c(qq,q)
ge.q <- intersect(msig.gen,gene);
si.g <- ed.mod[ed.mod[,6]%in%ge.q,c(1,6)] #sites & genes in mod i
site.gene.kegg <- apply(si.g,1,function(x)paste(x,collapse="-"))
site.gene.kegg <- paste(site.gene.kegg,collapse=",")
enri.si.ge <- c(enri.si.ge,site.gene.kegg);
}
}
res.pp <- data.frame(res,pp,p.adjust(pp,"BH"),qq,enri.si.ge);
if(nrow(res.pp)!=0 ){
if(length(which(p.adjust(pp,"BH")<0.05))>0){  
ind <- which(p.adjust(pp,"BH")<0.05) #fdr<0.05
res.pp <- res.pp[ind,]
com.go <- res.pp[order(res.pp[,3]),];
write.table(com.go,paste0("./mcode_net/enrich/kegg_enrich_mod",i,".txt"),sep="\t",row.names=F,quote=F);
}
}
}



enri.si.ge <- c();

for(g in 1:length(msig.go)){
msig.gen <- msig.go[[g]][-c(1:2)] 
q <- length(intersect(msig.gen,gene))  #intersect genes
if(q >= 3 ){
m <- length(msig.gen)
n <- length(bg) - m
k <- length(gene) 
p.value <- phyper(q-1,m,n,k,lower.tail=F)  #P（X>q） 
res <- c(res,go.name[g])
pp <- c(pp,p.value)
qq <- c(qq,q)
ge.q <- intersect(msig.gen,gene);
si.g <- ed.mod[ed.mod[,6]%in%ge.q,c(1,6)] #sites & genes in mod i
site.gene.kegg <- apply(si.g,1,function(x)paste(x,collapse="-"))
site.gene.kegg <- paste(site.gene.kegg,collapse=",")
enri.si.ge <- c(enri.si.ge,site.gene.kegg);
}
}
res.pp <- data.frame(res,pp,p.adjust(pp,"BH"),qq,enri.si.ge);
if(nrow(res.pp)!=0 ){
if(length(which(p.adjust(pp,"BH")<0.05))>0){  
ind <- which(p.adjust(pp,"BH")<0.05) #fdr<0.05
res.pp <- res.pp[ind,]
com.go <- res.pp[order(res.pp[,3]),];
write.table(com.go,paste0("./mcode_net/enrich/kegg_enrich_mod",i,".txt"),sep="\t",row.names=F,quote=F);
}
}
}
