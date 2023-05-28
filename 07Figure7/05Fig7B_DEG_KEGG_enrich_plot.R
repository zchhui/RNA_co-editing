##plot_enrich
rm(list=ls());
gc();
par(mfcol=c(1,2))
pa.res <- "D:/HFU/RNA_editing_combine_regulation/res/deg/";
setwd(pa.res);
data <- read.table("kegg_enrich_up.txt",stringsAsFactors=F,header=T,sep="\t") 
logp <- -log2(data[,3]) ###log2FDR
data <- cbind(data[,1],logp)
barplot(as.numeric(data[,2]),col="#DE203F",horiz=T,border="NA",xlab="Enrichment(-logFDR)")
text(3,seq(0.7,5.5,length=nrow(data)),data[,1],adj=0)

data <- read.table("kegg_enrich_dn.txt",stringsAsFactors=F,header=T,sep="\t") 
logp <- -log2(data[,3]) ###log2FDR
data <- cbind(data[,1],logp)
barplot(as.numeric(data[,2]),col="#67B740",horiz=T,border="NA",xlab="Enrichment(-logFDR)")
text(3,seq(0.7,24.7,length=nrow(data)),data[,1],adj=0)

