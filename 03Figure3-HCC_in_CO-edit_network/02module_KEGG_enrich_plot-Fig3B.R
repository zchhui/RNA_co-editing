###K3 mcode
setwd("D:/HFU/RNA_editing_combine_regulation/res/enrich/");
mod1 <- read.table("c2.cp.kegg.v6.2.symbols.gmt_mod1.txt",sep="\t",header=T,as.is=T);
mod2 <- read.table("c2.cp.kegg.v6.2.symbols.gmt_mod2.txt",sep="\t",header=T,as.is=T);
mod3 <- read.table("c2.cp.kegg.v6.2.symbols.gmt_mod3.txt",sep="\t",header=T,as.is=T);
par(mfrow=c(3,1))
barplot(abs(log2(mod1[,3])),horiz=T,col="pink",main="Enriched KEGG in module1",xlab="Enrichment(log2(adjust p value))");
text(0,seq(0.7,1.9,length=nrow(mod1)),mod1[,1],adj=0)
barplot(abs(log2(mod2[,3])),horiz=T,col="pink",main="Enriched KEGG in module2",xlab="Enrichment(log2(adjust p value))");
text(0,seq(0.7,4.3,length=nrow(mod2)),mod2[,1],adj=0)
barplot(abs(log2(mod3[,3])),horiz=T,col="pink",main="Enriched KEGG in module3",xlab="Enrichment(log2(adjust p value))");
text(0,seq(0.7,6.7,length=nrow(mod3)),mod3[,1],adj=0)



