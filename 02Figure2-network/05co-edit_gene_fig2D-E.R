######################################
rm(list=ls())
options(stringsAsFactors=F);
co <- read.table("D:/HFU/RNA_editing_combine_regulation/res/RADAR_CO_network_gene.txt",header=T,sep="\t");
ind <- which(co[,3]==co[,4]);
length(ind);
#5793

co1 <- sapply(co[,1],function(x)strsplit(x,";"));
co1 <- unlist(co1);
co1 <- matrix(co1,ncol=2,byrow=T);

co2 <- sapply(co[,2],function(x)strsplit(x,";"));
co2 <- unlist(co2);
co2 <- matrix(co2,ncol=2,byrow=T);

co.chr <- cbind(co1[,1],co2[,1]);
co.chr <- t(apply(co.chr,1,sort));
co.chrs <- apply(co.chr,1,function(x)paste(x,collapse="_"));

sam.ch <- which(co.chr[,1]==co.chr[,2]);
sam.Chr <- length(sam.ch)/nrow(co.chr);
sam.Chr
#[1] 0.1435901
sam.gene <- length(ind)/length(sam.ch);###同一个染色体上，44%左右同一个基因
sam.gene
# 0.4357
par(mfrow=c(1,2));
pie(c(same_chrosome=sam.Chr,different_chromosomes=1-sam.Chr))
pie(c(same_genes=sam.gene,different_genes=1-sam.gene))

setwd("D:/HFU/RNA_editing_combine_regulation/res/");

sam.chr.loc <- cbind(co1[sam.ch,2],co2[sam.ch,2]);
sam.chr.locs <- t(apply(sam.chr.loc,1,as.numeric));
sam.loc <- abs(sam.chr.locs[,1]-sam.chr.locs[,2]);
index <- which(sam.loc <= 10^7)##同一条染色体上距离小于10MB
length(index)/length(sam.loc) #0.5801

index2 <- which(sam.loc <= 10^3)##同一条染色体上距离小于1kb
length(index2)/length(sam.loc)#0.33458

s1 <- length(sam.loc[sam.loc<=10^7]);
s2 <- length(sam.loc[sam.loc>10^7 & sam.loc<= 2*10^7]);
s3 <- length(sam.loc[sam.loc>2*10^7 & sam.loc<=3*10^7]);
s4 <- length(sam.loc[sam.loc>3*10^7 & sam.loc<=4*10^7]);
s5 <- length(sam.loc[sam.loc>4*10^7 & sam.loc<=5*10^7]);
s6 <- length(sam.loc[sam.loc>5*10^7 & sam.loc<=6*10^7]);
s7 <- length(sam.loc[sam.loc>6*10^7 & sam.loc<=7*10^7]);
s8 <- length(sam.loc[sam.loc> 7 * 10^7]);
da <- c(s1/length(sam.loc),s2/length(sam.loc),s3/length(sam.loc),s4/length(sam.loc),s5/length(sam.loc),s6/length(sam.loc),s7/length(sam.loc),s8/length(sam.loc));
barplot(da,names.arg=F,ylim=c(0,0.5),ylab="Percentage of synergistic editing pairs");
mtext(at=seq(0.7,9.2,len=8),side=1,line=2,text=c("0-10Mb","10Mb-20Mb","20Mb-30Mb", "30Mb-40Mb","40Mb-50Mb","50Mb-60Mb","60Mb-70Mb",">70Mb"))
