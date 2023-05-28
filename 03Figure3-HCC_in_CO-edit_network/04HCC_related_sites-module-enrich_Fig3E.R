#Fig3E
rm(list=ls())
setwd("D:/HFU/RNA_editing_combine_regulation/data/");
mod <- read.table("../res/mcode_net/mcode_K3.txt",header=F,as.is=T,sep="\t");
module <- paste(mod[,5],collapse=", ")
modu <- strsplit(module,", ")[[1]];

de.hcc <- read.csv("../res/degree_HCC_Attr.csv",as.is=T,heade=T) #degree & HCC related nodes
ind.mod <- which(de.hcc[,1]%in%modu)
de.hcc <- data.frame(de.hcc,Mode_attr="Other")
de.hcc[ind.mod,"Mode_attr"] <- "Module_nodes" 

#module enrich
hcc <- which(de.hcc[,3]=="HCC")
hcc.modu <- which(de.hcc[,3]=="HCC" & de.hcc[,4]=="Module_nodes") #HCC & module nodes
nohcc <- which(de.hcc[,3]!="HCC") 
nohcc.modu <- which(de.hcc[,3]!="HCC"  &  de.hcc[,4]=="Module_nodes") 
hcc.N <- length(hcc)  
nohcc.N <- length(nohcc) 
hcc.modu.N <- length(hcc.modu)  
nohcc.modu.N <- length(nohcc.modu) 
z <- matrix(c(hcc.modu.N,nohcc.modu.N,(hcc.N-hcc.modu.N),(nohcc.N-nohcc.modu.N)),nrow=2)
hcc_modu <- fisher.test(z)
p.modu <- hcc_modu$p.value 
barplot(c(hcc.modu.N/hcc.N,nohcc.modu.N/nohcc.N),col=c("darkred","grey"),ylim=c(0,1),
        names.arg=c("HCC_related_editing_sites","Other nodes"),
        ylab="Percentage of module nodes",cex.lab=1.5)
text(c(1,0.6),paste0("Fisher's test","\n","p=",signif(p.modu,2)))
write.csv(de.hcc,"../res/degree_HCC_module_Attr.csv",row.names=F,quote=F)

#z
#   [,1] [,2]
#[1,]   111  88
#[2,] 1091 5859

