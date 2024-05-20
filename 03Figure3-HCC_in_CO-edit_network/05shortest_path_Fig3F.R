library("igraph")
rm(list=ls())
setwd("D:/HFU/RNA_editing_combine_regulation/res/");
net <- read.table("RADAR_CO_network.txt",header=T,as.is=T);
library(igraph)
g <- graph.data.frame(net[,1:2],directed=F)
de <- decompose.graph(g)[[1]]
deco.net <- get.data.frame(de)
##Max component
graph <- graph.data.frame(deco.net,directed=F)
net.node <- unique(c(deco.net[,1],deco.net[,2]))
hcc <- read.table("D:/HFU/RNA_editing_combine_regulation/data/RADAR_edit_type_hcc.txt",sep="\t",header=T,stringsAsFactors=F)
hcc <- intersect(hcc[,5],net.node)##hcc related nodes in Max component
D1 <- shortest.paths(graph, as.character(hcc), as.character(hcc))
num_pairs <- (choose(length(hcc),2))*2
cpl <- sum(D1)/ num_pairs

ind <- which(net.node%in%hcc)
oth <- net.node[-ind]
set.seed(2)
ran <- replicate(1000,{
                node.random <- sample(oth,length(ind),replace=F)
                D <- shortest.paths(graph, as.character(node.random), as.character(node.random))
                sum(D)/ num_pairs
                }
                 )
length(which(ran<cpl))/1000
write.table(ran,"D:/HFU/RNA_editing_combine_regulation/res/shortest_path/shortest_path_random1000.txt",sep="\t",col.names=F,row.names=F,quote=F)

p.value <-  length(which(ran<cpl))/1000  #0

plot(density(ran),xlim=c(2.5,4),col="grey",xlab="Mean shortest path length",cex=2,lwd=2.5,cex.axis=1.5,cex.lab=1.5,cex.main=2)
arrows(cpl,1,cpl,0, angle = 10, col = "darkred", lwd = 2)
legend("topright",c("Random nodes","HCC_realted_editing_sites"),col=c("grey","darkred"),lwd=2.5)
text(3,3,paste0("p.value=",p.value),cex=1.5)




