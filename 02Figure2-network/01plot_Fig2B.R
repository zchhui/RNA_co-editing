######################################R4.2.2
rm(list=ls())
library(igraph);
options(stringsAsFactors=F);
IN <- "D:/HFU/RNA_editing_combine_regulation/res/RADAR_CO_network.txt"
OUT <- "D:/HFU/RNA_editing_combine_regulation/res/RADAR_co_net_degree.txt"

co <- read.table(IN,header=T,sep="\t");
g <- graph.data.frame(co[,1:2], directed = F)
dg <- degree(g,  mode = "total",loops = TRUE, normalized = FALSE);
dg <- dg[rev(order(dg))];
dgs <- cbind(names(dg),dg);
write.table(dgs,OUT,sep="\t",col.names=F,row.names=F,quote=F);

#degree Fig2B
dg <- read.table(OUT,header=T,sep="\t");
degree <- table(dg[,2]);
x <- names(degree);
y <- degree;
log_x <- log10(as.numeric(x));
log_y <- log10(y);
res.lm <- lm(log_y~log_x)
cof <- summary(res.lm)
coef <- round(cof$coefficients[2,1],2) #-1.41
adj.r <- round(cof$adj.r.squared,2)  #0.91

plot(log_x,log_y,col="#9B181B",xaxt="n",yaxt="n",xlab="Degree",ylab="Number of nodes",cex=1.8,cex.lab=1.8,cex.main=2,main="Degree distribution of network");
axis(1,labels=c("0","10","100"),at=c(0,1,2),cex.axis=1.8);
axis(2,labels=c("0","10","100"),at=c(0,1,2),cex.axis=1.8);
abline(res.lm,col="#9B181B",cex=1.8);
text(2.0,2.2, expression(y == x^-1.41 ),cex=1.8)
text(2.0,2.0,expression(R^2 == 0.91),cex=1.8)

