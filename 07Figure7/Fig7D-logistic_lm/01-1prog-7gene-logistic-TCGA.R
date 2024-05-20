rm(list=ls())
pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/exp/LIHC.htseq_fpkm_symbol_90per_unique.txt",row.names=1,header=T,as.is=T);
gene <- read.table("co_tar_deg_gene.txt",header=T,sep="\t",stringsAsFactors=F); 
exp.g <- exp[gene[,1],]
exp.g <- t(exp.g)
colnames(exp.g) <- gsub("-",".",colnames(exp.g))

samp <- substr(rownames(exp.g),14,15);
group <-  ifelse(samp==11,"normal","cancer");
exp.g <- as.data.frame(exp.g);
exp.g$group <- group;
exp.g$group <- as.factor(exp.g$group)

log1 <- glm(group ~.,family=binomial(link='logit'),data=exp.g)
summary(log1)
log2 <- step(log1)
summary(log2)
pred <- predict(log2)
prob <- exp(pred)/(1+exp(pred))
yhat <- 1*(prob>0.5)
table(exp.g$group,yhat)
library(ROCR)
pred2 <- prediction(pred,exp.g$group)
N1 <- length(exp.g$group);
p1 <- performance(pred2,'auc')@y.values
p1 <- round(p1[[1]],2)
perf1 <- performance(pred2,'tpr','fpr')
colorss <- c("red","orange","blue","darkgreen","purple");
plot(perf1,col=colorss[1],type="l",lwd=2,main="TCGA")
