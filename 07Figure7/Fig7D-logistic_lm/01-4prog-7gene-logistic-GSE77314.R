pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
exp <- read.csv("D:/HFU/RNA_editing_combine_regulation/data/validate/GSE77314_expression.csv",header=T,as.is=T);
gene <- read.table("co_tar_deg_gene.txt",header=T,sep="\t",stringsAsFactors=F); 
exp.g <- exp[exp[,1]%in%gene[,1],]
rownames(exp.g) <- exp.g[,1];
exp.g <- t(exp.g[,-1])
group <-  ifelse(grepl("N",rownames(exp.g)),"normal","cancer");
exp.g <- as.data.frame(exp.g);
exp.g$group = group;
exp.g$group <- as.factor(exp.g$group);

log1 <- glm(group ~.,family=binomial(link='logit'),data=exp.g)
summary(log1)
log2 <- step(log1)
summary(log2)
pred<-predict(log2)
prob<-exp(pred)/(1+exp(pred))
yhat<-1*(prob>0.5)
table(exp.g$group,yhat)
library(ROCR)
pred2<-prediction(pred,exp.g$group)
N4 <- length(exp.g$group);
p4 <- performance(pred2,'auc')@y.values
p4 <- round(p4[[1]],2)
perf4<-performance(pred2,'tpr','fpr')
plot(perf4,col=colorss[4],type="l",lwd=2,main="GSE77314",add=T)

