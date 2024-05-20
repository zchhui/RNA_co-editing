pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/validate/GSE164359_gene_expression.txt",sep="\t",header=T,as.is=T);
sam <- read.table("D:/HFU/RNA_editing_combine_regulation/data/validate/GSE164359_sam_info.csv",header=T,as.is=T,sep=",");
gene <- read.table("co_tar_deg_gene.txt",header=T,sep="\t",stringsAsFactors=F); 

exp.g <- exp[exp[,1]%in%gene[,1],];
rownames(exp.g) <- exp.g[,1];
exp.g <- t(exp.g[,-1]);
rownames(exp.g) <- gsub("FPKM.","",rownames(exp.g))

group <-  ifelse(grepl("tumor",sam[,4]),"normal","cancer");
names(group) <- sam[,1];
exp.g <- as.data.frame(exp.g);

exp.g$group = group[rownames(exp.g)];
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
N5 <- length(exp.g$group);
p5 <- performance(pred2,'auc')@y.values
p5 <- round(p5[[1]],2)
perf5<-performance(pred2,'tpr','fpr')


sam.p <- sam[sam[,3]=="tumor type: primary",];
exp.g.p <- exp.g[rownames(exp.g)%in%sam.p[,1],];
log1 <- glm(group ~.,family=binomial(link='logit'),data=exp.g.p)
summary(log1)
log2 <- step(log1)
summary(log2)
pred<-predict(log2)
prob<-exp(pred)/(1+exp(pred))
yhat<-1*(prob>0.5)
table(exp.g.p$group,yhat)
library(ROCR)
pred2<-prediction(pred,exp.g.p$group)
N5.1 <- length(exp.g.p$group);
p5.1 <- performance(pred2,'auc')@y.values
p5.1 <- round(p5.1[[1]],2)
perf5.1 <- performance(pred2,'tpr','fpr')


sam.r <- sam[sam[,3]=="tumor type: Recurrent",];
exp.g.r <- exp.g[rownames(exp.g)%in%sam.r[,1],];
log1 <- glm(group ~.,family=binomial(link='logit'),data=exp.g.r)
summary(log1)
log2 <- step(log1)
summary(log2)
pred<-predict(log2)
prob<-exp(pred)/(1+exp(pred))
yhat<-1*(prob>0.5)
table(exp.g.r$group,yhat)
library(ROCR)
pred2<-prediction(pred,exp.g.r$group)
N5.2 <- length(exp.g.r$group);
p5.2 <- performance(pred2,'auc')@y.values
p5.2 <- round(p5.2[[1]],2)
perf5.2 <- performance(pred2,'tpr','fpr')
