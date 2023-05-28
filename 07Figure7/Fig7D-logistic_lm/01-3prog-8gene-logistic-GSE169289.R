pat <- "D:/HFU/RNA_editing_combine_regulation/res/";
options(stringsAsFactors = F)
setwd(pat);
exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/validate/GSE169289_Transcript_FPKM.txt",header=T,as.is=T);
id <- read.table("D:/HFU/ID/Gene_name2Refseq.txt",header=T,as.is=T,sep="\t");
gene <- read.table("8co_tar_deg_gene.txt",header=T,sep="\t",stringsAsFactors=F); 
#gene <- read.table("co_target3.txt",header=T,sep="\t",stringsAsFactors=F); 
#gene <- read.table("co_target2.txt",header=T,sep="\t",stringsAsFactors=F); 
id <- id[id[,1]%in%gene[,1],]
exp.1 <- exp[,1];
nch <- nchar(exp.1);
exp.1 <- substr(exp.1,1,(nch-2));
exp[,1] <- exp.1;


exp.g <- exp[exp[,1]%in%id[,2],];
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
N3 <- length(exp.g$group);
p3 <- performance(pred2,'auc')@y.values
p3 <- round(p3[[1]],2)
perf3<-performance(pred2,'tpr','fpr')
plot(perf3,col=colorss[3],type="l",lwd=2,main="GSE169289",add=T)

