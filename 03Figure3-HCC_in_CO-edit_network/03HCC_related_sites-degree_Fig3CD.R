#HCC
par(mfrow=c(1,3));
setwd("D:/HFU/RNA_editing_combine_regulation/data/");
hcc <- read.table("RADAR_edit_type_hcc.txt",sep="\t",header=T,stringsAsFactors=F)
de.g <- read.table("../res/RADAR_co_net_degree.txt",as.is=T,heade=F)
de.hcc <- data.frame(de.g,Node.att = "Other")
ind <- which(de.hcc[,1]%in%hcc[,5])
de.hcc[ind,ncol(de.hcc)] <- "HCC"
colnames(de.hcc) <- c("site","degree","Node.att")
table(de.hcc[,3])

w.p <- wilcox.test(degree~Node.att,data=de.hcc,alternative="greater")
p <- w.p$p.value;
boxplot(degree~Node.att,data=de.hcc,outline=F,names=c("HCC_related_editing_sites","Other nodes"),
        col=c("darkred","grey"),ylab="degree")
text(2,250,paste0("Wilcox rank-sum test","\n","p=",signif(p,2)));

#hub,top10% nodes
deg <- rev(sort(de.hcc[,2]))
nm <- round(nrow(de.hcc)*0.1)
hub <- de.hcc[de.hcc[,2]>=deg[nm],];


hcc.
