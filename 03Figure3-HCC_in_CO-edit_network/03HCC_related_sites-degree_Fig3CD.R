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


hcc.Web=intersect(de.hcc[,1],hcc[,5]) 
hcc.hub=intersect(hub[,1],hcc[,5]) 
nhcc.Web <- setdiff(de.hcc[,1],hcc[,5]);
nhcc.hub <- intersect(nhcc.Web,hub[,1]) 
hccweb.num=length(hcc.Web)  
nhccweb.num=length(nhcc.Web) 
hcc.hub.num=length(hcc.hub)  
nhcc.hub.num=length(nhcc.hub) 
z=matrix(c(hcc.hub.num,nhcc.hub.num,(hccweb.num-hcc.hub.num),
         (nhccweb.num-nhcc.hub.num)),nrow=2)
hcc_hub=fisher.test(z)
p.hub <- hcc_hub$p.value 
barplot(c(hcc.hub.num/hccweb.num,nhcc.hub.num/nhccweb.num),
        col=c("darkred","grey"),ylim=c(0,1),names.arg=c("HCC_related_editing_sites","Other nodes"),
        ylab="Percentage of hub nodes",cex.lab=1.5)
text(c(1,0.6),paste0("Fisher's test","\n","p=",signif(p.hub,2)))
write.csv(de.hcc,"../res/degree_HCC_Attr.csv",row.names=F,quote=F)
z
#     [,1] [,2]
#[1,]   96  103
#[2,]  619 6331
