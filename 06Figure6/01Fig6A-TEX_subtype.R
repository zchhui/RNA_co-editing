#Fig6A
rm(list=ls());
gc();
pa.res <- "D:/HFU/RNA_editing_combine_regulation/res/";
setwd(pa.res);
exp <- read.table("D:/HFU/RNA_editing_combine_regulation/data/exp/LIHC.htseq_fpkm_symbol_90per_unique.txt",row.names=1,header=T,as.is=T);
sam.gr <- read.table("HCC_sample_group_total.txt",header=T,as.is=T)
gene <- c("IL2","TNF","IL10","IFNG","TGFB1","TCF7","TOX","TBX21","CTLA4","TIGIT","LAG3","HAVCR2","PDCD1","CD274");
# exhaustic T cell related genes

exp.g <- exp[gene,];
low.r <- sam.gr[sam.gr[,2]=="Tum.S0",1]
high.r <- sam.gr[(sam.gr[,2]!="Tum.S0" & sam.gr[,2]!="Normal") ,1]
nor <- sam.gr[ sam.gr[,2]=="Normal" ,1]

pp1 <- c();
pp2 <- c();
pp3 <- c();
for(i in 1:nrow(exp.g)){
p1 <- wilcox.test(log2(as.numeric(exp.g[i,low.r])+0.05),log2(as.numeric(exp.g[i,high.r])+0.05))$p.value;
p2 <- wilcox.test(log2(as.numeric(exp.g[i,high.r])+0.05),log2(as.numeric(exp.g[i,nor])+0.05))$p.value;
p3 <- wilcox.test(log2(as.numeric(exp.g[i,low.r])+0.05),log2(as.numeric(exp.g[i,nor])+0.05))$p.value;
pp1 <- c(pp1,round(p1,4));
pp2 <- c(pp2,round(p2,4));
pp3 <- c(pp3,round(p3,4));
 }
p.val <- cbind(pp1,pp2,pp3);
rownames(p.val) <- gene;
ind <- apply(p.val,1,function(x)any(x<0.05));
index <- which(pp1<0.05);
dge <- gene[index];
p.val[index,];


library(pheatmap);
exp.gl <- apply(exp.g[index,],1,function(x)log2(x+0.05))
exp.gl <- t(exp.gl);

fac.grp <- c(rep("normal",length(nor)),rep("low.risk",length(low.r)),rep("high.risk",length(high.r)))
ann <- data.frame(group=factor(fac.grp));
rownames(ann) <- c(nor,low.r,high.r);
cols <- colorRampPalette(c("blue","white","red" ))(100);
ann_colors = list(
               group = c(high.risk="#619CFF",low.risk="#F8766D",normal="#00BA38")
                 )
pheatmap(exp.gl[,c(nor,high.r,low.r)], cluster_cols =F,annotation_col=ann,scale="row",color=cols,show_colnames =F,annotation_colors = ann_colors);

