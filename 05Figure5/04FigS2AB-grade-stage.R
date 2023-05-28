rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/HFU/RNA_editing_combine_regulation/res/")
getwd()
library(ggplot2)
library(dplyr)
library(ggalluvial)
data3<-read.table("annotation_Sample_TIMER.txt",header=T,sep="\t")
data.g <- data3[!is.na(data3[,"grade"]),];
syn.t12 <- which(data.g[,"group"]=="high.risk" & (data.g[,"grade"]=="G1" | data.g[,"grade"]=="G2"));
syn.t34 <- which(data.g[,"group"]=="high.risk" & (data.g[,"grade"]=="G3" | data.g[,"grade"]=="G4"));
Nosyn.t12 <- which(data.g[,"group"]=="Low.risk" & (data.g[,"grade"]=="G1" | data.g[,"grade"]=="G2"));
Nosyn.t34 <- which(data.g[,"group"]=="Low.risk" & (data.g[,"grade"]=="G3" | data.g[,"grade"]=="G4"));
m <- matrix(c(length(syn.t12),length(syn.t34),length(Nosyn.t12),length(Nosyn.t34)),ncol=2,byrow=T);
fisher.test(m)


data.st <- data3[!is.na(data3[,"stage.T"]),]
syn.t12 <- which(data.st[,"group"]=="high.risk" & (data.st[,"stage.T"]=="T1" | data.st[,"stage.T"]=="T2"));
syn.t34 <- which(data.st[,"group"]=="high.risk" & (data.st[,"stage.T"]=="T3" | data.st[,"stage.T"]=="T4"));
Nosyn.t12 <- which(data.st[,"group"]=="Low.risk" & (data.st[,"stage.T"]=="T1" | data.st[,"stage.T"]=="T2"));
Nosyn.t34 <- which(data.st[,"group"]=="Low.risk" & (data.st[,"stage.T"]=="T3" | data.st[,"stage.T"]=="T4"));
m <- matrix(c(length(syn.t12),length(syn.t34),length(Nosyn.t12),length(Nosyn.t34)),ncol=2,byrow=T);
fisher.test(m)


library(ggplot2)
#group

ggplot( data.g, aes( x = group,  fill = grade))+
  geom_bar( position = "fill")+
  labs(x="", y = "Percentage") +
  scale_fill_manual(values=c("G1"="#FDC6D5","G2" = "#E1B1D5","G3" = "#A19BCB","G4" = "#6D749B"))

##stage
ggplot( data.st, aes( x = group,  fill = stage.T))+
  geom_bar( position = "fill")+
  labs(x="", y = "Percentage") +
  scale_fill_manual(values=c("T1"="#a4cbec","T2" = "#69c2e3","T3" = "#7da2ce","T4" = "#1383c2"))

