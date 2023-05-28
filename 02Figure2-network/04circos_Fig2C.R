rm(list = ls())
library(circlize)
setwd("D:/HFU/RNA_editing_combine_regulation/res/")

data <- read.table("RADAR_CO_network.txt",header=T,as.is=T)
df1 <- strsplit(data[,1],";")
df1 <- matrix(unlist(df1),ncol=2,byrow=T)
df1 <- as.data.frame(df1)
df1[,2] <- as.numeric(df1[,2]) 
df1 <- df1[,c(1:2,2)] #chr, start end
colnames(df1) <- c("Chr","start","end")
df1 <- data.frame(df1,value=data[,3])


df2 <- strsplit(data[,2],";")
df2 <- matrix(unlist(df2),ncol=2,byrow=T)
df2 <- as.data.frame(df2)
df2[,2] <- as.numeric(df2[,2]) 
df2 <- df2[,c(1:2,2)]
colnames(df2) <- c("Chr","start","end")
df2 <- data.frame(df2,value=data[,3])

df <- rbind(df1, df2)


circos.clear()
# Hg38 human genome
pdf("circos_20230527.pdf",width=10, height=10) # pdf
circos.initializeWithIdeogram(species="hg38")
circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.genomicTrack(df, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, col = rand_color(1), cex = 0.5, ...)
})

circos.genomicLink(df1, df2,  # bed files,chr,start,end,value
                   col= rgb(229, 229, 229, 20, maxColorValue=255), lwd=1,lty=3 )
dev.off()

