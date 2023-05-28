###Plot Fig7D
##data saved from 1-5
colorss <- c("red","orange","blue","darkgreen","purple","#CCCCFF","#9999FF");
plot(perf1,col=colorss[1],type="l",lwd=2)
plot(perf2,col=colorss[2],type="l",lwd=2,add=T)
plot(perf3,col=colorss[3],type="l",lwd=2,add=T)
plot(perf4,col=colorss[4],type="l",lwd=2,add=T)
plot(perf5,col=colorss[5],type="l",lwd=2,add=T);
plot(perf5.1,col=colorss[6],type="l",lwd=2,add=T);
plot(perf5.2,col=colorss[7],type="l",lwd=2,add=T);
legend(0.7,0.3,c("TCGA","GSE65485","GSE169289","GSE77314","GSE164359-total","GSE164359-primary","GSE164359-recurent"),col=colorss,lty=1,text.col=colorss,lwd=2);
text(0.7,0.8,paste0("TCGA: n=",N1," AUC=",p1),col=colorss[1]);
text(0.7,0.75,paste0("GSE65485: n=",N2," AUC=",p2),col=colorss[2]);
text(0.7,0.7,paste0("GSE169289: n=",N3," AUC=",p3),col=colorss[3]);
text(0.7,0.65,paste0("GSE77314: n=",N4," AUC=",p4),col=colorss[4]);
text(0.7,0.6,paste0("GSE164359-total: n=",N5," AUC=",p5),col=colorss[5]);
text(0.7,0.55,paste0("GSE164359-primary: n=",N5.1," AUC=",p5.1),col=colorss[6]);
text(0.7,0.5,paste0("GSE164359-recurent: n=",N5.2," AUC=",p5.2),col=colorss[7]);
f=function(x){ y=x;return(y)}
curve(f(x),0,1,col="grey",lwd=2,lty=2,add=T)

