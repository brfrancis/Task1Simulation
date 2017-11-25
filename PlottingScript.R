setwd("C:\\Users/Ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/")
setwd("/home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/")

png("LRT_Timing.png",width=1200,height=900)
par(mar=c(5.1,5.1,4.1,2.1))
plot(seq(1,6,1),c(14/60,27/60,85/60,290/60,1220/60,123), ylim=c(0,150),xlim=c(0,7),ylab="Computation Time (minutes)",xlab="Number of Samples",xaxt="n",cex=3.5,cex.axis=2.0,cex.lab=2.0)
lines(seq(1,6,1),c(14/60,27/60,85/60,290/60,1220/60,123))
axis(1, at=1:6, labels=c(100,250,500,1000,2000,5000),cex=2.0,cex.axis=2.0,cex.lab=2.0)
dev.off()
#system("scp -P 8022 gaia:/work/projects/epipgx/users/LIV/Task1Surv/Simulation/temp/Results*_1000_0.4.txt /home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/")

datapos5pos5=read.csv("Results0.5_0.5_1000_0.4.txt",header=FALSE,sep=" ")
datapos5pos25=read.csv("Results0.5_0.25_1000_0.4.txt",header=FALSE, sep=" ")
datapos50=read.csv("Results0.5_0_1000_0.4.txt",header=FALSE, sep=" ")
datapos5neg25=read.csv("Results0.5_-0.25_1000_0.4.txt",header=FALSE, sep=" ")
datapos5neg5=read.csv("Results0.5_-0.5_1000_0.4.txt",header=FALSE, sep=" ")

datapos25pos5=read.csv("Results0.25_0.5_1000_0.4.txt",header=FALSE, sep=" ")
datapos25pos25=read.csv("Results0.25_0.25_1000_0.4.txt",header=FALSE, sep=" ")
datapos250=read.csv("Results0.25_0_1000_0.4.txt",header=FALSE, sep=" ")
datapos25neg25=read.csv("Results0.25_-0.25_1000_0.4.txt",header=FALSE, sep=" ")
datapos25neg5=read.csv("Results0.25_-0.5_1000_0.4.txt",header=FALSE, sep=" ")

data0pos5=read.csv("Results0_0.5_1000_0.4.txt",header=FALSE, sep=" ")
data0pos25=read.csv("Results0_0.25_1000_0.4.txt",header=FALSE, sep=" ")
data00=read.csv("Results0_0_1000_0.4.txt",header=FALSE, sep=" ")
data0neg25=read.csv("Results0_-0.25_1000_0.4.txt",header=FALSE, sep=" ")
data0neg5=read.csv("Results0_-0.5_1000_0.4.txt",header=FALSE, sep=" ")

dataneg25pos5=read.csv("Results-0.25_0.5_1000_0.4.txt",header=FALSE, sep=" ")
dataneg25pos25=read.csv("Results-0.25_0.25_1000_0.4.txt",header=FALSE, sep=" ")
dataneg250=read.csv("Results-0.25_0_1000_0.4.txt",header=FALSE, sep=" ")
dataneg25neg25=read.csv("Results-0.25_-0.25_1000_0.4.txt",header=FALSE, sep=" ")
dataneg25neg5=read.csv("Results-0.25_-0.5_1000_0.4.txt",header=FALSE, sep=" ")

dataneg5pos5=read.csv("Results-0.5_0.5_1000_0.4.txt",header=FALSE, sep=" ")
dataneg5pos25=read.csv("Results-0.5_0.25_1000_0.4.txt",header=FALSE, sep=" ")
dataneg50=read.csv("Results-0.5_0_1000_0.4.txt",header=FALSE, sep=" ")
dataneg5neg25=read.csv("Results-0.5_-0.25_1000_0.4.txt",header=FALSE, sep=" ")
dataneg5neg5=read.csv("Results-0.5_-0.5_1000_0.4.txt",header=FALSE, sep=" ")

png("LRT_Prop_Sig_1000_0_4_MAF.png",width=1200,height=900)
par(mfrow=c(2,3),mar=c(5.1,5.1,4.1,2.1))
plot(jitter(rep(1,6),factor=2.5),(colSums(datapos5neg5[seq(3,dim(datapos5neg5)[1],3),]<=0.05)/(dim(datapos5neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,6),col=c("blue","red","green","green"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.5",xaxt="n",pch=c(1,1,1,2,2,2),cex=1.5)
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,6),factor=2.5),(colSums(datapos5neg25[seq(3,dim(datapos5neg25)[1],3),]<=0.05)/(dim(datapos5neg25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(3,6),factor=2.5),(colSums(datapos50[seq(3,dim(datapos50)[1],3),]<=0.05)/(dim(datapos50)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(4,6),factor=2.5),(colSums(datapos5pos25[seq(3,dim(datapos5pos25)[1],3),]<=0.05)/(dim(datapos5pos25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(5,6),factor=2.5),(colSums(datapos5pos5[seq(3,dim(datapos5pos5)[1],3),]<=0.05)/(dim(datapos5pos5)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)


plot(jitter(rep(1,6),factor=2.5),(colSums(datapos25neg5[seq(3,dim(datapos25neg5)[1],3),]<=0.05)/(dim(datapos25neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,6),col=c("blue","red","green","green"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.25",xaxt="n",pch=c(1,1,1,2,2,2),cex=1.5)
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,6),factor=2.5),(colSums(datapos25neg25[seq(3,dim(datapos25neg25)[1],3),]<=0.05)/(dim(datapos25neg25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(3,6),factor=2.5),(colSums(datapos250[seq(3,dim(datapos250)[1],3),]<=0.05)/(dim(datapos250)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(4,6),factor=2.5),(colSums(datapos25pos25[seq(3,dim(datapos25pos25)[1],3),]<=0.05)/(dim(datapos25pos25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(5,6),factor=2.5),(colSums(datapos25pos5[seq(3,dim(datapos25pos5)[1],3),]<=0.05)/(dim(datapos25pos5)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)

plot(jitter(rep(1,6),factor=2.5),(colSums(data0neg5[seq(3,dim(data0neg5)[1],3),]<=0.05)/(dim(data0neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,6),col=c("blue","red","green","green"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0",xaxt="n",pch=c(1,1,1,2,2,2),cex=1.5)
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,6),factor=2.5),(colSums(data0neg25[seq(3,dim(data0neg25)[1],3),]<=0.05)/(dim(data0neg25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(3,6),factor=2.5),(colSums(data00[seq(3,dim(data00)[1],3),]<=0.05)/(dim(data00)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(4,6),factor=2.5),(colSums(data0pos25[seq(3,dim(data0pos25)[1],3),]<=0.05)/(dim(data0pos25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(5,6),factor=2.5),(colSums(data0pos5[seq(3,dim(data0pos5)[1],3),]<=0.05)/(dim(data0pos5)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)

plot(jitter(rep(1,6),factor=2.5),(colSums(dataneg25neg5[seq(3,dim(dataneg25neg5)[1],3),]<=0.05)/(dim(dataneg25neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,6),col=c("blue","red","green","green"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.25",xaxt="n",pch=c(1,1,1,2,2,2),cex=1.5)
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,6),factor=2.5),(colSums(dataneg25neg25[seq(3,dim(dataneg25neg25)[1],3),]<=0.05)/(dim(dataneg25neg25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(3,6),factor=2.5),(colSums(dataneg250[seq(3,dim(dataneg250)[1],3),]<=0.05)/(dim(dataneg250)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(4,6),factor=2.5),(colSums(dataneg25pos25[seq(3,dim(dataneg25pos25)[1],3),]<=0.05)/(dim(dataneg25pos25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(5,6),factor=2.5),(colSums(dataneg25pos5[seq(3,dim(dataneg25pos5)[1],3),]<=0.05)/(dim(dataneg25pos5)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)

plot(jitter(rep(1,6),factor=2.5),(colSums(dataneg5pos5[seq(3,dim(dataneg5neg5)[1],3),]<=0.05)/(dim(dataneg5neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,6),col=c("blue","red","green","green"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.5",xaxt="n",pch=c(1,1,1,2,2,2),cex=1.5)
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,6),factor=2.5),(colSums(dataneg5neg25[seq(3,dim(dataneg5neg25)[1],3),]<=0.05)/(dim(dataneg5neg25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(3,6),factor=2.5),(colSums(dataneg50[seq(3,dim(dataneg50)[1],3),]<=0.05)/(dim(dataneg50)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(4,6),factor=2.5),(colSums(dataneg5pos25[seq(3,dim(dataneg5pos25)[1],3),]<=0.05)/(dim(dataneg5pos25)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)
points(jitter(rep(5,6),factor=2.5),(colSums(dataneg5pos5[seq(3,dim(dataneg5pos5)[1],3),]<=0.05)/(dim(dataneg5pos5)[1]/3))[c(-4,-5)],col=c("blue","red","green","green"),pch=c(1,1,1,2,2,2),cex=1.5)

plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft",legend=c("Lu Survival","Residual Survival","Lu Susceptibility","Residual Susceptibility","Lu Pleiotropy","Residual Pleiotropy"),pch=c(1,2,1,2,1,2), pt.cex=2,cex=2,col=c("blue","blue","red","red","green","green"))

dev.off()


#system("scp -P 8022 gaia:/work/projects/epipgx/users/LIV/Task1Surv/Simulation/temp/Results*_2000_0.4.txt /home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/")

datapos5pos5=read.csv("Results0.5_0.5_2000_0.4.txt",header=FALSE,sep=" ")
datapos5pos25=read.csv("Results0.5_0.25_2000_0.4.txt",header=FALSE, sep=" ")
datapos50=read.csv("Results0.5_0_2000_0.4.txt",header=FALSE, sep=" ")
datapos5neg25=read.csv("Results0.5_-0.25_2000_0.4.txt",header=FALSE, sep=" ")
datapos5neg5=read.csv("Results0.5_-0.5_2000_0.4.txt",header=FALSE, sep=" ")

datapos25pos5=read.csv("Results0.25_0.5_2000_0.4.txt",header=FALSE, sep=" ")
datapos25pos25=read.csv("Results0.25_0.25_2000_0.4.txt",header=FALSE, sep=" ")
datapos250=read.csv("Results0.25_0_2000_0.4.txt",header=FALSE, sep=" ")
datapos25neg25=read.csv("Results0.25_-0.25_2000_0.4.txt",header=FALSE, sep=" ")
datapos25neg5=read.csv("Results0.25_-0.5_2000_0.4.txt",header=FALSE, sep=" ")

data0pos5=read.csv("Results0_0.5_2000_0.4.txt",header=FALSE, sep=" ")
data0pos25=read.csv("Results0_0.25_2000_0.4.txt",header=FALSE, sep=" ")
data00=read.csv("Results0_0_2000_0.4.txt",header=FALSE, sep=" ")
data0neg25=read.csv("Results0_-0.25_2000_0.4.txt",header=FALSE, sep=" ")
data0neg5=read.csv("Results0_-0.5_2000_0.4.txt",header=FALSE, sep=" ")

dataneg25pos5=read.csv("Results-0.25_0.5_2000_0.4.txt",header=FALSE, sep=" ")
dataneg25pos25=read.csv("Results-0.25_0.25_2000_0.4.txt",header=FALSE, sep=" ")
dataneg250=read.csv("Results-0.25_0_2000_0.4.txt",header=FALSE, sep=" ")
dataneg25neg25=read.csv("Results-0.25_-0.25_2000_0.4.txt",header=FALSE, sep=" ")
dataneg25neg5=read.csv("Results-0.25_-0.5_2000_0.4.txt",header=FALSE, sep=" ")

dataneg5pos5=read.csv("Results-0.5_0.5_2000_0.4.txt",header=FALSE, sep=" ")
dataneg5pos25=read.csv("Results-0.5_0.25_2000_0.4.txt",header=FALSE, sep=" ")
dataneg50=read.csv("Results-0.5_0_2000_0.4.txt",header=FALSE, sep=" ")
dataneg5neg25=read.csv("Results-0.5_-0.25_2000_0.4.txt",header=FALSE, sep=" ")
dataneg5neg5=read.csv("Results-0.5_-0.5_2000_0.4.txt",header=FALSE, sep=" ")

png("LRT_Prop_Sig_2000_0_4_NullEffects.png",width=1200,height=900)
par(mfrow=c(1,2),mar=c(5.1,5.1,4.1,2.1))
plot(jitter(rep(1,4),factor=2.5),(colSums(dataneg50[seq(3,dim(dataneg50)[1],3),]<=0.05)/(dim(dataneg50)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,6),col=c("blue","red","darkgreen","darkgreen"),ylab="Prop LRT p<0.05",xlab="Survival Effect",main = "Null Susceptibility Effect",xaxt="n",pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0,cex.main=2.5)
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5),cex.axis=2.0)
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,4),factor=2.5),(colSums(dataneg250[seq(3,dim(dataneg250)[1],3),]<=0.05)/(dim(dataneg250)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
points(jitter(rep(3,4),factor=2.5),(colSums(data00[seq(3,dim(data00)[1],3),]<=0.05)/(dim(data00)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
points(jitter(rep(4,4),factor=2.5),(colSums(datapos250[seq(3,dim(datapos250)[1],3),]<=0.05)/(dim(datapos250)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
points(jitter(rep(5,4),factor=2.5),(colSums(datapos50[seq(3,dim(datapos50)[1],3),]<=0.05)/(dim(datapos50)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
#legend("center",legend=c("Lu Survival","Lu Susceptibility","Lu Pleiotropy","SCOPA"),pch=c(19,19,19,17), pt.cex=2,cex=2,col=c("blue","red","darkgreen","darkgreen"))


plot(jitter(rep(1,4),factor=2.5),(colSums(data0neg5[seq(3,dim(data0neg5)[1],3),]<=0.05)/(dim(data0neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,6),col=c("blue","red","darkgreen","darkgreen"),xlab="Susceptibility Effect",ylab="",main = "Null Survival Effect",xaxt="n",pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0,cex.main=2.5)
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5),cex.axis=2.0)
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,4),factor=2.5),(colSums(data0neg25[seq(3,dim(data0neg25)[1],3),]<=0.05)/(dim(data0neg25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
points(jitter(rep(3,4),factor=2.5),(colSums(data00[seq(3,dim(data00)[1],3),]<=0.05)/(dim(data00)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
points(jitter(rep(4,4),factor=2.5),(colSums(data0pos25[seq(3,dim(data0pos25)[1],3),]<=0.05)/(dim(data0pos25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
points(jitter(rep(5,4),factor=2.5),(colSums(data0pos5[seq(3,dim(data0pos5)[1],3),]<=0.05)/(dim(data0pos5)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.0,cex.axis=2.0,cex.lab=2.0)
dev.off()

png("LRT_Prop_Sig_2000_0_4_MAF_5effects.png",width=1200,height=900)
par(mfrow=c(1,2),mar=c(5.1,5.1,5.1,2.1))
plot(jitter(rep(1,4),factor=2.5),(colSums(datapos5neg5[seq(3,dim(datapos5neg5)[1],3),]<=0.05)/(dim(datapos5neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,5),col=c("blue","red","darkgreen","darkgreen"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.5",xaxt="n",pch=c(19,19,19,17),cex=2.5,cex.axis=2.0,cex.lab=2.0,cex.main=2.5)
axis(1, at=1:4, labels=c(-0.5,-0.25,0.25,0.5),cex.axis=2.0)
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,4),factor=2.5),(colSums(datapos5neg25[seq(3,dim(datapos5neg25)[1],3),]<=0.05)/(dim(datapos5neg25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(3,4),factor=2.5),(colSums(datapos5pos25[seq(3,dim(datapos5pos25)[1],3),]<=0.05)/(dim(datapos5pos25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(4,4),factor=2.5),(colSums(datapos5pos5[seq(3,dim(datapos5pos5)[1],3),]<=0.05)/(dim(datapos5pos5)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)

plot(jitter(rep(1,4),factor=2.5),(colSums(dataneg5pos5[seq(3,dim(dataneg5neg5)[1],3),]<=0.05)/(dim(dataneg5neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,5),col=c("blue","red","darkgreen","darkgreen"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.5",xaxt="n",pch=c(19,19,19,17),cex=2.5,cex.axis=2.0,cex.lab=2.0,cex.main=2.5)
axis(1, at=1:4, labels=c(-0.5,-0.25,0.25,0.5),cex.axis=2.0)
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,4),factor=2.5),(colSums(dataneg5neg25[seq(3,dim(dataneg5neg25)[1],3),]<=0.05)/(dim(dataneg5neg25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(3,4),factor=2.5),(colSums(dataneg5pos25[seq(3,dim(dataneg5pos25)[1],3),]<=0.05)/(dim(dataneg5pos25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(4,4),factor=2.5),(colSums(dataneg5pos5[seq(3,dim(dataneg5pos5)[1],3),]<=0.05)/(dim(dataneg5pos5)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
dev.off()

png("LRT_Prop_Sig_2000_0_4_MAF_25effects.png",width=1200,height=900)
par(mfrow=c(1,2),mar=c(5.1,5.1,5.1,2.1))
plot(jitter(rep(1,4),factor=2.5),(colSums(datapos25neg5[seq(3,dim(datapos25neg5)[1],3),]<=0.05)/(dim(datapos25neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,5),col=c("blue","red","darkgreen","darkgreen"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.25",xaxt="n",pch=c(19,19,19,17),cex=2.5,cex.axis=2.0,cex.lab=2.0,cex.main=2.5)
axis(1, at=1:4, labels=c(-0.5,-0.25,0.25,0.5),cex.axis=2.0)
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,4),factor=2.5),(colSums(datapos25neg25[seq(3,dim(datapos25neg25)[1],3),]<=0.05)/(dim(datapos25neg25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(3,4),factor=2.5),(colSums(datapos25pos25[seq(3,dim(datapos25pos25)[1],3),]<=0.05)/(dim(datapos25pos25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(4,4),factor=2.5),(colSums(datapos25pos5[seq(3,dim(datapos25pos5)[1],3),]<=0.05)/(dim(datapos25pos5)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)

plot(jitter(rep(1,4),factor=2.5),(colSums(dataneg25neg5[seq(3,dim(dataneg25neg5)[1],3),]<=0.05)/(dim(dataneg25neg5)[1]/3))[c(-4,-5)], ylim=c(0,1),xlim=c(0,5),col=c("blue","red","darkgreen","darkgreen"),ylab="Prop LRT p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.25",xaxt="n",pch=c(19,19,19,17),cex=2.5,cex.axis=2.0,cex.lab=2.0,cex.main=2.5)
axis(1, at=1:4, labels=c(-0.5,-0.25,0.25,0.5),cex.axis=2.0)
abline(h=0.8,lty=2)
abline(h=0.05,lty=2)
points(jitter(rep(2,4),factor=2.5),(colSums(dataneg25neg25[seq(3,dim(dataneg25neg25)[1],3),]<=0.05)/(dim(dataneg25neg25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(3,4),factor=2.5),(colSums(dataneg25pos25[seq(3,dim(dataneg25pos25)[1],3),]<=0.05)/(dim(dataneg25pos25)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
points(jitter(rep(4,4),factor=2.5),(colSums(dataneg25pos5[seq(3,dim(dataneg25pos5)[1],3),]<=0.05)/(dim(dataneg25pos5)[1]/3))[c(-4,-5)],col=c("blue","red","darkgreen","darkgreen"),pch=c(19,19,19,17),cex=2.5)
dev.off()

png("LRT_Legend.png",width=900,height=900)
plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft",legend=c("Lu Survival","Lu Susceptibility","Lu Surv and Sus","SCOPA"),pch=c(19,19,19,17), pt.cex=2,cex=2,col=c("blue","red","darkgreen","darkgreen"))
dev.off()

# data100=read.csv("PvalueResults_100_0.4.txt",header=FALSE)
# data250=read.csv("PvalueResults_250_0.4.txt",header=FALSE)
# data500=read.csv("PvalueResults_500_0.4.txt",header=FALSE)
# data1000=read.csv("PvalueResults_1000_0.4.txt",header=FALSE)
# data2000=read.csv("PvalueResults_2000_0.4.txt",header=FALSE)
# 
# #Convert Function
# datatrans=function(dat){
#   Null1a=c((dat[1:6,2]+dat[1:6,3])/1000,(dat[1:6,6]+dat[1:6,7])/1000)
#   alp1=Null1a[c(1,2,7,3,4,5,11,6,12)]
#   Null1p=jitter(rep(0,12))#c((dat[1:6,2]+dat[1:6,3])/1000,(dat[1:6,6]+dat[1:6,7])/1000)
#   pow1=Null1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Lu1a=c((dat[7:12,2]+dat[7:12,3])/1000,(dat[7:12,6]+dat[7:12,7])/1000)
#   alp2=Lu1a[c(1,2,7,3,4,5,11,6,12)]
#   Lu1p=jitter(rep(0,12))#c((dat[7:12,2]+dat[7:12,3])/1000,(dat[7:12,6]+dat[7:12,7])/1000)
#   pow2=Lu1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Surv1a=c((dat[13:18,3])/500,(dat[13:18,6]+dat[13:18,7])/1000)
#   alp3=Surv1a[c(1,2,7,3,4,5,11,6,12)]
#   Surv1p=c((dat[13:18,2])/500,jitter(rep(0,6)))
#   pow3=Surv1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Sus1a=c((dat[19:24,2]+dat[19:24,3])/1000,(dat[19:24,7])/500)
#   alp4=Sus1a[c(1,2,7,3,4,5,11,6,12)]
#   Sus1p=c(jitter(rep(0,6),(dat[19:24,6])/500))
#   pow4=Sus1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Both1a=c((dat[25:30,3])/500,(dat[25:30,7])/500)
#   alp5=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c((dat[25:30,2])/500,(dat[25:30,6])/500)
#   pow5=Both1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   alp=c(alp1,alp2,alp3,alp4,alp5)
#   pow=c(pow1,pow2,pow3,pow4,pow5)
#   return(cbind(alp,pow))
# }
# 
# datar100=datatrans(data100)
# datar250=datatrans(data250)
# datar500=datatrans(data500)
# datar1000=datatrans(data1000)
# datar2000=datatrans(data2000)
# 
# 
# 
# 
# png("Sig_0_4_MAF.png",width=1200,height=900)
# par(mfrow=c(2,3))
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Sample Size",main = "\nNull model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
# points(jitter(rep(2,9),factor=2.5),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(5,9),factor=2.5),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# abline(h=0.05,lty=2)
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# st=10
# fin=18
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Sample Size",main = "Significance with 0.4 MAF\nClinical Only model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
# points(jitter(rep(2,9),factor=2.5),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(5,9),factor=2.5),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# abline(h=0.05,lty=2)
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# st=19
# fin=27
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Sample Size",main = "\nGenetic Survival model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
# points(jitter(rep(2,9),factor=2.5),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(5,9),factor=2.5),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# abline(h=0.05,lty=2)
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# st=28
# fin=36
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Sample Size",main = "\nGenetic Susceptibility model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
# points(jitter(rep(2,9),factor=2.5),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(5,9),factor=2.5),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# abline(h=0.05,lty=2)
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# st=37
# fin=45
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Sample Size",main = "\nGenetic Pleiotropy model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
# points(jitter(rep(2,9),factor=2.5),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# #points(jitter(rep(5,9),factor=2.5),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# abline(h=0.05,lty=2)
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# dev.off()
# 
# 
# png("Power_0_4_MAF.png",width=1600,height=450)
# par(mfrow=c(1,4))
# st=19
# fin=27
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,1),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Sample Size",main = "\nGenetic Survival model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:5, labels=c(100,250,500,1000,2000))
# points(jitter(rep(2,9),factor=2.5),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(5,9),factor=2.5),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# abline(h=0.8,lty=2)
# 
# st=28
# fin=36
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,1),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Sample Size",main = "Power with 0.4 MAF\nGenetic Susceptibility model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:5, labels=c(100,250,500,1000,2000))
# points(jitter(rep(2,9),factor=2.5),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(5,9),factor=2.5),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# abline(h=0.8,lty=2)
# 
# st=37
# fin=45
# plot(jitter(rep(1,4),factor=2.5),factor=2.5),datar100[st:fin,2], ylim=c(0,1),xlim=c(0,8),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Sample Size",main = "\nGenetic Pleiotropy model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
# axis(1, at=1:5, labels=c(100,250,500,1000,2000))
# points(jitter(rep(2,9),factor=2.5),datar250[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(3,9),factor=2.5),datar500[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# points(jitter(rep(4,9),factor=2.5),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# #points(jitter(rep(5,9),factor=2.5),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3))
# #legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# abline(h=0.8,lty=2)
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# dev.off()
# 
# 
# ####Vary parameters
# 
# datapos5pos5=read.csv("PvalueResults0.5_0.5_2000_0.4.txt",header=FALSE)
# datapos5pos25=read.csv("PvalueResults0.5_0.25_2000_0.4.txt",header=FALSE)
# datapos50=read.csv("PvalueResults0.5_0_2000_0.4.txt",header=FALSE)
# datapos5neg25=read.csv("PvalueResults0.5_-0.25_2000_0.4.txt",header=FALSE)
# datapos5neg5=read.csv("PvalueResults0.5_-0.5_2000_0.4.txt",header=FALSE)
# 
# datapos25pos5=read.csv("PvalueResults0.25_0.5_2000_0.4.txt",header=FALSE)
# datapos25pos25=read.csv("PvalueResults0.25_0.25_2000_0.4.txt",header=FALSE)
# datapos250=read.csv("PvalueResults0.25_0_2000_0.4.txt",header=FALSE)
# datapos25neg25=read.csv("PvalueResults0.25_-0.25_2000_0.4.txt",header=FALSE)
# datapos25neg5=read.csv("PvalueResults0.25_-0.5_2000_0.4.txt",header=FALSE)
# 
# data0pos5=read.csv("PvalueResults0_0.5_2000_0.4.txt",header=FALSE)
# data0pos25=read.csv("PvalueResults0_0.25_2000_0.4.txt",header=FALSE)
# data00=read.csv("PvalueResults0_0_2000_0.4.txt",header=FALSE)
# data0neg25=read.csv("PvalueResults0_-0.25_2000_0.4.txt",header=FALSE)
# data0neg5=read.csv("PvalueResults0_-0.5_2000_0.4.txt",header=FALSE)
# 
# dataneg25pos5=read.csv("PvalueResults-0.25_0.5_2000_0.4.txt",header=FALSE)
# dataneg25pos25=read.csv("PvalueResults-0.25_0.25_2000_0.4.txt",header=FALSE)
# dataneg250=read.csv("PvalueResults-0.25_0_2000_0.4.txt",header=FALSE)
# dataneg25neg25=read.csv("PvalueResults-0.25_-0.25_2000_0.4.txt",header=FALSE)
# dataneg25neg5=read.csv("PvalueResults-0.25_-0.5_2000_0.4.txt",header=FALSE)
# 
# dataneg5pos5=read.csv("PvalueResults-0.5_0.5_2000_0.4.txt",header=FALSE)
# dataneg5pos25=read.csv("PvalueResults-0.5_0.25_2000_0.4.txt",header=FALSE)
# dataneg50=read.csv("PvalueResults-0.5_0_2000_0.4.txt",header=FALSE)
# dataneg5neg25=read.csv("PvalueResults-0.5_-0.25_2000_0.4.txt",header=FALSE)
# dataneg5neg5=read.csv("PvalueResults-0.5_-0.5_2000_0.4.txt",header=FALSE)
# 
# 
# #Convert Function
# datatrans1<-function(dat){
#   Both1a=c((dat[1:6,3])/500,(dat[1:6,7])/500)
#   alp1=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c((dat[1:6,2])/500,(dat[1:6,6])/500)
#   pow1=Both1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Both1a=c((dat[7:12,3])/500,(dat[7:12,7])/500)
#   alp2=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c((dat[7:12,2])/500,(dat[7:12,6])/500)
#   pow2=Both1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Surv1a=c((dat[13:18,3])/500,(dat[13:18,6]+dat[13:18,7])/1000)
#   alp3=Surv1a[c(1,2,7,3,4,5,11,6,12)]
#   Surv1p=c((dat[13:18,2])/500,rep(0,6))
#   pow3=Surv1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Both1a=c((dat[19:24,3])/500,(dat[19:24,7])/500)
#   alp4=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c((dat[19:24,2])/500,(dat[19:24,6])/500)
#   pow4=Both1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   Both1a=c((dat[25:30,3])/500,(dat[25:30,7])/500)
#   alp5=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c((dat[25:30,2])/500,(dat[25:30,6])/500)
#   pow5=Both1p[c(1,2,7,3,4,5,11,6,12)]
#   
#   alp=c(alp1,alp2,alp3,alp4,alp5)
#   pow=c(pow1,pow2,pow3,pow4,pow5)
#   return(cbind(alp,pow))
# }
# 
# datatrans2<-function(dat){
#   Both1a=c((dat[1:6,3]+dat[1:6,2])/1000,(dat[1:6,7])/500)
#   alp1=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c(rep(0,6),(dat[1:6,6])/500)
#   pow1=Both1p[c(1,2,7,9,10,5,11,6,12)]
#   
#   Both1a=c((dat[7:12,3]+dat[7:12,2])/1000,(dat[7:12,7])/500)
#   alp2=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c(rep(0,6),(dat[7:12,6])/500)
#   pow2=Both1p[c(1,2,7,9,10,5,11,6,12)]
#   
#   Surv1a=c((dat[13:18,3]+dat[13:18,2])/1000,(dat[13:18,6]+dat[13:18,7])/1000)
#   alp3=Surv1a[c(1,2,7,3,4,5,11,6,12)]
#   Surv1p=c(rep(0,6),rep(0,6))
#   pow3=Surv1p[c(1,2,7,9,10,5,11,6,12)]
#   
#   Both1a=c((dat[19:24,3]+dat[19:24,2])/1000,(dat[19:24,7])/500)
#   alp4=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c(rep(0,6),(dat[19:24,6])/500)
#   pow4=Both1p[c(1,2,7,9,10,5,11,6,12)]
#   
#   Both1a=c((dat[25:30,3]+dat[25:30,2])/1000,(dat[25:30,7])/500)
#   alp5=Both1a[c(1,2,7,3,4,5,11,6,12)]
#   Both1p=c(rep(0,6),(dat[25:30,6])/500)
#   pow5=Both1p[c(1,2,7,9,10,5,11,6,12)]
#   
#   alp=c(alp1,alp2,alp3,alp4,alp5)
#   pow=c(pow1,pow2,pow3,pow4,pow5)
#   return(cbind(alp,pow))
# }
# 
# datatrans3<-function(dat){
#   Both1a=c((dat[1:6,2])/500,(dat[1:6,6])/500)
#   alp1=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(0,0,0,0,0,0,0,0,0)
#   Both1p=c((dat[1:6,2])/500,(dat[1:6,6])/500)
#   pow1=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(1,1,1,1,1,1,1,1,1)
#   
#   Both1a=c((dat[7:12,2])/500,(dat[7:12,6])/500)
#   alp2=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(0,0,0,0,0,0,0,0,0)
#   Both1p=c((dat[7:12,2])/500,(dat[7:12,6])/500)
#   pow2=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(1,1,1,1,1,1,1,1,1)
#   
#   Surv1a=c((dat[13:18,2])/500,(dat[13:18,6])/500)
#   alp3=Surv1a[c(1,2,7,3,4,5,11,6,12)]*c(0,0,1,0,0,0,1,0,1)
#   Surv1p=c((dat[13:18,2])/500,(dat[13:18,6])/500)
#   pow3=Surv1p[c(1,2,7,3,4,5,11,6,12)]*c(1,1,0,1,1,1,0,1,0)
#   
#   Both1a=c((dat[19:24,2])/500,(dat[19:24,6])/500)
#   alp4=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(0,0,0,0,0,0,0,0,0)
#   Both1p=c((dat[19:24,2])/500,(dat[19:24,6])/500)
#   pow4=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(1,1,1,1,1,1,1,1,1)
#   
#   Both1a=c((dat[25:30,2])/500,(dat[25:30,6])/500)
#   alp5=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(0,0,0,0,0,0,0,0,0)
#   Both1p=c((dat[25:30,2])/500,(dat[25:30,6])/500)
#   pow5=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(1,1,1,1,1,1,1,1,1)
#   
#   alp=c(alp1,alp2,alp3,alp4,alp5)
#   pow=c(pow1,pow2,pow3,pow4,pow5)
#   return(cbind(alp,pow))
# }
# 
# datatrans4<-function(dat){
#   Both1a=c(dat[1:6,2]/500,dat[1:6,6]/500)
#   alp1=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(1,1,0,0,0,1,0,1,0)
#   Both1p=c(dat[1:6,2]/500,dat[1:6,6]/500)
#   pow1=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(0,0,1,1,1,0,1,0,1)
#   
#   Both1a=c(dat[7:12,2]/500,dat[7:12,6]/500)
#   alp2=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(1,1,0,0,0,1,0,1,0)
#   Both1p=c(dat[7:12,2]/500,dat[7:12,6]/500)
#   pow2=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(0,0,1,1,1,0,1,0,1)
#   
#   Surv1a=c(dat[13:18,2]/500,dat[13:18,6]/500)
#   alp3=Surv1a[c(1,2,7,3,4,5,11,6,12)]*c(1,1,1,1,1,1,1,1,1)
#   Surv1p=c(dat[13:18,2]/500,dat[13:18,6]/500)
#   pow3=Surv1p[c(1,2,7,3,4,5,11,6,12)]*c(0,0,0,0,0,0,0,0,0)
#   
#   Both1a=c(dat[19:24,2]/500,dat[19:24,6]/500)
#   alp4=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(1,1,0,0,0,1,0,1,0)
#   Both1p=c(dat[19:24,2]/500,dat[19:24,6]/500)
#   pow4=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(0,0,1,1,1,0,1,0,1)
#   
#   Both1a=c(dat[25:30,2]/500,dat[25:30,6]/500)
#   alp5=Both1a[c(1,2,7,3,4,5,11,6,12)]*c(1,1,0,0,0,1,0,1,0)
#   Both1p=c(dat[25:30,2]/500,(dat[25:30,6])/500)
#   pow5=Both1p[c(1,2,7,3,4,5,11,6,12)]*c(0,0,1,1,1,0,1,0,1)
#   
#   alp=c(alp1,alp2,alp3,alp4,alp5)
#   pow=c(pow1,pow2,pow3,pow4,pow5)
#   return(cbind(alp,pow))
# }
# 
# datatrans5<-function(dat){
#   SNPs3=dat[seq(3,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   FULL3=dat[seq(15,9000,18),]
#   FIXd3=dat[seq(18,9000,18),]
#   alp=c(sum(SNPs3[,2]<0.05)/500,sum(SNPe3[,2]<0.05)/500,sum(SNPs3[,6]<0.05)/500,sum(PLEs3[,2]<0.05)/500,sum(PLEe3[,2]<0.05)/500,sum(FIXd3[,2]<0.05)/500,sum(FIXd3[,6]<0.05)/500,sum(FULL3[,2]<0.05)/500,sum(FULL3[,6]<0.05)/500)
#   return(alp)
# }
# 
# dataneg5=rbind(datatrans5(Rdataneg5neg5),datatrans5(Rdataneg5neg25),datatrans5(Rdataneg50),datatrans5(Rdataneg5pos25),datatrans5(Rdataneg5pos5))[,c(4,8,9)]
# dataneg25=rbind(datatrans5(Rdataneg25neg5),datatrans5(Rdataneg25neg25),datatrans5(Rdataneg250),datatrans5(Rdataneg25pos25),datatrans5(Rdataneg25pos5))[,c(4,8,9)]
# data0=rbind(datatrans5(Rdata0neg5),datatrans5(Rdata0neg25),datatrans5(Rdata00),datatrans5(Rdata0pos25),datatrans5(Rdata0pos5))[,c(4,8,9)]
# datapos25=rbind(datatrans5(Rdatapos25neg5),datatrans5(Rdatapos25neg25),datatrans5(Rdatapos250),datatrans5(Rdatapos25neg25),datatrans5(Rdatapos25pos5))[,c(4,8,9)]
# datapos5=rbind(datatrans5(Rdatapos5neg5),datatrans5(Rdatapos5neg25),datatrans5(Rdatapos50),datatrans5(Rdatapos5pos25),datatrans5(Rdatapos5pos5))[,c(4,8,9)]
# 
# 
# png("VaryPara_Prop_2000_0_4_MAF.png",width=1200,height=900)
# par(mfrow=c(2,3))
# plot(jitter(rep(1,3),factor = 3),dataneg5[1,], ylim = c(0,1),xlim=c(0,6),col=c("red","green","green"),ylab="Prop p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.5",xaxt="n",pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2);abline(h=0.8,lty=2)
# points(jitter(rep(2,3),factor = 2.5),dataneg5[2,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(3,3),factor = 2.5),dataneg5[3,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(4,3),factor = 2.5),dataneg5[4,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(5,3),factor = 2.5),dataneg5[5,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# 
# plot(jitter(rep(1,3),factor=2.5),dataneg25[1,], ylim = c(0,1),xlim=c(0,6),col=c("red","green","green"),ylab="Prop p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.25",xaxt="n",pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2);abline(h=0.8,lty=2)
# points(jitter(rep(2,3),factor = 2.5),dataneg25[2,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(3,3),factor = 2.5),dataneg25[3,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(4,3),factor = 2.5),dataneg25[4,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(5,3),factor = 2.5),dataneg25[5,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# 
# plot(jitter(rep(1,3),factor=2.5),data0[1,], ylim = c(0,1),xlim=c(0,6),col=c("red","green","green"),ylab="Prop p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0",xaxt="n",pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2);abline(h=0.8,lty=2)
# points(jitter(rep(2,3),factor = 2.5),data0[2,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(3,3),factor = 2.5),data0[3,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(4,3),factor = 2.5),data0[4,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(5,3),factor = 2.5),data0[5,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# 
# plot(jitter(rep(1,3),factor=2.5),datapos25[1,], ylim = c(0,1),xlim=c(0,6),col=c("red","green","green"),ylab="Prop p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.25",xaxt="n",pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2);abline(h=0.8,lty=2)
# points(jitter(rep(2,3),factor = 2.5),datapos25[2,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(3,3),factor = 2.5),datapos25[3,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(4,3),factor = 2.5),datapos25[4,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(5,3),factor = 2.5),datapos25[5,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# 
# plot(jitter(rep(1,3),factor=2.5),datapos5[1,], ylim = c(0,1),xlim=c(0,6),col=c("red","green","green"),ylab="Prop p<0.05",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.5",xaxt="n",pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2);abline(h=0.8,lty=2)
# points(jitter(rep(2,3),factor = 2.5),datapos5[2,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(3,3),factor = 2.5),datapos5[3,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(4,3),factor = 2.5),datapos5[4,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# points(jitter(rep(5,3),factor = 2.5),datapos5[5,],col=c("red","green","green"),pch=c(1,1,3),cex=c(1.5,1.5,1.5))
# 
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# legend("topleft",legend=c("Pleio(s)","Full(surv)","Full(sus)"),pch=c(1,2,1,3), pt.cex=2,cex=2,col=c("red","green","green"))
# 
# dev.off()
# 
# png("VaryPara_Sig_2000_0_4_MAF.png",width=1200,height=900)
# par(mfrow=c(2,3))
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),dataneg5[st:fin,1], ylim = c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),dataneg5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),dataneg5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),dataneg5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),dataneg5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),dataneg25[st:fin,1], ylim = c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),dataneg25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),dataneg25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),dataneg25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),dataneg25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),data0[st:fin,1], ylim = c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),data0[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),data0[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),data0[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),data0[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),datapos25[st:fin,1], ylim = c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),datapos25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),datapos25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),datapos25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),datapos25[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),datapos5[st:fin,1], ylim = c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.05,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),datapos5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),datapos5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),datapos5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),datapos5[st:fin,1],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# 
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3), pt.cex=2,cex=2,col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# dev.off()
# 
# 
# png("VaryPara_Power_2000_0_4_MAF.png",width=1200,height=900)
# par(mfrow=c(2,3))
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),dataneg5[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.8,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),dataneg25[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.8,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),data0[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.8,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),data0[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),data0[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),data0[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),data0[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),datapos25[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.8,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# st=1
# fin=9
# plot(jitter(rep(1,4),factor=2.5),datapos5[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
# abline(h=0.8,lty=2)
# st=10
# fin=18
# points(jitter(rep(2,9),factor=2.5),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=19
# fin=27
# points(jitter(rep(3,9),factor=2.5),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=28
# fin=36
# points(jitter(rep(4,9),factor=2.5),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# st=37
# fin=45
# points(jitter(rep(5,9),factor=2.5),factor=2.5),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"),pch=c(1,2,3,1,2,1,3,1,3),cex=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5))
# 
# 
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3), pt.cex=2,cex=2,col=c("blue","blue","blue","red","red","darkgray","darkgray","green","green"))
# 
# dev.off()
# 
# 
# Rdatapos5pos5=read.csv("Results0.5_0.5_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos5pos25=read.csv("Results0.5_0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos50=read.csv("Results0.5_0_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos5neg25=read.csv("Results0.5_-0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos5neg5=read.csv("Results0.5_-0.5_2000_0.4.txt",header=FALSE,sep=" ")
# 
# Rdatapos25pos5=read.csv("Results0.25_0.5_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos25pos25=read.csv("Results0.25_0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos250=read.csv("Results0.25_0_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos25neg25=read.csv("Results0.25_-0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdatapos25neg5=read.csv("Results0.25_-0.5_2000_0.4.txt",header=FALSE,sep=" ")
# 
# Rdata0pos5=read.csv("Results0_0.5_2000_0.4.txt",header=FALSE,sep=" ")
# Rdata0pos25=read.csv("Results0_0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdata00=read.csv("Results0_0_2000_0.4.txt",header=FALSE,sep=" ")
# Rdata0neg25=read.csv("Results0_-0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdata0neg5=read.csv("Results0_-0.5_2000_0.4.txt",header=FALSE,sep=" ")
# 
# Rdataneg25pos5=read.csv("Results-0.25_0.5_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg25pos25=read.csv("Results-0.25_0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg250=read.csv("Results-0.25_0_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg25neg25=read.csv("Results-0.25_-0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg25neg5=read.csv("Results-0.25_-0.5_2000_0.4.txt",header=FALSE,sep=" ")
# 
# Rdataneg5pos5=read.csv("Results-0.5_0.5_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg5pos25=read.csv("Results-0.5_0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg50=read.csv("Results-0.5_0_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg5neg25=read.csv("Results-0.5_-0.25_2000_0.4.txt",header=FALSE,sep=" ")
# Rdataneg5neg5=read.csv("Results-0.5_-0.5_2000_0.4.txt",header=FALSE,sep=" ")
# 
# Rdatatrans<-function(dat){
#   SNPs=dat[seq(1,9000,18),]
#   SNPs2=dat[seq(2,9000,18),]
#   SNPs3=dat[seq(3,9000,18),]
#   
#   SNPe=dat[seq(4,9000,18),]
#   SNPe2=dat[seq(5,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   
#   PLEs=dat[seq(7,9000,18),]
#   PLEs2=dat[seq(8,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   
#   PLEe=dat[seq(10,9000,18),]
#   PLEe2=dat[seq(11,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   
#   FIXd=dat[seq(13,9000,18),]
#   FIXd2=dat[seq(13,9000,18),]
#   FIXd3=dat[seq(13,9000,18),]
#   
#   FULL=dat[seq(16,9000,18),]
#   FULL2=dat[seq(17,9000,18),]
#   FULL3=dat[seq(18,9000,18),]
#   
#   par(mfrow=c(5,3))
#   plot(SNPe[,2],col=ifelse(SNPe3[,1]>0.05,"green","blue"),pch=2,main="Parameter Value",ylab = "SNPtest(e)")
#   plot(SNPe2[,2],col=ifelse(SNPe3[,1]>0.05,"green","blue"),pch=2,main="SE Value")
#   plot(SNPe3[,2],col=ifelse(SNPe3[,1]>0.05,"green","blue"),pch=2,main="P-value",ylim=c(0,1))
#   plot(SNPs[,2],col=ifelse(SNPs3[,1]>0.05,"green","blue"),ylab = "SNPtest(s)")
#   plot(SNPs2[,2],col=ifelse(SNPs3[,1]>0.05,"green","blue"))
#   plot(SNPs3[,2],col=ifelse(SNPs3[,1]>0.05,"green","blue"),ylim=c(0,1))
#   plot(SNPe[,6],col=ifelse(SNPe3[,5]>0.05,"green","blue"),pch=3,ylab = "SNPtest(sus.dev)")
#   plot(SNPe2[,6],col=ifelse(SNPe3[,5]>0.05,"green","blue"),pch=3)
#   plot(SNPe3[,6],col=ifelse(SNPe3[,5]>0.05,"green","blue"),pch=3,ylim=c(0,1))
#   plot(PLEe[,2],col=ifelse(PLEe3[,1]>0.05,"green","red"),pch=2,ylab = "Pleiotropy (e+sus.dev)")
#   plot(PLEe2[,2],col=ifelse(PLEe3[,1]>0.05,"green","red"),pch=2)
#   plot(PLEe3[,2],col=ifelse(PLEe3[,1]>0.05,"green","red"),pch=2,ylim=c(0,1))
#   plot(PLEs[,2],col=ifelse(PLEs3[,1]>0.05,"green","red"),ylab = "Pleiotropy (s+sus.dev)")
#   plot(PLEs2[,2],col=ifelse(PLEs3[,1]>0.05,"green","red"))
#   plot(PLEs3[,2],col=ifelse(PLEs3[,1]>0.05,"green","red"),ylim=c(0,1))
#   print(c(sum(SNPe3[,1]<0.05)/500,sum(SNPs3[,1]<0.05)/500,sum(SNPs3[,5]<0.05)/500,sum(PLEs3[,1]<0.05)/500,sum(PLEe3[,1]<0.05)/500))
# }
# 
# Rdatatrans2<-function(dat){
#   SNPs=dat[seq(1,9000,18),]
#   SNPs2=dat[seq(2,9000,18),]
#   SNPs3=dat[seq(3,9000,18),]
#   
#   SNPe=dat[seq(4,9000,18),]
#   SNPe2=dat[seq(5,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   
#   PLEs=dat[seq(7,9000,18),]
#   PLEs2=dat[seq(8,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   
#   PLEe=dat[seq(10,9000,18),]
#   PLEe2=dat[seq(11,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   
#   FIXd=dat[seq(13,9000,18),]
#   FIXd2=dat[seq(13,9000,18),]
#   FIXd3=dat[seq(13,9000,18),]
#   
#   FULL=dat[seq(16,9000,18),]
#   FULL2=dat[seq(17,9000,18),]
#   FULL3=dat[seq(18,9000,18),]
#   
#   par(mfrow=c(5,3))
#   plot(SNPe[,2],col=ifelse(SNPe3[,2]>0.05,"green","blue"),pch=2,main="Parameter Value",ylab = "SNPtest(e)")
#   plot(SNPe2[,2],col=ifelse(SNPe3[,2]>0.05,"green","blue"),pch=2,main="SE Value")
#   plot(SNPe3[,2],col=ifelse(SNPe3[,2]>0.05,"green","blue"),pch=2,main="P-value",ylim=c(0,1))
#   plot(SNPs[,2],col=ifelse(SNPs3[,2]>0.05,"green","blue"),ylab = "SNPtest(s)")
#   plot(SNPs2[,2],col=ifelse(SNPs3[,2]>0.05,"green","blue"))
#   plot(SNPs3[,2],col=ifelse(SNPs3[,2]>0.05,"green","blue"),ylim=c(0,1))
#   plot(SNPe[,6],col=ifelse(SNPe3[,6]>0.05,"green","blue"),pch=3,ylab = "SNPtest(sus.dev)")
#   plot(SNPe2[,6],col=ifelse(SNPe3[,6]>0.05,"green","blue"),pch=3)
#   plot(SNPe3[,6],col=ifelse(SNPe3[,6]>0.05,"green","blue"),pch=3,ylim=c(0,1))
#   plot(PLEe[,2],col=ifelse(PLEe3[,2]>0.05,"green","red"),pch=2,ylab = "Pleiotropy (e+sus.dev)")
#   plot(PLEe2[,2],col=ifelse(PLEe3[,2]>0.05,"green","red"),pch=2)
#   plot(PLEe3[,2],col=ifelse(PLEe3[,2]>0.05,"green","red"),pch=2,ylim=c(0,1))
#   plot(PLEs[,2],col=ifelse(PLEs3[,2]>0.05,"green","red"),ylab = "Pleiotropy (s+sus.dev)")
#   plot(PLEs2[,2],col=ifelse(PLEs3[,2]>0.05,"green","red"))
#   plot(PLEs3[,2],col=ifelse(PLEs3[,2]>0.05,"green","red"),ylim=c(0,1))
#   print(c(sum(SNPe3[,2]<0.05)/500,sum(SNPs3[,2]<0.05)/500,sum(SNPs3[,6]<0.05)/500,sum(PLEs3[,2]<0.05)/500,sum(PLEe3[,2]<0.05)/500))
# }
# 
# Rdatatrans3<-function(dat){
#   SNPs=dat[seq(1,9000,18),]
#   SNPs2=dat[seq(2,9000,18),]
#   SNPs3=dat[seq(3,9000,18),]
#   
#   SNPe=dat[seq(4,9000,18),]
#   SNPe2=dat[seq(5,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   
#   PLEs=dat[seq(7,9000,18),]
#   PLEs2=dat[seq(8,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   
#   PLEe=dat[seq(10,9000,18),]
#   PLEe2=dat[seq(11,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   
#   FIXd=dat[seq(13,9000,18),]
#   FIXd2=dat[seq(13,9000,18),]
#   FIXd3=dat[seq(13,9000,18),]
#   
#   FULL=dat[seq(16,9000,18),]
#   FULL2=dat[seq(17,9000,18),]
#   FULL3=dat[seq(18,9000,18),]
#   
#   par(mfrow=c(5,3))
#   plot(SNPe[,2],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="Parameter Value",ylab = "SNPtest(e)")
#   plot(SNPe2[,2],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="SE Value")
#   plot(SNPe3[,2],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="P-value",ylim=c(0,1))
#   plot(SNPs[,2],col=ifelse(SNPs3[,3]>0.05,"green","blue"),ylab = "SNPtest(s)")
#   plot(SNPs2[,2],col=ifelse(SNPs3[,3]>0.05,"green","blue"))
#   plot(SNPs3[,2],col=ifelse(SNPs3[,3]>0.05,"green","blue"),ylim=c(0,1))
#   plot(SNPe[,6],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3,ylab = "SNPtest(sus.dev)")
#   plot(SNPe2[,6],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3)
#   plot(SNPe3[,6],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3,ylim=c(0,1))
#   plot(PLEe[,2],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2,ylab = "Pleiotropy (e+sus.dev)")
#   plot(PLEe2[,2],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2)
#   plot(PLEe3[,2],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2,ylim=c(0,1))
#   plot(PLEs[,2],col=ifelse(PLEs3[,3]>0.05,"green","red"),ylab = "Pleiotropy (s+sus.dev)")
#   plot(PLEs2[,2],col=ifelse(PLEs3[,3]>0.05,"green","red"))
#   plot(PLEs3[,2],col=ifelse(PLEs3[,3]>0.05,"green","red"),ylim=c(0,1))
#   print(c(sum(SNPe3[,3]<0.05)/500,sum(SNPs3[,3]<0.05)/500,sum(SNPs3[,7]<0.05)/500,sum(PLEs3[,3]<0.05)/500,sum(PLEe3[,3]<0.05)/500))
# }
# 
# Rdatatranst<-function(dat){
#   SNPs=dat[seq(1,9000,18),]
#   SNPs2=dat[seq(2,9000,18),]
#   SNPs3=dat[seq(3,9000,18),]
#   
#   SNPe=dat[seq(4,9000,18),]
#   SNPe2=dat[seq(5,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   
#   PLEs=dat[seq(7,9000,18),]
#   PLEs2=dat[seq(8,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   
#   PLEe=dat[seq(10,9000,18),]
#   PLEe2=dat[seq(11,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   
#   FIXd=dat[seq(13,9000,18),]
#   FIXd2=dat[seq(13,9000,18),]
#   FIXd3=dat[seq(13,9000,18),]
#   
#   FULL=dat[seq(16,9000,18),]
#   FULL2=dat[seq(17,9000,18),]
#   FULL3=dat[seq(18,9000,18),]
#   
#   par(mfrow=c(5,3))
#   plot(SNPe[,1],col=ifelse(SNPe3[,1]>0.05,"green","blue"),pch=2,main="Parameter Value",ylab = "SNPtest(e)")
#   plot(SNPe2[,1],col=ifelse(SNPe3[,1]>0.05,"green","blue"),pch=2,main="SE Value")
#   plot(SNPe3[,1],col=ifelse(SNPe3[,1]>0.05,"green","blue"),pch=2,main="P-value",ylim=c(0,1))
#   plot(SNPs[,1],col=ifelse(SNPs3[,1]>0.05,"green","blue"),ylab = "SNPtest(s)")
#   plot(SNPs2[,1],col=ifelse(SNPs3[,1]>0.05,"green","blue"))
#   plot(SNPs3[,1],col=ifelse(SNPs3[,1]>0.05,"green","blue"),ylim=c(0,1))
#   plot(SNPe[,5],col=ifelse(SNPe3[,5]>0.05,"green","blue"),pch=3,ylab = "SNPtest(sus.dev)")
#   plot(SNPe2[,5],col=ifelse(SNPe3[,5]>0.05,"green","blue"),pch=3)
#   plot(SNPe3[,5],col=ifelse(SNPe3[,5]>0.05,"green","blue"),pch=3,ylim=c(0,1))
#   plot(PLEe[,1],col=ifelse(PLEe3[,1]>0.05,"green","red"),pch=2,ylab = "Pleiotropy (e+sus.dev)")
#   plot(PLEe2[,1],col=ifelse(PLEe3[,1]>0.05,"green","red"),pch=2)
#   plot(PLEe3[,1],col=ifelse(PLEe3[,1]>0.05,"green","red"),pch=2,ylim=c(0,1))
#   plot(PLEs[,1],col=ifelse(PLEs3[,1]>0.05,"green","red"),ylab = "Pleiotropy (s+sus.dev)")
#   plot(PLEs2[,1],col=ifelse(PLEs3[,1]>0.05,"green","red"))
#   plot(PLEs3[,1],col=ifelse(PLEs3[,1]>0.05,"green","red"),ylim=c(0,1))
#   print(c(sum(SNPe3[,1]<0.05)/500,sum(SNPs3[,1]<0.05)/500,sum(SNPs3[,5]<0.05)/500,sum(PLEs3[,1]<0.05)/500,sum(PLEe3[,1]<0.05)/500))
# }
# 
# Rdatatrans2t<-function(dat){
#   SNPs=dat[seq(1,9000,18),]
#   SNPs2=dat[seq(2,9000,18),]
#   SNPs3=dat[seq(3,9000,18),]
#   
#   SNPe=dat[seq(4,9000,18),]
#   SNPe2=dat[seq(5,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   
#   PLEs=dat[seq(7,9000,18),]
#   PLEs2=dat[seq(8,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   
#   PLEe=dat[seq(10,9000,18),]
#   PLEe2=dat[seq(11,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   
#   FIXd=dat[seq(13,9000,18),]
#   FIXd2=dat[seq(13,9000,18),]
#   FIXd3=dat[seq(13,9000,18),]
#   
#   FULL=dat[seq(16,9000,18),]
#   FULL2=dat[seq(17,9000,18),]
#   FULL3=dat[seq(18,9000,18),]
#   
#   par(mfrow=c(5,3))
#   plot(SNPe[,2],col=ifelse(SNPe3[,2]>0.05,"green","blue"),pch=2,main="Parameter Value",ylab = "SNPtest(e)")
#   plot(SNPe2[,2],col=ifelse(SNPe3[,2]>0.05,"green","blue"),pch=2,main="SE Value")
#   plot(SNPe3[,2],col=ifelse(SNPe3[,2]>0.05,"green","blue"),pch=2,main="P-value",ylim=c(0,1))
#   plot(SNPs[,2],col=ifelse(SNPs3[,2]>0.05,"green","blue"),ylab = "SNPtest(s)")
#   plot(SNPs2[,2],col=ifelse(SNPs3[,2]>0.05,"green","blue"))
#   plot(SNPs3[,2],col=ifelse(SNPs3[,2]>0.05,"green","blue"),ylim=c(0,1))
#   plot(SNPe[,6],col=ifelse(SNPe3[,6]>0.05,"green","blue"),pch=3,ylab = "SNPtest(sus.dev)")
#   plot(SNPe2[,6],col=ifelse(SNPe3[,6]>0.05,"green","blue"),pch=3)
#   plot(SNPe3[,6],col=ifelse(SNPe3[,6]>0.05,"green","blue"),pch=3,ylim=c(0,1))
#   plot(PLEe[,2],col=ifelse(PLEe3[,2]>0.05,"green","red"),pch=2,ylab = "Pleiotropy (e+sus.dev)")
#   plot(PLEe2[,2],col=ifelse(PLEe3[,2]>0.05,"green","red"),pch=2)
#   plot(PLEe3[,2],col=ifelse(PLEe3[,2]>0.05,"green","red"),pch=2,ylim=c(0,1))
#   plot(PLEs[,2],col=ifelse(PLEs3[,2]>0.05,"green","red"),ylab = "Pleiotropy (s+sus.dev)")
#   plot(PLEs2[,2],col=ifelse(PLEs3[,2]>0.05,"green","red"))
#   plot(PLEs3[,2],col=ifelse(PLEs3[,2]>0.05,"green","red"),ylim=c(0,1))
#   print(c(sum(SNPe3[,2]<0.05)/500,sum(SNPs3[,2]<0.05)/500,sum(SNPs3[,6]<0.05)/500,sum(PLEs3[,2]<0.05)/500,sum(PLEe3[,2]<0.05)/500))
# }
# 
# Rdatatrans3t<-function(dat){
#   SNPs=dat[seq(1,9000,18),]
#   SNPs2=dat[seq(2,9000,18),]
#   SNPs3=dat[seq(3,9000,18),]
#   
#   SNPe=dat[seq(4,9000,18),]
#   SNPe2=dat[seq(5,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   
#   PLEs=dat[seq(7,9000,18),]
#   PLEs2=dat[seq(8,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   
#   PLEe=dat[seq(10,9000,18),]
#   PLEe2=dat[seq(11,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   
#   FIXd=dat[seq(13,9000,18),]
#   FIXd2=dat[seq(13,9000,18),]
#   FIXd3=dat[seq(13,9000,18),]
#   
#   FULL=dat[seq(16,9000,18),]
#   FULL2=dat[seq(17,9000,18),]
#   FULL3=dat[seq(18,9000,18),]
#   
#   par(mfrow=c(5,3))
#   plot(SNPe[,3],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="Parameter Value",ylab = "SNPtest(e)")
#   plot(SNPe2[,3],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="SE Value")
#   plot(SNPe3[,3],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="P-value",ylim=c(0,1))
#   plot(SNPs[,3],col=ifelse(SNPs3[,3]>0.05,"green","blue"),ylab = "SNPtest(s)")
#   plot(SNPs2[,3],col=ifelse(SNPs3[,3]>0.05,"green","blue"))
#   plot(SNPs3[,3],col=ifelse(SNPs3[,3]>0.05,"green","blue"),ylim=c(0,1))
#   plot(SNPe[,7],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3,ylab = "SNPtest(sus.dev)")
#   plot(SNPe2[,7],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3)
#   plot(SNPe3[,7],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3,ylim=c(0,1))
#   plot(PLEe[,3],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2,ylab = "Pleiotropy (e+sus.dev)")
#   plot(PLEe2[,3],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2)
#   plot(PLEe3[,3],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2,ylim=c(0,1))
#   plot(PLEs[,3],col=ifelse(PLEs3[,3]>0.05,"green","red"),ylab = "Pleiotropy (s+sus.dev)")
#   plot(PLEs2[,3],col=ifelse(PLEs3[,3]>0.05,"green","red"))
#   plot(PLEs3[,3],col=ifelse(PLEs3[,3]>0.05,"green","red"),ylim=c(0,1))
#   print(c(sum(SNPe3[,3]<0.05)/500,sum(SNPs3[,3]<0.05)/500,sum(SNPs3[,7]<0.05)/500,sum(PLEs3[,3]<0.05)/500,sum(PLEe3[,3]<0.05)/500))
# }
# par(mar = rep(2, 4))
# Rdatatrans(Rdataneg25neg5)
# Rdatatrans2(Rdataneg25neg5)
# Rdatatrans3(Rdataneg25neg5)
# Rdatatranst(Rdataneg25neg5)
# Rdatatrans2t(Rdatapos5pos25)
# Rdatatrans3t(Rdataneg25neg5)
# 
# Rdatatrans(Rdataneg25neg5)
# 
# Rdatatrans(Rdataneg25neg25)
# Rdatatrans(Rdatapos25pos25)
# Rdatatrans(Rdataneg25neg5)
# Rdatatrans(Rdatapos25pos5)
# Rdatatrans(Rdatapos25pos25)
# 
# Rdatatrans3t<-function(dat){
#   SNPs=dat[seq(1,9000,18),]
#   SNPs2=dat[seq(2,9000,18),]
#   SNPs3=dat[seq(3,9000,18),]
#   
#   SNPe=dat[seq(4,9000,18),]
#   SNPe2=dat[seq(5,9000,18),]
#   SNPe3=dat[seq(6,9000,18),]
#   
#   PLEs=dat[seq(7,9000,18),]
#   PLEs2=dat[seq(8,9000,18),]
#   PLEs3=dat[seq(9,9000,18),]
#   
#   PLEe=dat[seq(10,9000,18),]
#   PLEe2=dat[seq(11,9000,18),]
#   PLEe3=dat[seq(12,9000,18),]
#   
#   FIXd=dat[seq(13,9000,18),]
#   FIXd2=dat[seq(13,9000,18),]
#   FIXd3=dat[seq(13,9000,18),]
#   
#   FULL=dat[seq(16,9000,18),]
#   FULL2=dat[seq(17,9000,18),]
#   FULL3=dat[seq(18,9000,18),]
#   
#   par(mfrow=c(5,3))
#   plot(SNPe[,3],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="Parameter Value",ylab = "SNPtest(e)")
#   plot(SNPe2[,3],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="SE Value")
#   plot(SNPe3[,3],col=ifelse(SNPe3[,3]>0.05,"green","blue"),pch=2,main="P-value",ylim=c(0,1))
#   plot(SNPs[,3],col=ifelse(SNPs3[,3]>0.05,"green","blue"),ylab = "SNPtest(s)")
#   plot(SNPs2[,3],col=ifelse(SNPs3[,3]>0.05,"green","blue"))
#   plot(SNPs3[,3],col=ifelse(SNPs3[,3]>0.05,"green","blue"),ylim=c(0,1))
#   plot(SNPe[,7],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3,ylab = "SNPtest(sus.dev)")
#   plot(SNPe2[,7],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3)
#   plot(SNPe3[,7],col=ifelse(SNPe3[,7]>0.05,"green","blue"),pch=3,ylim=c(0,1))
#   plot(PLEe[,3],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2,ylab = "Pleiotropy (e+sus.dev)")
#   plot(PLEe2[,3],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2)
#   plot(PLEe3[,3],col=ifelse(PLEe3[,3]>0.05,"green","red"),pch=2,ylim=c(0,1))
#   plot(PLEs[,3],col=ifelse(PLEs3[,3]>0.05,"green","red"),ylab = "Pleiotropy (s+sus.dev)")
#   plot(PLEs2[,3],col=ifelse(PLEs3[,3]>0.05,"green","red"))
#   plot(PLEs3[,3],col=ifelse(PLEs3[,3]>0.05,"green","red"),ylim=c(0,1))
#   print(c(sum(SNPe3[,3]<0.05)/500,sum(SNPs3[,3]<0.05)/500,sum(SNPs3[,7]<0.05)/500,sum(PLEs3[,3]<0.05)/500,sum(PLEe3[,3]<0.05)/500))
# }

