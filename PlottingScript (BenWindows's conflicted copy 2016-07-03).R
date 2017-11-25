setwd("C:\\Users/Ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/")
setwd("/home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/")
data100=read.csv("PvalueResults_100_0.4.txt",header=FALSE)
data250=read.csv("PvalueResults_250_0.4.txt",header=FALSE)
data500=read.csv("PvalueResults_500_0.4.txt",header=FALSE)
data1000=read.csv("PvalueResults_1000_0.4.txt",header=FALSE)
data2000=read.csv("PvalueResults_2000_0.4.txt",header=FALSE)

#Convert Function
datatrans<-function(dat){
Null1a=c((dat[1:6,2]+dat[1:6,3])/1000,(dat[1:6,6]+dat[1:6,7])/1000)
alp1=Null1a[c(1,2,7,3,4,5,11,6,12)]
Null1p=rep(0,12)#c((dat[1:6,2]+dat[1:6,3])/1000,(dat[1:6,6]+dat[1:6,7])/1000)
pow1=Null1p[c(1,2,7,3,4,5,11,6,12)]

Lu1a=c((dat[7:12,2]+dat[7:12,3])/1000,(dat[7:12,6]+dat[7:12,7])/1000)
alp2=Lu1a[c(1,2,7,3,4,5,11,6,12)]
Lu1p=rep(0,12)#c((dat[7:12,2]+dat[7:12,3])/1000,(dat[7:12,6]+dat[7:12,7])/1000)
pow2=Lu1p[c(1,2,7,3,4,5,11,6,12)]

Surv1a=c((dat[13:18,3])/500,(dat[13:18,6]+dat[13:18,7])/1000)
alp3=Surv1a[c(1,2,7,3,4,5,11,6,12)]
Surv1p=c((dat[13:18,2])/500,rep(0,6))
pow3=Surv1p[c(1,2,7,3,4,5,11,6,12)]

Sus1a=c((dat[19:24,2]+dat[19:24,3])/1000,(dat[19:24,7])/500)
alp4=Sus1a[c(1,2,7,3,4,5,11,6,12)]
Sus1p=c(rep(0,6),(dat[19:24,6])/500)
pow4=Sus1p[c(1,2,7,3,4,5,11,6,12)]

Both1a=c((dat[25:30,3])/500,(dat[25:30,7])/500)
alp5=Both1a[c(1,2,7,3,4,5,11,6,12)]
Both1p=c((dat[25:30,2])/500,(dat[25:30,6])/500)
pow5=Both1p[c(1,2,7,3,4,5,11,6,12)]

alp=c(alp1,alp2,alp3,alp4,alp5)
pow=c(pow1,pow2,pow3,pow4,pow5)
return(cbind(alp,pow))
}

datar100=datatrans(data100)
datar250=datatrans(data250)
datar500=datatrans(data500)
datar1000=datatrans(data1000)
datar2000=datatrans(data2000)




png("Sig_0_4_MAF.png",width=1200,height=900)
par(mfrow=c(2,3))
st=1
fin=9
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Sample Size",main = "\nNull model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
abline(h=0.05,lty=2)
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

st=10
fin=18
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Sample Size",main = "Significance with 0.4 MAF\nClinical Only model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
abline(h=0.05,lty=2)
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

st=19
fin=27
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Sample Size",main = "\nGenetic Survival model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
abline(h=0.05,lty=2)
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

st=28
fin=36
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Sample Size",main = "\nGenetic Susceptibility model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
abline(h=0.05,lty=2)
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

st=37
fin=45
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,0.4),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Sample Size",main = "\nGenetic Pleiotropy model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:7, labels=c(100,250,500,1000,2000,5000,10000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
#points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
abline(h=0.05,lty=2)
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

dev.off()


png("Power_0_4_MAF.png",width=1600,height=450)
par(mfrow=c(1,4))
st=19
fin=27
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,1),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Sample Size",main = "\nGenetic Survival model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(100,250,500,1000,2000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))
abline(h=0.8,lty=2)

st=28
fin=36
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,1),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Sample Size",main = "Power with 0.4 MAF\nGenetic Susceptibility model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(100,250,500,1000,2000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))
abline(h=0.8,lty=2)

st=37
fin=45
plot(rep(1,9),datar100[st:fin,2], ylim=c(0,1),xlim=c(0,8),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Sample Size",main = "\nGenetic Pleiotropy model",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(100,250,500,1000,2000))
points(rep(2,9),datar250[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(3,9),datar500[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
points(rep(4,9),datar1000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
#points(rep(5,9),datar2000[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
#legend("topright",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))
abline(h=0.8,lty=2)

plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

dev.off()


####Vary parameters

datapos5pos5=read.csv("PvalueResults0.5_0.5_2000_0.4.txt",header=FALSE)
datapos5pos25=read.csv("PvalueResults0.5_0.25_2000_0.4.txt",header=FALSE)
datapos50=read.csv("PvalueResults0.5_0_2000_0.4.txt",header=FALSE)
datapos5neg25=read.csv("PvalueResults0.5_-0.25_2000_0.4.txt",header=FALSE)
datapos5neg5=read.csv("PvalueResults0.5_-0.5_2000_0.4.txt",header=FALSE)

datapos25pos5=read.csv("PvalueResults0.25_0.5_2000_0.4.txt",header=FALSE)
datapos25pos25=read.csv("PvalueResults0.25_0.25_2000_0.4.txt",header=FALSE)
datapos250=read.csv("PvalueResults0.25_0_2000_0.4.txt",header=FALSE)
datapos25neg25=read.csv("PvalueResults0.25_-0.25_2000_0.4.txt",header=FALSE)
datapos25neg5=read.csv("PvalueResults0.25_-0.5_2000_0.4.txt",header=FALSE)

data0pos5=read.csv("PvalueResults0_0.5_2000_0.4.txt",header=FALSE)
data0pos25=read.csv("PvalueResults0_0.25_2000_0.4.txt",header=FALSE)
data00=read.csv("PvalueResults0_0_2000_0.4.txt",header=FALSE)
data0neg25=read.csv("PvalueResults0_-0.25_2000_0.4.txt",header=FALSE)
data0neg5=read.csv("PvalueResults0_-0.5_2000_0.4.txt",header=FALSE)

dataneg25pos5=read.csv("PvalueResults-0.25_0.5_2000_0.4.txt",header=FALSE)
dataneg25pos25=read.csv("PvalueResults-0.25_0.25_2000_0.4.txt",header=FALSE)
dataneg250=read.csv("PvalueResults-0.25_0_2000_0.4.txt",header=FALSE)
dataneg25neg25=read.csv("PvalueResults-0.25_-0.25_2000_0.4.txt",header=FALSE)
dataneg25neg5=read.csv("PvalueResults-0.25_-0.5_2000_0.4.txt",header=FALSE)

dataneg5pos5=read.csv("PvalueResults-0.5_0.5_2000_0.4.txt",header=FALSE)
dataneg5pos25=read.csv("PvalueResults-0.5_0.25_2000_0.4.txt",header=FALSE)
dataneg50=read.csv("PvalueResults-0.5_0_2000_0.4.txt",header=FALSE)
dataneg5neg25=read.csv("PvalueResults-0.5_-0.25_2000_0.4.txt",header=FALSE)
dataneg5neg5=read.csv("PvalueResults-0.5_-0.5_2000_0.4.txt",header=FALSE)


#Convert Function
datatrans1<-function(dat){
  Both1a=c((dat[1:6,3])/500,(dat[1:6,7])/500)
  alp1=Both1a[c(1,2,7,3,4,5,11,6,12)]
  Both1p=c((dat[25:30,2])/500,(dat[25:30,6])/500)
  pow1=Both1p[c(1,2,7,3,4,5,11,6,12)]
  
  Both1a=c((dat[7:12,3])/500,(dat[7:12,7])/500)
  alp2=Both1a[c(1,2,7,3,4,5,11,6,12)]
  Both1p=c((dat[7:12,2])/500,(dat[7:12,6])/500)
  pow2=Both1p[c(1,2,7,3,4,5,11,6,12)]
  
  Surv1a=c((dat[13:18,3])/500,(dat[13:18,6]+dat[13:18,7])/1000)
  alp3=Surv1a[c(1,2,7,3,4,5,11,6,12)]
  Surv1p=c((dat[13:18,2])/500,rep(0,6))
  pow3=Surv1p[c(1,2,7,3,4,5,11,6,12)]
  
  Both1a=c((dat[19:24,3])/500,(dat[19:24,7])/500)
  alp4=Both1a[c(1,2,7,3,4,5,11,6,12)]
  Both1p=c((dat[19:24,2])/500,(dat[19:24,6])/500)
  pow4=Both1p[c(1,2,7,3,4,5,11,6,12)]
  
  Both1a=c((dat[25:30,3])/500,(dat[25:30,7])/500)
  alp5=Both1a[c(1,2,7,3,4,5,11,6,12)]
  Both1p=c((dat[25:30,2])/500,(dat[25:30,6])/500)
  pow5=Both1p[c(1,2,7,3,4,5,11,6,12)]
  
  alp=c(alp1,alp2,alp3,alp4,alp5)
  pow=c(pow1,pow2,pow3,pow4,pow5)
  return(cbind(alp,pow))
}

dataneg5=datatrans1(rbind(dataneg5pos5,dataneg5pos25,dataneg50,dataneg5neg25,dataneg5neg5))
dataneg25=datatrans1(rbind(dataneg25pos5,dataneg25pos25,dataneg250,dataneg25neg25,dataneg25neg25))
data0=datatrans1(rbind(data0pos5,data0pos25,data00,data0neg25,data00))
datapos25=datatrans1(rbind(datapos25pos5,datapos25pos25,datapos250,datapos25neg25,datapos25pos25))
datapos5=datatrans1(rbind(datapos5pos5,datapos5pos25,datapos50,datapos5neg25,datapos5pos5))

png("Sig_0_4_MAF.png",width=1200,height=900)
par(mfrow=c(2,3))
st=1
fin=9
plot(rep(1,9),dataneg5[st:fin,2], ylim=c(0,0.1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),dataneg25[st:fin,2], ylim=c(0,0.1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),data0[st:fin,2], ylim=c(0,0.1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),datapos25[st:fin,2], ylim=c(0,0.1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),datapos5[st:fin,2], ylim=c(0,0.1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="Alpha",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))



plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

dev.off()


png("Power_0_4_MAF.png",width=1600,height=450)
par(mfrow=c(2,3))
st=1
fin=9
plot(rep(1,9),dataneg5[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),dataneg5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),dataneg25[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),dataneg25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),data0[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = 0",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),data0[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),datapos25[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.25",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),datapos25[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))

st=1
fin=9
plot(rep(1,9),datapos5[st:fin,2], ylim=c(0,1),xlim=c(0,6),col=c("blue","blue","blue","red","red","gray","gray","black","black"),ylab="1-Beta",xlab="Susceptibility Effect",main = "\nSurvival Effect = -0.5",xaxt="n",pch=c(1,2,3,1,2,1,3,1,3))
axis(1, at=1:5, labels=c(-0.5,-0.25,0,0.25,0.5))
abline(h=0.05,lty=2)
st=10
fin=18
points(rep(2,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=19
fin=27
points(rep(3,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=28
fin=36
points(rep(4,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))
st=37
fin=45
points(rep(5,9),datapos5[st:fin,2],col=c("blue","blue","blue","red","red","gray","gray","black","black"),pch=c(1,2,3,1,2,1,3,1,3))



plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft",legend=c("SNPtest(s)","SNPtest(e)","SNPtest(sus)","Pleio(s)","Pleio(e)","Fix(surv)","Fix(sus)","Full(surv)","Full(sus)"),pch=c(1,2,3,1,2,1,3,1,3),col=c("blue","blue","blue","red","red","gray","gray","black","black"))

dev.off()
