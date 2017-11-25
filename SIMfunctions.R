#system("scp -r -P 8022 /home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/SIMfunctions.R gaia:/work/projects/epipgx/users/LIV/Task1Surv/Simulation/")

library(MASS)
library(survival)
library(foreach)
library(doParallel)

res_coll = function(n,maf){
  Sim=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/NullSimResults_",n,"_",maf,".txt",sep=""))
  size=dim(Sim)[1]
  p1=rbind(apply(Sim[seq(3,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(6,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(9,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(12,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(18,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(15,size,18),],2,function(x){sum(x<0.05)}))
  biasu=rbind(apply(Sim[seq(1,size,18),],2,function(x){mean(x)}),apply(Sim[seq(4,size,18),],2,function(x){mean(x)}),apply(Sim[seq(7,size,18),],2,function(x){mean(x)}),apply(Sim[seq(10,size,18),],2,function(x){mean(x)}),apply(Sim[seq(16,size,18),],2,function(x){mean(x)}),apply(Sim[seq(13,size,18),],2,function(x){mean(x)}))
  bias1=cbind(biasu[,1]-0,biasu[,2]-0,biasu[,3]-0,biasu[,4]-0,biasu[,5]+0,biasu[,6]-0,biasu[,7]-0)
  sd1=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)}))
  se1=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)/sqrt(n)}))
  
  Sim=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/LuSimResults_",n,"_",maf,".txt",sep=""))#beta0 = c(1.0, 0, 0),gamma0 = c(0.5,-0.5,0,0)
  size=dim(Sim)[1]
  p2=rbind(apply(Sim[seq(3,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(6,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(9,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(12,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(18,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(15,size,18),],2,function(x){sum(x<0.05)}))
  biasu=rbind(apply(Sim[seq(1,size,18),],2,function(x){mean(x)}),apply(Sim[seq(4,size,18),],2,function(x){mean(x)}),apply(Sim[seq(7,size,18),],2,function(x){mean(x)}),apply(Sim[seq(10,size,18),],2,function(x){mean(x)}),apply(Sim[seq(16,size,18),],2,function(x){mean(x)}),apply(Sim[seq(13,size,18),],2,function(x){mean(x)}))
  bias2=cbind(biasu[,1]-1,biasu[,2]-0,biasu[,3]-0,biasu[,4]-0.5,biasu[,5]+0.5,biasu[,6]-0,biasu[,7]-0)
  sd2=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)}))
  se2=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)/sqrt(n)}))
  
  Sim=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Lu_SurvResults_",n,"_",maf,".txt",sep=""))#beta0 = c(1.0, 0.5, 0),gamma0 = c(0.5,-0.5,0,0)
  size=dim(Sim)[1]
  p3=rbind(apply(Sim[seq(3,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(6,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(9,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(12,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(18,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(15,size,18),],2,function(x){sum(x<0.05)}))
  biasu=rbind(apply(Sim[seq(1,size,18),],2,function(x){mean(x)}),apply(Sim[seq(4,size,18),],2,function(x){mean(x)}),apply(Sim[seq(7,size,18),],2,function(x){mean(x)}),apply(Sim[seq(10,size,18),],2,function(x){mean(x)}),apply(Sim[seq(16,size,18),],2,function(x){mean(x)}),apply(Sim[seq(13,size,18),],2,function(x){mean(x)}))
  bias3=cbind(biasu[,1]-1,biasu[,2]-0.5,biasu[,3]-0,biasu[,4]-0.5,biasu[,5]+0.5,biasu[,6]-0,biasu[,7]-0)
  sd3=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)}))
  se3=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)/sqrt(n)}))
  
  Sim=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Lu_SusResults_",n,"_",maf,".txt",sep=""))#beta0 = c(1.0, 0, 0),gamma0 = c(0.5,-0.5,0.5,0)
  size=dim(Sim)[1]
  p4=rbind(apply(Sim[seq(3,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(6,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(9,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(12,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(18,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(15,size,18),],2,function(x){sum(x<0.05)}))
  biasu=rbind(apply(Sim[seq(1,size,18),],2,function(x){mean(x)}),apply(Sim[seq(4,size,18),],2,function(x){mean(x)}),apply(Sim[seq(7,size,18),],2,function(x){mean(x)}),apply(Sim[seq(10,size,18),],2,function(x){mean(x)}),apply(Sim[seq(16,size,18),],2,function(x){mean(x)}),apply(Sim[seq(13,size,18),],2,function(x){mean(x)}))
  bias4=cbind(biasu[,1]-1,biasu[,2]-0,biasu[,3]-0,biasu[,4]-0.5,biasu[,5]+0.5,biasu[,6]-0.5,biasu[,7]-0)
  sd4=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)}))
  se4=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)/sqrt(n)}))
  
  Sim=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Lu_BothResults_",n,"_",maf,".txt",sep=""))#beta0 = c(1.0, -0.5, 0),gamma0 = c(0.5,-0.5,0.5,0)
  size=dim(Sim)[1]
  p5=rbind(apply(Sim[seq(3,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(6,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(9,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(12,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(18,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(15,size,18),],2,function(x){sum(x<0.05)}))
  biasu=rbind(apply(Sim[seq(1,size,18),],2,function(x){mean(x)}),apply(Sim[seq(4,size,18),],2,function(x){mean(x)}),apply(Sim[seq(7,size,18),],2,function(x){mean(x)}),apply(Sim[seq(10,size,18),],2,function(x){mean(x)}),apply(Sim[seq(16,size,18),],2,function(x){mean(x)}),apply(Sim[seq(13,size,18),],2,function(x){mean(x)}))
  bias5=cbind(biasu[,1]-1,biasu[,2]+0.5,biasu[,3]-0,biasu[,4]-0.5,biasu[,5]+0.5,biasu[,6]-0.5,biasu[,7]-0)
  sd5=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)}))
  se5=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)/sqrt(n)}))
  
  p=rbind(p1,p2,p3,p4,p5)
  bias=rbind(bias1,bias2,bias3,bias4,bias5)
  sd=rbind(sd1,sd2,sd3,sd4,sd5)
  
  write.table(p,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/PvalueResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F,sep=",")
  write.table(bias,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/BiasResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F,sep=",")
  write.table(sd,paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SDResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F,sep=",")
}
res_coll2 = function(n,maf,survb,susg){
  Sim=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Results",survb,"_",susg,"_",n,"_",maf,".txt",sep=""))
  size=dim(Sim)[1]
  p1=rbind(apply(Sim[seq(3,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(6,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(9,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(12,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(18,size,18),],2,function(x){sum(x<0.05)}),apply(Sim[seq(15,size,18),],2,function(x){sum(x<0.05)}))
  biasu=rbind(apply(Sim[seq(1,size,18),],2,function(x){mean(x)}),apply(Sim[seq(4,size,18),],2,function(x){mean(x)}),apply(Sim[seq(7,size,18),],2,function(x){mean(x)}),apply(Sim[seq(10,size,18),],2,function(x){mean(x)}),apply(Sim[seq(16,size,18),],2,function(x){mean(x)}),apply(Sim[seq(13,size,18),],2,function(x){mean(x)}))
  bias1=cbind(biasu[,1]-1,biasu[,2]-survb,biasu[,3]-0,biasu[,4]-0.5,biasu[,5]+0.5,biasu[,6]-susg,biasu[,7]-0)
  sd1=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)}))
  se1=rbind(apply(Sim[seq(2,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(5,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(8,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(11,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(17,size,18),],2,function(x){mean(x)/sqrt(n)}),apply(Sim[seq(14,size,18),],2,function(x){mean(x)/sqrt(n)}))

  write.table(p1,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/PvalueResults",survb,"_",susg,"_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F,sep=",")
  write.table(bias1,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/BiasResults",survb,"_",susg,"_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F,sep=",")
  write.table(sd1,paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SDResults",survb,"_",susg,"_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F,sep=",")
}

s_ext=function(s,e,delta){
  s = s[dense_rank(e)]
  s[e>max(e[delta==1])] = 0
  return(s)
}

expit = function(x) exp(x)/(1+exp(x))

pb.f = function(n,gamma,z){
  x = cbind(rep(1,n),z)
  a = exp(x%*%gamma)
  return(a/(1+a))
}

fpb.f = function(n,gamma,gamma.clin,z){
  a = exp((z*gamma)+gamma.clin)
  return(a/(1+a))
}



eta_prob = function(gamma,z,delta,eS){
  n = length(delta)
  pb = pb.f(n,gamma,z)
  w = delta + (1-delta)*(pb*eS)/(1-pb+pb*eS)
  return(w)
}

feta_prob = function(gamma,gamma.clin,z,delta,eS){
  n = length(delta)
  pb = fpb.f(n,gamma,gamma.clin,z)
  w = delta + (1-delta)*(pb*eS)/(1-pb+pb*eS)
  return(w)
}

post.lik1 = function(gamma,w,z){
  n = length(w)
  pb = pb.f(n,gamma,z)
  a=-sum(w*log(pb)+(1-w)*log(1-pb))
  return(a)
}

fpost.lik1 = function(gamma,gamma.clin,w,z){
  n = length(w)
  pb = fpb.f(n,gamma,gamma.clin,z)
  a=-sum(w*log(pb)+(1-w)*log(1-pb))
  return(a)
}

post.lik11 = function(gamma1,gamma2,w,z){
  n = length(w)
  gamma=c(gamma1,gamma2)
  pb = pb.f(n,gamma,z)
  a=-sum(w*log(pb)+(1-w)*log(1-pb))
  return(a)
} 

post.lik12 = function(gamma2,gamma1,w,z){
  n = length(w)
  gamma=c(gamma1,gamma2)
  pb = pb.f(n,gamma,z)
  a=-sum(w*log(pb)+(1-w)*log(1-pb))
  return(a)
} 

post.lik2 = function(beta,obs.time,z,delta,w,h){
  n = length(obs.time)
  if(length(beta)==1){
    x = log(obs.time) - beta*z}
  if(length(beta)>1){
    x = log(obs.time) - z%*%beta}
  xm = matrix(rep(x,n),nrow=n)
  dxm = xm - t(xm)
  f = apply(dnorm(dxm,sd=h)*delta,2,mean)
  S = apply(pnorm(dxm,sd=h)*w,2,mean)
  a = -sum(delta*(log(f)-log(S)))
  return(a)
} 

fpost.lik2 = function(beta,beta.clin,obs.time,z,delta,w,h){
  n = length(obs.time)
  if(length(beta)==1){
    x = log(obs.time) - beta*z-beta.clin}
  if(length(beta)>1){
    x = log(obs.time) - z%*%beta-beta.clin}
  xm = matrix(rep(x,n),nrow=n)
  dxm = xm - t(xm)
  f = apply(dnorm(dxm,sd=h)*delta,2,mean)
  S = apply(pnorm(dxm,sd=h)*w,2,mean)
  a = -sum(delta*(log(f)-log(S)))
  return(a)
} 

post.lik1n = function(gamma1,gamma2,k,w,z){
  n = length(w)
  gamma = numeric(length(z[1,])+1)
  gamma[k] = gamma2
  gamma[-k] = gamma1
  pb = pb.f(n,gamma,z)
  a=-sum(w*log(pb)+(1-w)*log(1-pb))
  return(a)
} 

post.lik2n = function(beta1,beta2,k,obs.time,z,delta,w,h){
  beta = numeric(length(z[1,]))
  beta[k] = beta2
  beta[-k] = beta1
  n = length(obs.time)
  x = log(obs.time) - z%*%beta
  xm = matrix(rep(x,n),nrow=n)
  dxm = xm - t(xm)
  f = apply(dnorm(dxm,sd=h)*delta,2,mean)
  S = apply(pnorm(dxm,sd=h)*w,2,mean)
  a = -sum(delta*(log(f)-log(S)))
  return(a)
} 

lam.f = function(x,e,delta,w,h){
  f = mean(dnorm(e-x,sd=h)*delta)
  S = mean(pnorm(e-x,sd=h)*w)
  lambda = f/S
  return(lambda)
}

dlam.f = function(x,e,delta,w,h){
  f = mean(dnorm(e-x,sd=h)*delta)
  S = mean(pnorm(e-x,sd=h)*w)
  fd = mean(dnorm(e-x,sd=h)*(e-x)*delta/h^2)
  Sd = mean(dnorm(e-x,sd=h)*w)
  dlambda = (fd*S-f*Sd)/(S*S)
  return(dlambda)
} 

Lam.f = function(e,delta,w,h){
  n=length(w)
  Lambda = numeric(n)
  temp = numeric(n)
  se = sort(e)
  a = se[1] - 2.0
  b = se[1]
  fa = lam.f(a,e,delta,w,h)
  fb = lam.f(b,e,delta,w,h)
  temp[1] = area(lam.f,a,b,e=e,delta=delta,w=w,h=h,fa=fa,fb=fb)
  for(i in 2:n){
    a = se[i-1]
    b = se[i]
    fa = lam.f(a,e,delta,w,h)
    fb = lam.f(b,e,delta,w,h)
    temp[i] = area(lam.f,a,b,e=e,delta=delta,w=w,h=h,fa=fa,fb=fb)
  }
  Lambda = cumsum(temp)
  return(Lambda)
}

pl.lik = function(beta,gamma,obs.time,z1,z2,delta,w,h){
  n = length(obs.time)
  lam = numeric(n)
    pb = pb.f(n,gamma,z2)
    a = w*log(pb)+(1-w)*log(1-pb)
  if(length(beta)==1){
    e = log(obs.time) - beta*z1}
  if(length(beta)>1){
    e = log(obs.time) - z1%*%beta}
  for(i in 1:n){
    x = e[i]
    lam[i] = lam.f(x,e,delta,w,h)
  }
  Lam = Lam.f(e,delta,w,h)
  Lam = Lam[dense_rank(e)]
  a = a+delta*log(lam)-w*Lam
  return(a)
}

fpl.lik = function(beta,gamma,beta.clin,gamma.clin,obs.time,z,delta,w,h){
  n = length(obs.time)
  lam = numeric(n)
  pb = fpb.f(n,gamma,gamma.clin,z)
  a = w*log(pb)+(1-w)*log(1-pb)
  if(length(beta)==1){
    e = log(obs.time) - beta*z - beta.clin}
  if(length(beta)>1){
    e = log(obs.time) - z%*%beta - beta.clin}
  for(i in 1:n){
    x = e[i]
    lam[i] = lam.f(x,e,delta,w,h)
  }
  Lam = Lam.f(e,delta,w,h)
  Lam = Lam[dense_rank(e)]
  a = a+delta*log(lam)-w*Lam
  return(a)
}


pl.score = function(beta,gamma,obs.time,z,delta,w,h){
  n = length(obs.time)
  smat = matrix(0,nrow=n,ncol=3)
  dlam = numeric(n)
  lam = numeric(n)
  pb = pb.f(n,gamma,z)
  smat[,2] = w*(1-pb)-(1-w)*pb
  smat[,3] = z*(w*(1-pb)-(1-w)*pb)
  e = log(obs.time) - beta*z
  for(i in 1:n){
    x = e[i]
    lam[i] = lam.f(x,e,delta,w,h)
    dlam[i] = dlam.f(x,e,delta,w,h)
  }
  smat[,1] = -z*(delta*dlam/lam-w*lam)
  return(smat)
} 

SNPgen <- function(n=1000, prop_SNP1=c(0.5625,0.375,0.0625), prop_SNP2=c(0.5625,0.375,0.0625)){
  snp1 <- sample(size=n, c("AA","Aa","aa"), prob=prop_SNP1, replace=TRUE) 
  snp2 <- sample(size=n, c("BB","Bb","bb"), prob=prop_SNP2, replace=TRUE)
  SNP1=rep(0, n); SNP1[snp1=="Aa"]=1; SNP1[snp1=="aa"]=2 # additive 0=AA; 1=Aa; 2=aa
  SNP2=rep(0, n); SNP2[snp2=="Bb"]=1; SNP2[snp2=="bb"]=2 # additive  0=BB; 1=Bb; 2=bb
  SNP1a=rep(NA, n); SNP1a[snp1=="AA"]=1; SNP1a[snp1=="Aa"]=0; SNP1a[snp1=="aa"]=0 # additive 0=AA; 1=Aa; 2=aa
  SNP1b=rep(NA, n); SNP1b[snp1=="AA"]=0; SNP1b[snp1=="Aa"]=1; SNP1b[snp1=="aa"]=0 # additive 0=AA; 1=Aa; 2=aa
  SNP1c=rep(NA, n); SNP1c[snp1=="AA"]=0; SNP1c[snp1=="Aa"]=0; SNP1c[snp1=="aa"]=1 # additive 0=AA; 1=Aa; 2=aa
  SNP2a=rep(NA, n); SNP2a[snp2=="BB"]=1; SNP2a[snp2=="Bb"]=0; SNP2a[snp2=="bb"]=0 # additive 0=AA; 1=Aa; 2=aa
  SNP2b=rep(NA, n); SNP2b[snp2=="BB"]=0; SNP2b[snp2=="Bb"]=1; SNP2b[snp2=="bb"]=0 # additive 0=AA; 1=Aa; 2=aa
  SNP2c=rep(NA, n); SNP2c[snp2=="BB"]=0; SNP2c[snp2=="Bb"]=0; SNP2c[snp2=="bb"]=1 # additive 0=AA; 1=Aa; 2=aa
  out <- data.frame(snp1=snp1, SNP1=SNP1, snp2=snp2, SNP2=SNP2, SNP1a=SNP1a, SNP1b=SNP1b, SNP1c=SNP1c, SNP2a=SNP2a, SNP2b=SNP2b, SNP2c=SNP2c)
}

dense_rank <-   # copied from pkg::dplyr
  function (x) 
  {   r <- rank(x)
  match(r, sort(unique(r)))
  }
SimData<-function(n,alpha0=alpha0,beta0=beta0,gamma0=gamma0,tau=tau,maf=maf,it=i){
  set.seed(as.numeric(Sys.time())*((1/n)/it))
    obs.time = numeric(n)
    delta = numeric(n)
    z2=c(0,0)
    z3=c(0,0)
    delta=c(1,0)
    while(sum(z2[delta==1])==0|sum(z2[delta==0])==0|sum(z3[delta==1])==0|sum(z3[delta==0])==0){
    while(sum(z2)==0|sum(z3)==0){
    SNPdat <- SNPgen(n, prop_SNP1=c((1-maf-(maf/2.5)),maf,maf/2.5), prop_SNP2=c((1-maf-(maf/2.5)),maf,maf/2.5)) 
    z1 = rbinom(n,1,0.5)
    z2 = SNPdat[,2]
    z3 = SNPdat[,4]
    }
    delta=numeric(n)
    gen = SNPdat[,5:10]
    z=cbind(z1,z2,z3)
    cen.time = tau*runif(n)
    #exterme value error distribution
    temp = runif(n)
    err = log(-log(1-temp))
    death.time = exp(alpha0+beta0[1]*z[,1]+beta0[2]*z[,2]+beta0[3]*z[,3]+0.5*err)
    #cure probability and indicator
    cure.prob = exp(gamma0[1]+gamma0[2]*z[,1]+gamma0[3]*z[,2]+gamma0[4]*z[,3])/(1+exp(gamma0[1]+gamma0[2]*z[,1]+gamma0[3]*z[,2]+gamma0[4]*z[,3]))
    eta = rbinom(n,1,cure.prob)  
    #observed event time and censoring indicator
    obs.time[eta==0] = cen.time[eta==0]
    obs.time[eta==1] = apply(cbind(death.time[eta==1],cen.time[eta==1]),1,min)
    delta[eta==1] = as.numeric(death.time[eta==1]<=cen.time[eta==1])}
    return(cbind(z,delta,obs.time,gen))
}

LuMixSing<-function(obs.time,delta,z,clin,s,it=i){
  beta=ifelse(clin[3,1]<0.05,clin[1,1],0)
  gamma=c(0,0)
  gamma[1]=ifelse(clin[3,2]<0.05,clin[1,2],0)
  gamma[2]=ifelse(clin[3,3]<0.05,clin[1,3],0)
  gamma=as.matrix(gamma)
  
  x = cbind(rep(1,n),z[,1])
  #gamma[pvalue[2:3]>0.05]<-0
  neta = pb.f(n,gamma,z[,1])
  w = as.vector(eta_prob(gamma,z[,1],delta,s))
  sus.dev= w-neta
  
  #beta[pvalue[1]>0.05]<-0
  #rc=-log(s+0.0001)
  #mart2=delta-w*rc
  #surv.dev=sign(mart2)*sqrt(-2*(mart2+delta*log(delta-mart2)))
  e = log(obs.time)-beta*z[,1]
  
  sample1<-c("ID_1","ID_2","missing","sus.dev","e","s")
  sample2<-c(0,0,0,"P","P","P")
  samplel<-cbind(seq(1,n),seq(1,n),rep(0,n),sus.dev,e,s)
  sample<-rbind(sample1,sample2,samplel)
  write.table(sample,file = paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".sample",sep=""),sep=" ",quote=F,row.names = F,col.names = F)
  
  gen1a<-c("SNP1","rs1",1000,"A","C")
  gen1<-apply(z[,4:6], 1, paste, collapse=" ")
  gen2a<-c("SNP2","rs2",2000,"A","C")
  gen2<-apply(z[,7:9], 1, paste, collapse=" ")
  gen<-rbind(c(gen1a,gen1),c(gen2a,gen2))
  write.table(gen,file = paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".gen",sep=""),sep=" ",quote=F,row.names = F,col.names = F)
  system(paste("/work/projects/epipgx/users/LIV/Task1Surv/PLEIOTROPY --print_all --betas --pheno_name s --pheno_name sus.dev -o /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,"_1 -g /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".gen -s /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".sample",sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
  #system(paste("/work/projects/epipgx/users/LIV/Task1Surv/PLEIOTROPY --print_all --betas --pheno_name e --pheno_name sus.dev -o /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,"_2 -g /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".gen -s /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".sample",sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
  pleio1=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,"_1.result",sep=""),header=T)
  betas1=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,"_1.betas",sep=""),header=T)
  #pleio2=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,"_2.result",sep=""),header=T)
  #betas2=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,"_2.betas",sep=""),header=T)
  
  #system(paste("/work/projects/epipgx/tools/packages/snptest_v2.5_linux_x86_64_static/snptest_v2.5 -data /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".gen /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".sample -frequentist 1 -use_raw_phenotypes -o /work/projects/epipgx/users/LIV/Task1Surv/Simulation/SimulationSNPs",it," -pheno s -method expected",sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
  #system(paste("/work/projects/epipgx/tools/packages/snptest_v2.5_linux_x86_64_static/snptest_v2.5 -data /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".gen /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".sample -frequentist 1 -use_raw_phenotypes  -o /work/projects/epipgx/users/LIV/Task1Surv/Simulation/SimulationSNPsd",it," -pheno sus.dev -method expected",sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
  #system(paste("/work/projects/epipgx/tools/packages/snptest_v2.5_linux_x86_64_static/snptest_v2.5 -data /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".gen /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,".sample -frequentist 1 -use_raw_phenotypes  -o /work/projects/epipgx/users/LIV/Task1Surv/Simulation/SimulationSNPe",it," -pheno e -method expected",sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
  #SNPtest1=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SimulationSNPs",it,"",sep=""),header=T)
  #SNPtest2=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SimulationSNPsd",it,"",sep=""),header=T)
  #SNPtest3=read.table(paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SimulationSNPe",it,"",sep=""),header=T)
  
  #sus.SNP1.fit =  c(SNPtest2[1,23],SNPtest2[1,24],SNPtest2[1,21])
  #sus.SNP2.fit =  c(SNPtest2[2,23],SNPtest2[2,24],SNPtest2[2,21])
  #surv.SNP1.fit1 = c(SNPtest1[1,23],SNPtest1[1,24],SNPtest1[1,21])
  #surv.SNP2.fit1 = c(SNPtest1[2,23],SNPtest1[2,24],SNPtest1[2,21])
  #surv.SNP1.fit2 = c(SNPtest3[1,23],SNPtest3[1,24],SNPtest3[1,21])
  #surv.SNP2.fit2 = c(SNPtest3[2,23],SNPtest3[2,24],SNPtest3[2,21])
  #surv.SNP1.fit3 = c(betas1[2,7],betas1[2,8],pleio1[1,18])
  surv.SNP2.fit3 = c(betas1[6,7],betas1[6,8],pleio1[4,18])
  #sus.SNP1.fit3 =  c(betas1[1,7],betas1[1,8],pleio1[1,18])
  sus.SNP2.fit3 =  c(betas1[5,7],betas1[5,8],pleio1[4,18])
  #surv.SNP1.fit4 = c(betas2[2,7],betas2[2,8],pleio2[1,18])
  #surv.SNP2.fit4 = c(betas2[6,7],betas2[6,8],pleio2[4,18])
  #sus.SNP1.fit4 =  c(betas2[1,7],betas2[1,8],pleio2[1,18])
  #sus.SNP2.fit4 =  c(betas2[5,7],betas2[5,8],pleio2[4,18])
  
  surv.over.fit3 = c(pleio1[6,15],pleio1[6,16],pleio1[6,18])
  sus.over.fit3 =  c(pleio1[5,15],pleio1[5,16],pleio1[5,18])
  both.over.fit3 =  c(pleio1[4,15],pleio1[4,16],pleio1[4,18])
  
  #surv.over.fit4 = c(pleio2[6,15],pleio2[6,16],pleio2[6,18])
  #sus.over.fit4 =  c(pleio2[5,15],pleio2[5,16],pleio2[5,18])
  #both.over.fit4 =  c(pleio2[4,15],pleio2[4,16],pleio2[4,18])
  
  #combo1=1-pchisq(c(betas1[1,7],betas1[2,7])%*%c(betas1[1,8],betas1[2,8])),df=1)
  #combo2=1-pchisq((betas1[5,7]*betas1[6,7])/(betas1[5,8]*betas1[6,8]),df=1)
  
  clinpara=clin
  #parau1=cbind(clinpara,surv.SNP1.fit1,surv.SNP2.fit1,sus.SNP1.fit,sus.SNP2.fit)
  #parau2=cbind(clinpara,surv.SNP1.fit2,surv.SNP2.fit2,sus.SNP1.fit,sus.SNP2.fit)
  parau3=cbind(clinpara,surv.SNP2.fit3,sus.SNP2.fit3)#
  likeu3=cbind(surv.over.fit3,sus.over.fit3,both.over.fit3)
  #parau4=cbind(clinpara,surv.SNP1.fit4,surv.SNP2.fit4,sus.SNP1.fit4,sus.SNP2.fit4)
  #likeu4=cbind(surv.over.fit4,sus.over.fit4,both.over.fit4)
  
  parau=parau3#rbind(parau3,parau4)
  para=parau[,c(1,4,2,3,5)]
  like=likeu3#rbind(likeu3,likeu4)
  system(paste("rm -f /work/projects/epipgx/users/LIV/Task1Surv/Simulation/Simulation",it,"*",sep=""))
  ParaList<-list("para"=para,"like"=like)
  return(ParaList)
} 

LuMixFullPrev<-function(obs.time,delta,z){
  n=dim(z)[1]  
  z=as.matrix(z)
  lm.fit = lm(log(obs.time[delta==1]) ~ z[delta==1,])
  beta.ini = lm.fit$coef[-1]
  g.fit = glm(delta ~ z, family = "binomial")
  gamma.ini = g.fit$coef
  e = log(obs.time)-z%*%beta.ini 
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s.ini = l$surv
  s.ini = s_ext(s.ini,e,delta)
  
  sigma = sqrt(var(log(obs.time[delta==1])-z[delta==1,]%*%beta.ini))
  h.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
    beta = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-z%*%beta
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s = s_ext(s.fit,e,delta)
    
    dx = max(c(abs(beta-beta.ini),abs(gamma-gamma.ini)))
    gamma.ini = gamma
    beta.ini = beta
    s.ini = s
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  
 gamma.ini2=gamma[1:3] 
 beta.ini2=beta[1:2] 
 dx = 1.0
 iter = 1
 max.iter=1000  
 
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w2 = as.vector(eta_prob(gamma.ini2,z[,1:2],delta,s.ini))
    
    ####M-strep
    gamma2 = optim(gamma.ini2,post.lik1,w=w2,z=z[,1:2],control=list(reltol=0.00001,maxit=1000))$par
    beta2 = optim(beta.ini2,post.lik2,obs.time=obs.time,z=z[,1:2],delta=delta,w=w2,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e2 = log(obs.time)-z[,1:2]%*%beta2
    index = seq(1,n,1)
    index1 = index[delta[order(e2)]==1]
    Lambda2 = Lam.f(e2,delta,w2,h.opt)
    s.fit2 = exp(-Lambda2)
    s2 = s_ext(s.fit2,e2,delta)
    
    dx = max(c(abs(beta2-beta.ini2),abs(gamma2-gamma.ini2)))
    gamma.ini2 = gamma2
    beta.ini2 = beta2
    s.ini2 = s2
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  
 gamma.ini3=gamma[1:3] 
 beta.ini3=beta 
 dx = 1.0
 iter = 1
 max.iter=1000  
 
 while(iter <= max.iter && dx >= 0.01){
   ####E-step
   w3 = as.vector(eta_prob(gamma.ini3,z[,1:2],delta,s.ini))
   
   ####M-strep
   gamma3 = optim(gamma.ini3,post.lik1,w=w3,z=z[,1:2],control=list(reltol=0.00001,maxit=1000))$par
   beta3 = optim(beta.ini3,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w3,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
   
   e3 = log(obs.time)-z%*%beta3
   index = seq(1,n,1)
   index1 = index[delta[order(e3)]==1]
   Lambda3 = Lam.f(e3,delta,w3,h.opt)
   s.fit3 = exp(-Lambda3)
   s3 = s_ext(s.fit3,e3,delta)
   
   dx = max(c(abs(beta3-beta.ini3),abs(gamma3-gamma.ini3)))
   gamma.ini3 = gamma3
   beta.ini3 = beta3
   s.ini3 = s3
   iter = iter + 1  
   #cat(iter,dx,beta,gamma,"\n")
 }
 
 gamma.ini4=gamma
 beta.ini4=beta[1:2] 
 dx = 1.0
 iter = 1
 max.iter=1000  
 
 while(iter <= max.iter && dx >= 0.01){
   ####E-step
   w4 = as.vector(eta_prob(gamma.ini4,z,delta,s.ini))
   
   ####M-strep
   gamma4 = optim(gamma.ini4,post.lik1,w=w4,z=z,control=list(reltol=0.00001,maxit=1000))$par
   beta4 = optim(beta.ini4,post.lik2,obs.time=obs.time,z=z[,1:2],delta=delta,w=w4,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
   
   e4 = log(obs.time)-z[,1:2]%*%beta4
   index = seq(1,n,1)
   index1 = index[delta[order(e4)]==1]
   Lambda4 = Lam.f(e4,delta,w4,h.opt)
   s.fit4 = exp(-Lambda4)
   s4 = s_ext(s.fit4,e4,delta)
   
   dx = max(c(abs(beta4-beta.ini4),abs(gamma4-gamma.ini4)))
   gamma.ini4 = gamma4
   beta.ini4 = beta4
   s.ini4 = s4
   iter = iter + 1  
   #cat(iter,dx,beta,gamma,"\n")
 }
 
  ########variance estimation based on profile likelihood 
  ###fix beta
  h.n = 2/n
  pl.s1=sapply(1:length(beta),function(i){
    beta.ini[i] = beta[i] + h.n
    beta.ini[-i] = beta[-i] 
    gamma.ini = gamma
    beta.temp = beta.ini[-i]
    beta.new = numeric(length(beta))
    beta.new[i] = beta.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.01){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
      ####M-strep
      gamma.new = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
      beta.new[-i] = optim(beta.temp,post.lik2n,beta2=beta.ini[i],k=i,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
      e = log(obs.time)-z%*%beta.ini
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta)
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      beta.ini = beta.new
      beta.temp = beta.new[-i]
      s.ini = s.new
      iter = iter + 1   
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }  
    w = as.vector(eta_prob(gamma.new,z,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)
  })
  
  pl.s2=sapply(1:length(beta),function(i){
    beta.ini[i] = beta[i] - h.n
    beta.ini[-i] = beta[-i] 
    gamma.ini = gamma
    beta.temp = beta.ini[-i]
    beta.new = numeric(length(beta))
    beta.new[i] = beta.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.01){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
      ####M-strep
      gamma.new = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
      beta.new[-i] = optim(beta.temp,post.lik2n,beta2=beta.ini[i],k=i,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
      e = log(obs.time)-z%*%beta.ini
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta)
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      beta.ini = beta.new
      beta.temp = beta.new[-i]
      s.ini = s.new
      iter = iter + 1   
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }  
    w = as.vector(eta_prob(gamma.new,z,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)
  })
  
  pl.s3=sapply(1:length(gamma),function(i){
    beta.ini = beta
    gamma.ini[i] = gamma[i] + h.n
    gamma.ini[-i] = gamma[-i]
    gamma.temp = gamma.ini[-i]
    gamma.new = numeric(length(gamma))
    gamma.new[i] = gamma.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.012){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
      
      ####M-strep
      gamma.new[-i] = optim(gamma.temp,post.lik1n,gamma2=gamma.ini[i],k=i,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
      beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
      
      e = log(obs.time)-z%*%beta.ini
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta) 
      
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      beta.ini = beta.new
      gamma.temp = gamma.new[-i]
      s.ini = s.new
      iter = iter + 1  
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }
    w = as.vector(eta_prob(gamma.new,z,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)})
  
  pl.s4=sapply(1:length(gamma),function(i){
    beta.ini = beta
    gamma.ini[i] = gamma[i] - h.n
    gamma.ini[-i] = gamma[-i]
    gamma.temp = gamma.ini[-i]
    gamma.new = numeric(length(gamma))
    gamma.new[i] = gamma.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.012){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
      
      ####M-strep
      gamma.new[-i] = optim(gamma.temp,post.lik1n,gamma2=gamma.ini[i],k=i,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
      beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
      
      e = log(obs.time)-z%*%beta.ini
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta)
      
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      beta.ini = beta.new
      gamma.temp = gamma.new[-i]
      s.ini = s.new
      iter = iter + 1  
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }
    w = as.vector(eta_prob(gamma.new,z,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)})

  smat = cbind(sapply(1:length(beta),function(i){((pl.s1[,i]-pl.s2[,i])/(2*h.n))}),sapply(1:length(gamma),function(i){((pl.s3[,i]-pl.s4[,i])/(2*h.n))}))
  est.var = solve(t(smat)%*%(smat))
  err=as.matrix(sqrt(diag(est.var)))
  coef=as.matrix(c(beta,gamma))
  pvalue=sapply(1:length(coef),function(i){2*(1-pnorm(abs(coef[i]/err[i])))})
  w = as.vector(eta_prob(gamma,z,delta,s))
  w3 = as.vector(eta_prob(gamma3,z[,1:2],delta,s3))
  w4 = as.vector(eta_prob(gamma4,z,delta,s4))
  w2 = as.vector(eta_prob(gamma2,z[,1:2],delta,s2))

  alt.l=sum(pl.lik(beta,gamma,obs.time,z,delta,w,h.opt))
  susnull.l=sum(pl.lik(beta3,gamma3,obs.time,z,delta,w3,h.opt))
  survnull.l=sum(pl.lik(beta4,gamma4,obs.time,z,delta,w4,h.opt))
  null.l=sum(pl.lik(beta2,gamma2,obs.time,z,delta,w2,h.opt))
  
  D.sus=2*(alt.l-susnull.l)
  D.surv=2*(alt.l-survnull.l)
  D.both=2*(alt.l-null.l)
  
  p.sus=pchisq(D.sus,df=1,lower.tail = FALSE)
  p.surv=pchisq(D.surv,df=1,lower.tail = FALSE)
  p.both=pchisq(D.both,df=2,lower.tail = FALSE)

  return(cbind(rbind(t(coef),t(err),t(pvalue)),c(alt.l,survnull.l,p.surv),c(alt.l,susnull.l,p.sus),c(alt.l,null.l,p.both)))#
}
  

ahmc.data<-function(b,n,a,dist,link){
  #b<-c(2,-1,1);n<-500;dist<-"lognormal";a<-40
  z<-array(0,c(n,2))
  y<-x<-censor<-Status<-Time<-array(0,n)
  time<-array(Inf,n)
  for (i in 1:n){
    while (time[i]==Inf){
      x[i]<-rnorm(1,0,1);z[i,]<-c(1,rbinom(1,1,0.5))
      #x[i]<-rbinom(1,1,0.5);z[i,]<-c(1,x[i])
      temp1<-z[i,]%*%b[-3]
      if (link=="logit") {p1<-exp(temp1)/(1+exp(temp1))}
      if (link=="probit") {p1<-pnorm(temp1)}
      if (link=="cloglog") {p1<-1-exp(-exp(temp1))}
      y[i]<-rbinom(1,1,p1)
      temp2<-exp(x[i]*b[3])
      temp3<-1-exp(log(1-runif(1,0,1))*temp2)
      time[i]<-qlnorm(temp3,0,1)/temp2
    }}
  censor<-runif(n,0,a)
  Status<-as.numeric(time<censor)
  for (i in 1:n){
    if (y[i]==0){Status[i]<-0;Time[i]<-censor[i]}
    else {Time[i]<-Status[i]*time[i]+(1-Status[i])*censor[i]}}
  censor_rate<-1-sum(Status)/n
  list(Time=Time,Status=Status,x=x,z=z,crate=censor_rate)
}

ahmc.apply<-function(simu){
  source("/home/ben/Dropbox/Epilepsy/Simulation/ahmc.data.R")
  source("/home/ben/Dropbox/Epilepsy/Simulation/ahmc.profile.em.fun.R") 
  tau0<-1e-4
  
  finish<-FALSE
  while (!finish){
    
    #generate a simulation data set
    true<-c(2,-1,1)
    truelink<-"logit"
    data<-ahmc.data(true,n,a,dist,truelink)
    Time<-data$Time;Status<-data$Status;x<-data$x;z<-data$z;n<-length(Time)
    
    plink1<-ahmc.profile.em.fun(Time,Status,x,z,tau0,"logit")
    if (is.character(plink1)|plink1$tau>tau0) {finish<-FALSE;next} 
    cat(simu,"profile=",formatC(c(plink1$est,plink1$cure),digits=3,format="f"),"\n",file=estimatef1,append=TRUE)
    
    plink2<-ahmc.profile.em.fun(Time,Status,x,z,tau0,"probit")
    if (is.character(plink2)|plink2$tau>tau0) {finish<-FALSE;next}  
    cat(simu,"profile=",formatC(c(plink2$est,plink2$cure),digits=3,format="f"),"\n",file=estimatef2,append=TRUE)
    
    plink3<-ahmc.profile.em.fun(Time,Status,x,z,tau0,"cloglog")
    if (is.character(plink3)|plink3$tau>tau0) {finish<-FALSE;next}  
    cat(simu,"profile=",formatC(c(plink3$est,plink3$cure),digits=3,format="f"),"\n",file=estimatef3,append=TRUE)
    
    finish<-TRUE
  }
  list(plink1=plink1,plink2=plink2,plink3=plink3)
}

ahmc.profile.em.fun<-function(Time,Status,x,z,tau0,link){
  
  #estimate survival function
  survival<-function(beta,w){
    ltime<-Time*exp(x*beta)
    e<-exp(-x*beta)
    
    h0_profile<-function(t){
      temp1<-temp2<-array(0,length(t))
      for (i in 1:length(t)){
        res1<-(log(ltime)-log(t[i]))/an1
        res2<-(log(ltime)-log(t[i]))/an2
        temp1[i]<-sum(Status*dnorm(res1)/(an1*t[i]))
        temp2[i]<-sum(w*pnorm(res2)*e)}
      temp1/temp2} 
    
    survival<-array(0,n)
    t.max<-max(subset(cbind(ltime,Status),Status==1)[,1])
    #out<-gauss.quad(20,"legendre")
    for (i in 1:n){
      hazard<-integrate(h0_profile,lower=0,ltime[i],stop.on.error=FALSE)$value
      #hazard<-sum(h0_profile((out$nodes*ltime[i]+ltime[i])/2)*ltime[i]/2*out$weight)
      survival[i]<-ifelse(ltime[i]>t.max,0,exp(-hazard))
    }
    list(survival=survival,h0=h0_profile(ltime))}
  
  #initial value
  ##added b coeff instead of "true"
  true<-c(12,20,3)
  b<-true;w<-Status;s<-array(0,n)
  update_b<-array(0,length(true))
  #bandwidths
  sigma<-sqrt(var(log(Time)))
  an1<-an2<-sigma*n^(-1/4)
  
  lc1<-function(b,w){
    temp<-c(z%*%b)
    if (link=="logit") p1<-exp(temp)/(1+exp(temp))
    if (link=="probit") p1<-pnorm(temp)
    if (link=="cloglog") p1<-1-exp(-exp(temp))
    p0<-1-p1
    -sum(w*log(p1)+(1-w)*log(p0))/n}
  
  lc2<-function(beta,w){
    ltime<-log(Time)+x*beta;sum1<-array(0,n)
    for (i in 1:n){
      res1<-(ltime-ltime[i])/an1
      res2<-(ltime-ltime[i])/an2
      sum1[i]<-log((sum(Status*dnorm(res1))/an1)/(sum(w*exp(-x*beta)*pnorm(res2))))}  
    -sum((-ltime+sum1)*Status)/n}
  
  tau<-1;count<-1
  while (tau>tau0 & count<1000){ 
    #update_b[-3]<-as.numeric(glm(w~z[,-1],family=quasibinomial(link=link))$coef) 
    #lc1_fit<-try(optim(b[-3],lc1,w=w,method="Nelder-Mead"))
    #if (is.character(lc1_fit)) {tau<-1;break}
    #update_b[-3]<-lc1_fit$par
    #lc2_fit<-try(optim(b[3],lc2,w=w,method="Nelder-Mead"))
    #if (is.character(lc2_fit)) {tau<-1;break}
    #update_b[3]<-lc2_fit$par
    lc<-function(b){lc1(b[-3],w)+lc2(b[3],w)}
    lc_fit<-try(optim(b,lc,method="Nelder-Mead"))
    if (is.character(lc_fit)) {tau<-1;break}
    update_b<-lc_fit$par
    tau1<-sum((update_b-b)^2);b<-update_b
    
    survival_fit<-try(survival(b[3],w))
    if (is.character(survival_fit)) {tau<-1;break}
    update_s<-survival_fit$survival^(exp(-x*b[3]))
    tau2<-mean((update_s-s)^2);s<-update_s 
    
    temp<-c(z%*%b[-3])
    if (link=="logit") p1<-exp(temp)/(1+exp(temp))
    if (link=="probit") p1<-pnorm(temp)
    if (link=="cloglog") p1<-1-exp(-exp(temp))
    p0<-1-p1
    w<-((s*p1)/(p0+p1*s))*(1-Status)+Status 
    
    tau<-max(tau1);count<-count+1}
  
  
  cure<-c(0,0)
  if (link=="logit") {cure[1]<-exp(b[1])/(1+exp(b[1]));cure[2]<-exp(b[1]+b[2])/(1+exp(b[1]+b[2]))}
  if (link=="probit") {cure[1]<-pnorm(b[1]);cure[2]<-pnorm(b[1]+b[2])}
  if (link=="cloglog") {cure[1]<-1-exp(-exp(b[1]));cure[2]<-1-exp(-exp(b[1]+b[2]))}
  cure<-(c(1,1)-cure)*100
  
  list(est=b,cure=cure,tau=tau,count=count)} 

LuMixFix<-function(obs.time,delta,z){
  n=dim(z)[1]  
  lm.fit = lm(log(obs.time[delta==1])~z[delta==1,1])
  beta.ini = lm.fit$coef[2]
  g.fit = glm(delta ~ z[,1], family = "binomial")
  gamma.ini = g.fit$coef
  e = log(obs.time)-beta.ini*z[,1]  
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s.ini = l$surv
  s.ini = s_ext(s.ini,e,delta)
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta.ini*z[delta==1,1]))
  h.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.ini))
    
    ####M-strep
    gamma = optim(gamma.ini,post.lik1,w=w,z=z[,1],control=list(reltol=0.00001,maxit=1000))$par
    beta = optim(beta.ini,post.lik2,obs.time=obs.time,z=z[,1],delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta*z[,1]
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s = s_ext(s.fit,e,delta)
    
    dx = max(c(abs(beta-beta.ini),abs(gamma-gamma.ini)))
    gamma.ini = gamma
    beta.ini = beta
    s.ini = s
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  
  ########variance estimation based on profile likelihood 
  ###fix beta
  h.n = 2/n
  
  beta.ini = beta + h.n
  gamma.ini = gamma
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.ini))
    
    ####M-strep
    gamma.new = optim(gamma.ini,post.lik1,w=w,z=z[,1],control=list(reltol=0.00001,maxit=1000))$par
    
    e = log(obs.time)-beta.ini*z[,1]
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(gamma.new-gamma.ini))
    gamma.ini = gamma.new
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(eta_prob(gamma.new,z[,1],delta,s.new))
  pl.s1 = pl.lik(beta.ini,gamma.new,obs.time,z[,1],delta,w,h.opt)
  
  beta.ini = beta - h.n
  gamma.ini = gamma
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.ini))
    
    ####M-strep
    gamma.new = optim(gamma.ini,post.lik1,w=w,z=z[,1],control=list(reltol=0.00001,maxit=1000))$par
    
    e = log(obs.time)-beta.ini*z[,1]
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(gamma.new-gamma.ini))
    gamma.ini = gamma.new
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(eta_prob(gamma.new,z[,1],delta,s.new))
  pl.s2 = pl.lik(beta.ini,gamma.new,obs.time,z[,1],delta,w,h.opt)
  
  ##fix gamma[1]
  beta.ini = beta
  gamma.ini[1] = gamma[1] + h.n
  gamma.ini[2] = gamma[2]
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.ini))
    
    ####M-strep
    gamma.new[2] = optim(gamma.ini[2],post.lik12,gamma1=gamma.ini[1],w=w,z=z[,1],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z[,1],delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z[,1]
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[2]-gamma.ini[2]))
    beta.ini = beta.new
    gamma.ini[2] = gamma.new[2]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.new))
  pl.s3 = pl.lik(beta.new,gamma.ini,obs.time,z[,1],delta,w,h.opt)
  
  beta.ini = beta
  gamma.ini[1] = gamma[1] - h.n
  gamma.ini[2] = gamma[2]
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.ini))
    
    ####M-strep
    gamma.new[2] = optim(gamma.ini[2],post.lik12,gamma1=gamma.ini[1],w=w,z=z[,1],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z[,1],delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z[,1]
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[2]-gamma.ini[2]))
    beta.ini = beta.new
    gamma.ini[2] = gamma.new[2]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.new))
  pl.s4 = pl.lik(beta.new,gamma.ini,obs.time,z[,1],delta,w,h.opt)
  
  ##fix gamma[2]
  beta.ini = beta
  gamma.ini[1] = gamma[1]
  gamma.ini[2] = gamma[2] + h.n
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.ini))
    
    ####M-strep
    gamma.new[1] = optim(gamma.ini[1],post.lik11,gamma2=gamma.ini[2],w=w,z=z[,1],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z[,1],delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z[,1]
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[1]-gamma.ini[1]))
    beta.ini = beta.new
    gamma.ini[1] = gamma.new[1]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.new))
  pl.s5 = pl.lik(beta.new,gamma.ini,obs.time,z[,1],delta,w,h.opt)
  
  beta.ini = beta
  gamma.ini[1] = gamma[1]
  gamma.ini[2] = gamma[2] - h.n
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.ini))
    
    ####M-strep
    gamma.new[1] = optim(gamma.ini[1],post.lik11,gamma2=gamma.ini[2],w=w,z=z[,1],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z[,1],delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z[,1]
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[1]-gamma.ini[1]))
    beta.ini = beta.new
    gamma.ini[1] = gamma.new[1]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z[,1],delta,s.new))
  pl.s6 = pl.lik(beta.new,gamma.ini,obs.time,z[,1],delta,w,h.opt)
  
  smat = cbind(cbind((pl.s1-pl.s2)/(2*h.n),(pl.s3-pl.s4)/(2*h.n)),(pl.s5-pl.s6)/(2*h.n))
  est.var = solve(t(smat)%*%(smat))
  err=as.matrix(sqrt(diag(est.var)))
  coef=as.matrix(c(beta,gamma))
  pvalue=sapply(1:length(coef),function(i){2*(1-pnorm(abs(coef[i]/err[i])))})
  
  x = cbind(rep(1,n),z[,1])
  gamma[pvalue[2:3]>0.05]<-0
  gamma.clin=x%*%gamma
  
  beta[pvalue[1]>0.05]<-0
  beta.clin=beta*z[,1]
  
  #####Use fixed values in another full model######

  
  n=dim(z)[1]  
  lm.fit = lm(log(obs.time[delta==1])~z[delta==1,2])
  beta1.ini = lm.fit$coef[2]
  g.fit = glm(delta ~ z[,2], family = "binomial")
  gamma1.ini = g.fit$coef[2]
  e = log(obs.time)-beta1.ini*z[,2]-beta.clin
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s1.ini = l$surv
  s1.ini = s_ext(s1.ini,e,delta)
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta1.ini*z[delta==1,2]-beta.clin[delta==1]))
  h1.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z[,2],delta,s1.ini))
    
    ####M-strep
    gamma1 = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z[,2],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta1 = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z[,2],delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1*z[,2]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1 = s_ext(s1.fit,e,delta)
    
    dx = max(c(abs(beta1-beta1.ini),abs(gamma1-gamma1.ini)))
    gamma1.ini = gamma1
    beta1.ini = beta1
    s1.ini = s1
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  
  ########variance estimation based on profile likelihood 
  ###fix beta
  h1.n = 2/n
  
  beta1.ini = beta1 + h1.n
  gamma1.ini = gamma1
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z[,2],delta,s1.ini))
    
    ####M-strep
    gamma1.new = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z[,2],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.ini*z[,2]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta)
    
    dx = max(abs(gamma1.new-gamma1.ini))
    gamma1.ini = gamma1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(feta_prob(gamma1.new,gamma.clin,z[,2],delta,s1.new))
  f1pl.s1 = fpl.lik(beta1.ini,gamma1.new,beta.clin,gamma.clin,obs.time,z[,2],delta,w,h1.opt)
  
  beta1.ini = beta1 - h.n
  gamma1.ini = gamma1
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z[,2],delta,s1.ini))
    
    ####M-strep
    gamma1.new = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z[,2],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.ini*z[,2]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta) 
    
    dx = max(abs(gamma1.new-gamma1.ini))
    gamma1.ini = gamma1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(feta_prob(gamma1.new,gamma.clin,z[,2],delta,s1.new))
  f1pl.s2 = fpl.lik(beta1.ini,gamma1.new,beta.clin,gamma.clin,obs.time,z[,2],delta,w,h1.opt)
  
  ##fix gamma[1]
  beta1.ini = beta1
  gamma1.ini = gamma1 + h.n
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z[,2],delta,s1.ini))
    
    ####M-strep
    beta1.new = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z[,2],delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.new*z[,2]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta) 
    
    dx = max(abs(beta1.new-beta1.ini))
    beta1.ini = beta1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(feta_prob(gamma1.ini,gamma.clin,z[,2],delta,s1.new))
  f1pl.s3 = fpl.lik(beta1.new,gamma1.ini,beta.clin,gamma.clin,obs.time,z[,2],delta,w,h1.opt)
  
  beta1.ini = beta1
  gamma1.ini = gamma1 - h.n
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z[,2],delta,s1.ini))
    
    ####M-strep
    beta.new = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z[,2],delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.new*z[,2]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta) 
    
    dx = max(abs(beta1.new-beta1.ini))
    beta1.ini = beta1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(feta_prob(gamma1.ini,gamma.clin,z[,2],delta,s1.new))
  f1pl.s4 = fpl.lik(beta1.new,gamma1.ini,beta.clin,gamma.clin,obs.time,z[,2],delta,w,h1.opt)
  
  fsmat1 = cbind(cbind((f1pl.s1-f1pl.s2)/(2*h.n),(f1pl.s3-f1pl.s4)/(2*h.n)))
  est.var1 = solve(t(fsmat1)%*%(fsmat1))
  err1=as.matrix(sqrt(diag(est.var1)))
  coef1=as.matrix(c(beta1,gamma1))
  pvalue1=sapply(1:length(coef1),function(i){2*(1-pnorm(abs(coef1[i]/err1[i])))})
  
  #####Use fixed values in another full model######Go for 2!!!
  
  n=dim(z)[1]  
  lm.fit = lm(log(obs.time[delta==1])~z[delta==1,3])
  beta2.ini = lm.fit$coef[2]
  g.fit = glm(delta ~ z[,3], family = "binomial")
  gamma2.ini = g.fit$coef[2]
  e = log(obs.time)-beta2.ini*z[,3]  
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s2.ini = l$surv
  s2.ini = s_ext(s2.ini,e,delta) 
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta2.ini*z[delta==1,3]-beta.clin[delta==1]))
  h2.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma2.ini,gamma.clin,z[,3],delta,s2.ini))
    
    ####M-strep
    gamma2 = optim(gamma2.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z[,3],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta2 = optim(beta2.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z[,3],delta=delta,w=w,h=h2.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta2*z[,3]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h2.opt)
    s2.fit = exp(-Lambda)
    s2 = s_ext(s2.fit,e,delta) 
    
    dx = max(c(abs(beta2-beta2.ini),abs(gamma2-gamma2.ini)))
    gamma2.ini = gamma2
    beta2.ini = beta2
    s2.ini = s2
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  
  ########variance estimation based on profile likelihood 
  ###fix beta
  h2.n = 2/n
  
  beta2.ini = beta2 + h2.n
  gamma2.ini = gamma2
  s2.ini = s2
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma2.ini,gamma.clin,z[,3],delta,s2.ini))
    
    ####M-strep
    gamma2.new = optim(gamma2.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z[,3],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta2.ini*z[,3]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h2.opt)
    s2.fit = exp(-Lambda)
    s2.new = s_ext(s2.fit,e,delta) 
    
    dx = max(abs(gamma2.new-gamma2.ini))
    gamma2.ini = gamma2.new
    s2.ini = s2.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(feta_prob(gamma2.new,gamma.clin,z[,3],delta,s2.new))
  f2pl.s1 = fpl.lik(beta2.ini,gamma2.new,beta.clin,gamma.clin,obs.time,z[,3],delta,w,h2.opt)
  
  beta2.ini = beta2 - h.n
  gamma2.ini = gamma2
  s2.ini = s2
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma2.ini,gamma.clin,z[,3],delta,s2.ini))
    
    ####M-strep
    gamma2.new = optim(gamma2.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z[,3],control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta2.ini*z[,3]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s2.fit = exp(-Lambda)
    s2.new = s_ext(s2.fit,e,delta) 
    
    dx = max(abs(gamma2.new-gamma2.ini))
    gamma2.ini = gamma2.new
    s2.ini = s2.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(feta_prob(gamma2.new,gamma.clin,z[,3],delta,s2.new))
  f2pl.s2 = fpl.lik(beta2.ini,gamma2.new,beta.clin,gamma.clin,obs.time,z[,3],delta,w,h2.opt)
  
  ##fix gamma[1]
  beta2.ini = beta2
  gamma2.ini = gamma2 + h.n
  s2.ini = s2
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma2.ini,gamma.clin,z[,3],delta,s2.ini))
    
    ####M-strep
    beta2.new = optim(beta2.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z[,3],delta=delta,w=w,h=h2.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta2.new*z[,3]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s2.fit = exp(-Lambda)
    s2.new = s_ext(s2.fit,e,delta) 
    
    dx = max(abs(beta2.new-beta2.ini))
    beta2.ini = beta2.new
    s2.ini = s2.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(feta_prob(gamma2.ini,gamma.clin,z[,3],delta,s2.new))
  f2pl.s3 = fpl.lik(beta2.new,gamma2.ini,beta.clin,gamma.clin,obs.time,z[,3],delta,w,h2.opt)
  
  beta2.ini = beta2
  gamma2.ini = gamma2 - h.n
  s2.ini = s2
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma2.ini,gamma.clin,z[,3],delta,s2.ini))
    
    ####M-strep
    beta.new = optim(beta2.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z[,3],delta=delta,w=w,h=h2.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta2.new*z[,3]-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s2.fit = exp(-Lambda)
    s2.new = s_ext(s2.fit,e,delta) 
    
    dx = max(abs(beta2.new-beta2.ini))
    beta2.ini = beta2.new
    s2.ini = s2.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(feta_prob(gamma2.ini,gamma.clin,z[,3],delta,s2.new))
  f2pl.s4 = fpl.lik(beta2.new,gamma2.ini,beta.clin,gamma.clin,obs.time,z[,3],delta,w,h2.opt)
  
  fsmat2 = cbind(cbind((f2pl.s1-f2pl.s2)/(2*h.n),(f2pl.s3-f2pl.s4)/(2*h.n)))
  est.var2 = solve(t(fsmat2)%*%(fsmat2))
  err2=as.matrix(sqrt(diag(est.var2)))
  coef2=as.matrix(c(beta2,gamma2))
  pvalue2=sapply(1:length(coef2),function(i){2*(1-pnorm(abs(coef2[i]/err2[i])))})
  
  w = as.vector(eta_prob(gamma,z,delta,s))
  w3 = as.vector(eta_prob(gamma3,z[,1:2],delta,s3))
  w4 = as.vector(eta_prob(gamma4,z,delta,s4))
  w2 = as.vector(eta_prob(gamma2,z[,1:2],delta,s2))
  
  alt.l=sum(pl.lik(beta,gamma,obs.time,z,delta,w,h.opt))
  susnull.l=sum(pl.lik(beta3,gamma3,obs.time,z,delta,w3,h.opt))
  survnull.l=sum(pl.lik(beta4,gamma4,obs.time,z,delta,w4,h.opt))
  null.l=sum(pl.lik(beta2,gamma2,obs.time,z,delta,w2,h.opt))
  
  D.sus=2*(alt.l-susnull.l)
  D.surv=2*(alt.l-survnull.l)
  D.both=2*(alt.l-null.l)
  
  p.sus=pchisq(D.sus,df=1,lower.tail = FALSE)
  p.surv=pchisq(D.surv,df=1,lower.tail = FALSE)
  p.both=pchisq(D.both,df=2,lower.tail = FALSE)
  
  return(cbind(rbind(t(coef),t(err),t(pvalue)),c(alt.l,survnull.l,p.surv),c(alt.l,susnull.l,p.sus),c(alt.l,null.l,p.both)))
  
  
  clinpara=rbind(t(coef),t(err),t(pvalue))
  genpara1=rbind(t(coef1),t(err1),t(pvalue1))
  genpara2=rbind(t(coef2),t(err2),t(pvalue2))
  parau1=cbind(clinpara,genpara1,genpara2)
  para=parau1[,c(1,4,6,2,3,5,7)]
  return(para)
}

LuMixFixRe<-function(obs.time,delta,beta.clin,gamma.clin,z){
  n=length(z) 
  lm.fit = lm(log(obs.time)[delta==1]~z[delta==1])
  beta1.ini = lm.fit$coef[2]
  g.fit = glm(delta ~ z, family = "binomial")
  gamma1.ini = g.fit$coef[2]
  e = log(obs.time)-beta1.ini*z-beta.clin
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s1.ini = l$surv
  s1.ini = s_ext(s1.ini,e,delta) 
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta1.ini*z[delta==1]-beta.clin[delta==1]))
  h1.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.ini))
    
    ####M-strep
    gamma1 = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta1 = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z,delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1*z-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1 = s_ext(s1.fit,e,delta) 
    
    dx = max(c(abs(beta1-beta1.ini),abs(gamma1-gamma1.ini)))
    gamma1.ini = gamma1
    beta1.ini = beta1
    s1.ini = s1
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  
  ########variance estimation based on profile likelihood 
  ###fix beta
  h.n = 2/n
  
  beta1.ini = beta1 + h.n
  gamma1.ini = gamma1
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.ini))
    
    ####M-strep
    gamma1.new = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.ini*z-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta)  
    
    dx = max(abs(gamma1.new-gamma1.ini))
    gamma1.ini = gamma1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(feta_prob(gamma1.new,gamma.clin,z,delta,s1.new))
  f1pl.s1 = fpl.lik(beta1.ini,gamma1.new,beta.clin,gamma.clin,obs.time,z,delta,w,h1.opt)
  
  beta1.ini = beta1 - h.n
  gamma1.ini = gamma1
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.ini))
    
    ####M-strep
    gamma1.new = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.ini*z-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta) 
    
    dx = max(abs(gamma1.new-gamma1.ini))
    gamma1.ini = gamma1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(feta_prob(gamma1.new,gamma.clin,z,delta,s1.new))
  f1pl.s2 = fpl.lik(beta1.ini,gamma1.new,beta.clin,gamma.clin,obs.time,z,delta,w,h1.opt)
  
  ##fix gamma[1]
  beta1.ini = beta1
  gamma1.ini = gamma1 + h.n
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.ini))
    
    ####M-strep
    beta1.new = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z,delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.new*z-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta) 
    
    dx = max(abs(beta1.new-beta1.ini))
    beta1.ini = beta1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.new))
  f1pl.s3 = fpl.lik(beta1.new,gamma1.ini,beta.clin,gamma.clin,obs.time,z,delta,w,h1.opt)
  
  beta1.ini = beta1
  gamma1.ini = gamma1 - h.n
  s1.ini = s1
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.ini))
    
    ####M-strep
    beta.new = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z,delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1.new*z-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1.new = s_ext(s1.fit,e,delta) 
    
    dx = max(abs(beta1.new-beta1.ini))
    beta1.ini = beta1.new
    s1.ini = s1.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.new))
  f1pl.s4 = fpl.lik(beta1.new,gamma1.ini,beta.clin,gamma.clin,obs.time,z,delta,w,h1.opt)
  
  fsmat1 = cbind(cbind((f1pl.s1-f1pl.s2)/(2*h.n),(f1pl.s3-f1pl.s4)/(2*h.n)))
  est.var1 = solve(t(fsmat1)%*%(fsmat1))
  err1=as.matrix(sqrt(diag(est.var1)))
  coef1=as.matrix(c(beta1,gamma1))
  pvalue1=sapply(1:length(coef1),function(i){2*(1-pnorm(abs(coef1[i]/err1[i])))})
  
  genpara2=rbind(t(coef1),t(err1),t(pvalue1))
  genpara2
  return(genpara2)
}

LuMixGENOSing<-function(obs.time,delta,z){
  n=length(z)
  lm.fit = lm(log(obs.time[delta==1])~z[delta==1])
  beta.ini = lm.fit$coef[2]
  g.fit = glm(delta ~ z, family = "binomial")
  gamma.ini = g.fit$coef
  e = log(obs.time)-beta.ini*z  
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s.ini = l$surv
  s.ini = s_ext(s.ini,e,delta)
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta.ini*z[delta==1]))
  h.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
    beta = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta*z
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s = s_ext(s.fit,e,delta) 
    
    dx = max(c(abs(beta-beta.ini),abs(gamma-gamma.ini)))
    gamma.ini = gamma
    beta.ini = beta
    s.ini = s
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  
  ########variance estimation based on profile likelihood 
  ###fix beta
  h.n = 2/n
  
  beta.ini = beta + h.n
  gamma.ini = gamma
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma.new = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
    
    e = log(obs.time)-beta.ini*z
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(gamma.new-gamma.ini))
    gamma.ini = gamma.new
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(eta_prob(gamma.new,z,delta,s.new))
  pl.s1 = pl.lik(beta.ini,gamma.new,obs.time,z,delta,w,h.opt)
  
  beta.ini = beta - h.n
  gamma.ini = gamma
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma.new = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
    
    e = log(obs.time)-beta.ini*z
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(gamma.new-gamma.ini))
    gamma.ini = gamma.new
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.ini,gamma.new,"\n")
  }
  w = as.vector(eta_prob(gamma.new,z,delta,s.new))
  pl.s2 = pl.lik(beta.ini,gamma.new,obs.time,z,delta,w,h.opt)
  
  ##fix gamma[1]
  beta.ini = beta
  gamma.ini[1] = gamma[1] + h.n
  gamma.ini[2] = gamma[2]
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma.new[2] = optim(gamma.ini[2],post.lik12,gamma1=gamma.ini[1],w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[2]-gamma.ini[2]))
    beta.ini = beta.new
    gamma.ini[2] = gamma.new[2]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z,delta,s.new))
  pl.s3 = pl.lik(beta.new,gamma.ini,obs.time,z,delta,w,h.opt)
  
  beta.ini = beta
  gamma.ini[1] = gamma[1] - h.n
  gamma.ini[2] = gamma[2]
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma.new[2] = optim(gamma.ini[2],post.lik12,gamma1=gamma.ini[1],w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta) 
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[2]-gamma.ini[2]))
    beta.ini = beta.new
    gamma.ini[2] = gamma.new[2]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z,delta,s.new))
  pl.s4 = pl.lik(beta.new,gamma.ini,obs.time,z,delta,w,h.opt)
  
  ##fix gamma[2]
  beta.ini = beta
  gamma.ini[1] = gamma[1]
  gamma.ini[2] = gamma[2] + h.n
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma.new[1] = optim(gamma.ini[1],post.lik11,gamma2=gamma.ini[2],w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[1]-gamma.ini[1]))
    beta.ini = beta.new
    gamma.ini[1] = gamma.new[1]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z,delta,s.new))
  pl.s5 = pl.lik(beta.new,gamma.ini,obs.time,z,delta,w,h.opt)
  
  beta.ini = beta
  gamma.ini[1] = gamma[1]
  gamma.ini[2] = gamma[2] - h.n
  s.ini = s
  dx = 1.0
  iter = 1
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
    
    ####M-strep
    gamma.new[1] = optim(gamma.ini[1],post.lik11,gamma2=gamma.ini[2],w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.new*z
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s.new = s_ext(s.fit,e,delta)
    
    dx = max(abs(beta.new-beta.ini),abs(gamma.new[1]-gamma.ini[1]))
    beta.ini = beta.new
    gamma.ini[1] = gamma.new[1]
    s.ini = s.new
    iter = iter + 1  
    #cat(iter,dx,beta.new,gamma.ini,"\n")
  }
  w = as.vector(eta_prob(gamma.ini,z,delta,s.new))
  pl.s6 = pl.lik(beta.new,gamma.ini,obs.time,z,delta,w,h.opt)
  
  smat = cbind(cbind((pl.s1-pl.s2)/(2*h.n),(pl.s3-pl.s4)/(2*h.n)),(pl.s5-pl.s6)/(2*h.n))
  est.var = solve(t(smat)%*%(smat))
  err=as.matrix(sqrt(diag(est.var)))
  coef=as.matrix(c(beta,gamma))
  pvalue=sapply(1:length(coef),function(i){2*(1-pnorm(abs(coef[i]/err[i])))})
  
  clinpara=rbind(t(coef),t(err),t(pvalue))
  return(clinpara)
}

LuMixFullPara<-function(obs.time,delta,z){
    n=dim(z)[1]  
    z=as.matrix(z)
    lm.fit = lm(log(obs.time[delta==1]) ~ z[delta==1,])
    beta.ini = lm.fit$coef[-1]
    g.fit = glm(delta ~ z, family = "binomial")
    gamma.ini = g.fit$coef
    e = log(obs.time)-z%*%beta.ini 
    
    l = survfit(Surv(e,delta)~1)
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    s.ini = l$surv
    s.ini = s_ext(s.ini,e,delta)
    
    sigma = sqrt(var(log(obs.time[delta==1])-z[delta==1,]%*%beta.ini))
    h.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
    
    dx = 1.0
    iter = 1
    max.iter=1000
    while(iter <= max.iter && dx >= 0.01){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
      
      ####M-strep
      gamma = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
      beta = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
      
      e = log(obs.time)-z%*%beta
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s = s_ext(s.fit,e,delta)
      
      dx = max(c(abs(beta-beta.ini),abs(gamma-gamma.ini)))
      gamma.ini = gamma
      beta.ini = beta
      s.ini = s
      iter = iter + 1  
      #cat(iter,dx,beta,gamma,"\n")
    }
    
    ########variance estimation based on profile likelihood 
    ###fix beta
    h.n = 2/n
    
    pl.s1=foreach(i=1:length(beta), .combine='cbind', .packages=c("MASS","survival")) %dopar% {
      source("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SIMfunctions.R")
      beta.ini[i] = beta[i] + h.n
      beta.ini[-i] = beta[-i] 
      gamma.ini = gamma
      beta.temp = beta.ini[-i]
      beta.new = numeric(length(beta))
      beta.new[i] = beta.ini[i]
      s.ini = s
      dx = 1.0
      iter = 1
      while(iter <= max.iter && dx >= 0.01){
        ####E-step
        w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
        ####M-strep
        gamma.new = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
        beta.new[-i] = optim(beta.temp,post.lik2n,beta2=beta.ini[i],k=i,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
        e = log(obs.time)-z%*%beta.ini
        index = seq(1,n,1)
        index1 = index[delta[order(e)]==1]
        Lambda = Lam.f(e,delta,w,h.opt)
        s.fit = exp(-Lambda)
        s.new = s_ext(s.fit,e,delta)
        dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
        gamma.ini = gamma.new
        beta.ini = beta.new
        beta.temp = beta.new[-i]
        s.ini = s.new
        iter = iter + 1   
        #cat(iter,dx,beta.ini,gamma.new,"\n")
      }  
      w = as.vector(eta_prob(gamma.new,z,delta,s.new))
      pl=pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)
      return(pl)
    }
    
    pl.s2=foreach(i=1:length(beta), .combine='cbind', .packages=c("MASS","survival")) %dopar% {
      source("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SIMfunctions.R")
      beta.ini[i] = beta[i] - h.n
      beta.ini[-i] = beta[-i] 
      gamma.ini = gamma
      beta.temp = beta.ini[-i]
      beta.new = numeric(length(beta))
      beta.new[i] = beta.ini[i]
      s.ini = s
      dx = 1.0
      iter = 1
      while(iter <= max.iter && dx >= 0.01){
        ####E-step
        w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
        ####M-strep
        gamma.new = optim(gamma.ini,post.lik1,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
        beta.new[-i] = optim(beta.temp,post.lik2n,beta2=beta.ini[i],k=i,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
        e = log(obs.time)-z%*%beta.ini
        index = seq(1,n,1)
        index1 = index[delta[order(e)]==1]
        Lambda = Lam.f(e,delta,w,h.opt)
        s.fit = exp(-Lambda)
        s.new = s_ext(s.fit,e,delta)
        dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
        gamma.ini = gamma.new
        beta.ini = beta.new
        beta.temp = beta.new[-i]
        s.ini = s.new
        iter = iter + 1   
        #cat(iter,dx,beta.ini,gamma.new,"\n")
      }  
      w = as.vector(eta_prob(gamma.new,z,delta,s.new))
      pl=pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)
      return(pl)
    }
    
    pl.s3=foreach(i=1:length(gamma), .combine='cbind', .packages=c("MASS","survival")) %dopar% {
      source("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SIMfunctions.R")
      beta.ini = beta
      gamma.ini[i] = gamma[i] + h.n
      gamma.ini[-i] = gamma[-i]
      gamma.temp = gamma.ini[-i]
      gamma.new = numeric(length(gamma))
      gamma.new[i] = gamma.ini[i]
      s.ini = s
      dx = 1.0
      iter = 1
      while(iter <= max.iter && dx >= 0.012){
        ####E-step
        w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
        
        ####M-strep
        gamma.new[-i] = optim(gamma.temp,post.lik1n,gamma2=gamma.ini[i],k=i,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
        beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
        
        e = log(obs.time)-z%*%beta.ini
        index = seq(1,n,1)
        index1 = index[delta[order(e)]==1]
        Lambda = Lam.f(e,delta,w,h.opt)
        s.fit = exp(-Lambda)
        s.new = s_ext(s.fit,e,delta) 
        
        dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
        gamma.ini = gamma.new
        beta.ini = beta.new
        gamma.temp = gamma.new[-i]
        s.ini = s.new
        iter = iter + 1  
        #cat(iter,dx,beta.ini,gamma.new,"\n")
      }
      w = as.vector(eta_prob(gamma.new,z,delta,s.new))
      pl=pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)
      return(pl)
    }
    
    pl.s4=foreach(i=1:length(gamma), .combine='cbind', .packages=c("MASS","survival")) %dopar% {
      source("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SIMfunctions.R")
      beta.ini = beta
      gamma.ini[i] = gamma[i] - h.n
      gamma.ini[-i] = gamma[-i]
      gamma.temp = gamma.ini[-i]
      gamma.new = numeric(length(gamma))
      gamma.new[i] = gamma.ini[i]
      s.ini = s
      dx = 1.0
      iter = 1
      while(iter <= max.iter && dx >= 0.012){
        ####E-step
        w = as.vector(eta_prob(gamma.ini,z,delta,s.ini))
        
        ####M-strep
        gamma.new[-i] = optim(gamma.temp,post.lik1n,gamma2=gamma.ini[i],k=i,w=w,z=z,control=list(reltol=0.00001,maxit=1000))$par
        beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000))$par
        
        e = log(obs.time)-z%*%beta.ini
        index = seq(1,n,1)
        index1 = index[delta[order(e)]==1]
        Lambda = Lam.f(e,delta,w,h.opt)
        s.fit = exp(-Lambda)
        s.new = s_ext(s.fit,e,delta)
        
        dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
        gamma.ini = gamma.new
        beta.ini = beta.new
        gamma.temp = gamma.new[-i]
        s.ini = s.new
        iter = iter + 1  
        #cat(iter,dx,beta.ini,gamma.new,"\n")
      }
      w = as.vector(eta_prob(gamma.new,z,delta,s.new))
      pl=pl.lik(beta.new,gamma.new,obs.time,z,delta,w,h.opt)
      return(pl)
    }
    
    smat = cbind(sapply(1:length(beta),function(i){((pl.s1[,i]-pl.s2[,i])/(2*h.n))}),sapply(1:length(gamma),function(i){((pl.s3[,i]-pl.s4[,i])/(2*h.n))}))
    est.var = solve(t(smat)%*%(smat))
    err=as.matrix(sqrt(diag(est.var)))
    coef=as.matrix(c(beta,gamma))
    pvalue=sapply(1:length(coef),function(i){2*(1-pnorm(abs(coef[i]/err[i])))})
    return(rbind(t(coef),t(err),t(pvalue)))
}

LuMixFullNew<-function(obs.time,delta,betaz,gammaz,sd){
  n=length(delta)
  z1=as.matrix(betaz)
  z2=as.matrix(gammaz)
  if(is.vector(betaz)){lm.fit = lm(log(obs.time[delta==1]) ~ z1[delta==1])}
  if(!is.vector(betaz)){lm.fit = lm(log(obs.time[delta==1]) ~ z1[delta==1,])}
  beta.ini = lm.fit$coef[-1]
  g.fit = glm(delta ~ z2, family = "binomial")
  gamma.ini = g.fit$coef
    if(length(beta.ini)==1){
      e = log(obs.time)-beta.ini*z1
    }
    if(length(beta.ini)>1){
      e = log(obs.time)-z1%*%beta.ini 
    }
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s.ini = l$surv
  s.ini = s_ext(s.ini,e,delta)
  
  if(length(beta.ini)==1){
    sigma = sqrt(var(log(obs.time[delta==1])-beta.ini*z1[delta==1]))
  }
  if(length(beta.ini)>1){
    sigma = sqrt(var(log(obs.time[delta==1])-z1[delta==1,]%*%beta.ini))
  }
  
  h.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(gamma.ini,z2,delta,s.ini))
    
    ####M-strep
    gamma = optim(gamma.ini,post.lik1,w=w,z=z2,control=list(reltol=0.00001,maxit=1000))$par
    beta = optim(beta.ini,post.lik2,obs.time=obs.time,z=z1,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    if(length(beta)==1){
      e = log(obs.time)-beta*z1
    }
    if(length(beta.ini)>1){
      e = log(obs.time)-z1%*%beta 
    }
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h.opt)
    s.fit = exp(-Lambda)
    s = s_ext(s.fit,e,delta)
    
    dx = max(c(abs(beta-beta.ini),abs(gamma-gamma.ini)))
    gamma.ini = gamma
    beta.ini = beta
    s.ini = s
    iter = iter + 1  
    #cat(iter,dx,beta,gamma,"\n")
  }
  w = as.vector(eta_prob(gamma,z2,delta,s))
  pl.lik=sum(pl.lik(beta,gamma,obs.time,z1,z2,delta,w,h.opt))
  if(sd==0){ParaList<-list("beta"=beta,"gamma"=gamma,"w"=w,"pl.lik"=pl.lik)}
  if(sd==1){para=LuMixFullSD(beta,gamma,s,h.opt,obs.time,delta,z1,z2)
            ParaList<-list("beta"=beta,"gamma"=gamma,"w"=w,"pl.lik"=pl.lik,"para"=para,"s"=s)}
  return(ParaList)}




LuMixFullLRT<-function(obs.time,delta,z){
  cat("Start Alt...",sep="\n")
  Alt=LuMixFullNew(obs.time,delta,z,z,1)
  cat("Start Null...",sep="\n")
  Null=LuMixFullNew(obs.time,delta,z[,1],z[,1],1)
  cat("Start SusNull...",sep="\n")
  SusNull=LuMixFullNew(obs.time,delta,z,z[,1],1)
  cat("Start SurvNull...",sep="\n")
  SurvNull=LuMixFullNew(obs.time,delta,z[,1],z,1)
  
  #Section for Alt SD calculation 
  #para=LuMixFullSD(Alt$beta,Alt$gamma,Alt$s,Alt$h.opt,obs.time,delta,cz,gz)
  
  #Alt.Fix=LuMixFixAlt(obs.time,delta,Null$beta,Null$gamma,cz,gz)
  #SusNull.Fix=LuMixFixSusNull(obs.time,delta,Null$beta,Null$gamma,cz,gz)
  #SurvNull.Fix=LuMixFixSurvNull(obs.time,delta,Null$beta,Null$gamma,cz,gz)
  
  #w.alt.fix = as.vector(eta_prob(Alt.Fix$gamma,z,delta,Alt.Fix$s))
  #w.null.fix = as.vector(eta_prob(Null$gamma,cz,delta,Null$s))
  #w.susn.fix = as.vector(eta_prob(SusNull.Fix$gamma,cz,delta,SusNull.Fix$s))
  #w.survn.fix = as.vector(eta_prob(SurvNull.Fix$gamma,z,delta,SurvNull.Fix$s))
  
  #alt.l.fix=sum(pl.lik(Alt.Fix$beta,Alt.Fix$gamma,obs.time,z,delta,w.alt.fix,Alt.Fix$h.opt))
  #susnull.l.fix=sum(pl.lik(SusNull.Fix$beta,SusNull.Fix$gamma,obs.time,z,delta,w.susn.fix,SusNull.Fix$h.opt))
  #survnull.l.fix=sum(pl.lik(SurvNull.Fix$beta,SurvNull.Fix$gamma,obs.time,z,delta,w.survn.fix,SurvNull.Fix$h.opt))
  
  D.sus=2*(Alt$pl.lik-SusNull$pl.lik)
  D.surv=2*(Alt$pl.lik-SurvNull$pl.lik)
  D.both=2*(Alt$pl.lik-Null$pl.lik)
  
  #D.sus.fix=2*(alt.l.fix-susnull.l.fix)
  #D.surv.fix=2*(alt.l.fix-survnull.l.fix)
  #D.both.fix=2*(alt.l.fix-null.l)
  
  p.sus=pchisq(D.sus,df=1,lower.tail = FALSE)
  p.surv=pchisq(D.surv,df=1,lower.tail = FALSE)
  p.both=pchisq(D.both,df=2,lower.tail = FALSE)
  
  #p.sus.fix=pchisq(D.sus.fix,df=1,lower.tail = FALSE)
  #p.surv.fix=pchisq(D.surv.fix,df=1,lower.tail = FALSE)
  #p.both.fix=pchisq(D.both.fix,df=2,lower.tail = FALSE)
  
  likelihoods=cbind(c(Alt$pl.lik,SurvNull$pl.lik,p.surv),c(Alt$pl.lik,SusNull$pl.lik,p.sus),c(Alt$pl.lik,Null$pl.lik,p.both))#,c(Alt$pl.lik.fix,SurvNull$pl.lik.fix,p.surv.fix),c(Alt$pl.lik.fix,SusNull$pl.lik.fix,p.sus.fix),c(Alt$pl.lik.fix,Null$pl.lik,p.both.fix)
  betas=cbind(c(SurvNull$beta,NA),c(SusNull$beta),c(Null$beta,NA),c(Alt$beta))#,c(SurvNull.Fix$beta,NA),c(SusNull.Fix$beta),c(Null$beta,NA),c(Alt.Fix$beta)
  gammas=cbind(c(SurvNull$gamma),c(SusNull$gamma,NA),c(Null$gamma,NA),c(Alt$gamma))#,c(SurvNull.Fix$gamma),c(SusNull.Fix$gamma,NA),c(Null$gamma,NA),c(Alt.Fix$gamma)
  ParaList<-list("para"=Alt$para,"likelihoods"=likelihoods,"gammas"=gammas,"betas"=betas,"s"=Null$s,"nullpara"=Null$para)
  return(ParaList)
  } 

LuMixFixAlt<-function(obs.time,delta,cbeta,cgamma,cz,z){
  beta.clin=cbeta*cz
  gamma.clin=cgamma*cz
  n=length(z)  
  lm.fit = lm((log(obs.time[delta==1])-beta.clin[delta==1])~z[delta==1])
  beta1.ini = lm.fit$coef[2]
  g.fit = glm(delta ~ z, family = "binomial")
  gamma1.ini = g.fit$coef[2]
  e = log(obs.time)-beta1.ini*z-beta.clin
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s1.ini = l$surv
  s1.ini = s_ext(s1.ini,e,delta)
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta1.ini*z[delta==1]-beta.clin[delta==1]))#
  h1.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.ini))
    
    ####M-strep
 
    gamma1 = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta1 = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z,delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1*z-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1 = s_ext(s1.fit,e,delta)
    
    dx = max(c(abs(beta1-beta1.ini),abs(gamma1-gamma1.ini)))
    gamma1.ini = gamma1
    beta1.ini = beta1
    s1.ini = s1
    iter = iter + 1  
  }
  ParaList<-list("beta"=c(cbeta,beta1),"gamma"=c(cgamma,gamma1),"s"=s1,"h.opt"=h1.opt)
  return(ParaList)
}

LuMixFixSurvNull<-function(obs.time,delta,cbeta,cgamma,cz,z){
  beta.clin=cbeta*cz
  gamma.clin=cgamma*cz
  n=length(z)  
  #lm.fit = lm((log(obs.time[delta==1])-beta.clin[delta==1])~z[delta==1])
  #beta1.ini = lm.fit$coef[2]
  g.fit = glm(delta ~ z, family = "binomial")
  gamma1.ini = g.fit$coef[2]
  e = log(obs.time)-beta.clin
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s1.ini = l$surv
  s1.ini = s_ext(s1.ini,e,delta)
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta.clin[delta==1]))#
  h1.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(feta_prob(gamma1.ini,gamma.clin,z,delta,s1.ini))
    
    ####M-strep
    
    gamma1 = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    #beta1 = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z,delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1 = s_ext(s1.fit,e,delta)
    
    dx = abs(gamma1-gamma1.ini)
gamma1.ini = gamma1
s1.ini = s1
iter = iter + 1  
  }
  ParaList<-list("beta"=cbeta,"gamma"=c(cgamma,gamma1),"s"=s1,"h.opt"=h1.opt)
  return(ParaList)
}

LuMixFixSusNull<-function(obs.time,delta,cbeta,cgamma,cz,z){
  beta.clin=cbeta*cz
  gamma.clin=cgamma*cz
  n=length(z)  
  lm.fit = lm((log(obs.time[delta==1])-beta.clin[delta==1])~z[delta==1])
  beta1.ini = lm.fit$coef[2]
  #g.fit = glm(delta ~ z, family = "binomial")
  #gamma1.ini = g.fit$coef[2]
  e = log(obs.time)-beta1.ini*z-beta.clin
  
  l = survfit(Surv(e,delta)~1)
  index = seq(1,n,1)
  index1 = index[delta[order(e)]==1]
  s1.ini = l$surv
  s1.ini = s_ext(s1.ini,e,delta)
  
  sigma = sqrt(var(log(obs.time[delta==1])-beta1.ini*z[delta==1]-beta.clin[delta==1]))#
  h1.opt = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
  
  dx = 1.0
  iter = 1
  max.iter=1000
  while(iter <= max.iter && dx >= 0.01){
    ####E-step
    w = as.vector(eta_prob(cgamma,cz,delta,s1.ini))
    
    ####M-strep
    
    #gamma1 = optim(gamma1.ini,fpost.lik1,gamma.clin=gamma.clin,w=w,z=z,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    beta1 = optim(beta1.ini,fpost.lik2,beta.clin=beta.clin,obs.time=obs.time,z=z,delta=delta,w=w,h=h1.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
    
    e = log(obs.time)-beta1*z-beta.clin
    index = seq(1,n,1)
    index1 = index[delta[order(e)]==1]
    Lambda = Lam.f(e,delta,w,h1.opt)
    s1.fit = exp(-Lambda)
    s1 = s_ext(s1.fit,e,delta)
    
    dx = abs(beta1-beta1.ini)
    #gamma1.ini = gamma1
    beta1.ini = beta1
    s1.ini = s1
    iter = iter + 1  
  }
  ParaList<-list("beta"=c(cbeta,beta1),"gamma"=cgamma,"s"=s1,"h.opt"=h1.opt)
  return(ParaList)
}


LuMixFullSD=function(beta,gamma,s,h.opt,obs.time,delta,betaz,gammaz){
  z1=betaz
  z2=gammaz
  max.iter=1000
  n=length(obs.time)
  h.n = 2/n
  cat("Start pl.s1...",sep="\n")
  pl.s1=sapply(1:length(beta),function(i){
    beta.ini=rep(NA,length(beta))
    gamma.ini=rep(NA,length(gamma))
    beta.ini[i] = beta[i] + h.n
    beta.ini[-i] = beta[-i] 
    gamma.ini = gamma
    beta.temp = beta.ini[-i]
    beta.new = numeric(length(beta))
    beta.new[i] = beta.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.01){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z2,delta,s.ini))
      ####M-strep
      gamma.new = optim(gamma.ini,post.lik1,w=w,z=z2,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
      if(length(beta.new>1)){
      beta.new[-i] = optim(beta.temp,post.lik2n,beta2=beta.ini[i],k=i,obs.time=obs.time,z=z1,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par}
    
      if(length(beta.new)==1){
        e = log(obs.time)-beta.new*z1}
      if(length(beta.new)>1){
        e = log(obs.time)-z1%*%beta.new}
      
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta)
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      if(length(beta.new)>1){beta.ini = beta.new
      beta.temp = beta.new[-i]}
      if(length(beta.new)==1){beta.ini = beta.new}
      s.ini = s.new
      iter = iter + 1   
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }  
    w = as.vector(eta_prob(gamma.new,z2,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z1,z2,delta,w,h.opt)
  })
  cat("Start pl.s2...",sep="\n")
  pl.s2=sapply(1:length(beta),function(i){
    beta.ini=rep(NA,length(beta))
    gamma.ini=rep(NA,length(gamma))
    beta.ini[i] = beta[i] - h.n
    beta.ini[-i] = beta[-i] 
    gamma.ini = gamma
    beta.temp = beta.ini[-i]
    beta.new = numeric(length(beta))
    beta.new[i] = beta.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.01){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z2,delta,s.ini))
      ####M-strep
      gamma.new = optim(gamma.ini,post.lik1,w=w,z=z2,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
      if(length(beta.new>1)){
        beta.new[-i] = optim(beta.temp,post.lik2n,beta2=beta.ini[i],k=i,obs.time=obs.time,z=z1,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par}
      
      if(length(beta.new)==1){
        e = log(obs.time)-beta.new*z1
      }
      if(length(beta.new)>1){
        e = log(obs.time)-z1%*%beta.new 
      }
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta)
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      if(length(beta.new)>1){beta.ini = beta.new
      beta.temp = beta.new[-i]}
      if(length(beta.new)==1){beta.ini = beta.new}
      s.ini = s.new
      iter = iter + 1   
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }  
    w = as.vector(eta_prob(gamma.new,z2,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z1,z2,delta,w,h.opt)
  })
  cat("Start pl.s3...",sep="\n")
  pl.s3=sapply(1:length(gamma),function(i){
    beta.ini=rep(NA,length(beta))
    gamma.ini=rep(NA,length(gamma))
    beta.ini = beta
    gamma.ini[i] = gamma[i] + h.n
    gamma.ini[-i] = gamma[-i]
    gamma.temp = gamma.ini[-i]
    gamma.new = numeric(length(gamma))
    gamma.new[i] = gamma.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.012){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z2,delta,s.ini))
      
      ####M-strep
      gamma.new[-i] = optim(gamma.temp,post.lik1n,gamma2=gamma.ini[i],k=i,w=w,z=z2,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
      beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z1,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
      
      if(length(beta.new)==1){
        e = log(obs.time)-beta.new*z1
      }
      if(length(beta.new)>1){
        e = log(obs.time)-z1%*%beta.new 
      }
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta) 
      
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      beta.ini = beta.new
      gamma.temp = gamma.new[-i]
      s.ini = s.new
      iter = iter + 1  
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }
    w = as.vector(eta_prob(gamma.new,z2,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z1,z2,delta,w,h.opt) 
  })
  
  cat("Start pl.s4...",sep="\n")
  pl.s4=sapply(1:length(gamma),function(i){
    beta.ini=rep(NA,length(beta))
    gamma.ini=rep(NA,length(gamma))
    beta.ini = beta
    gamma.ini[i] = gamma[i] - h.n
    gamma.ini[-i] = gamma[-i]
    gamma.temp = gamma.ini[-i]
    gamma.new = numeric(length(gamma))
    gamma.new[i] = gamma.ini[i]
    s.ini = s
    dx = 1.0
    iter = 1
    while(iter <= max.iter && dx >= 0.012){
      ####E-step
      w = as.vector(eta_prob(gamma.ini,z2,delta,s.ini))
      
      ####M-strep
      gamma.new[-i] = optim(gamma.temp,post.lik1n,gamma2=gamma.ini[i],k=i,w=w,z=z2,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
      beta.new = optim(beta.ini,post.lik2,obs.time=obs.time,z=z1,delta=delta,w=w,h=h.opt,control=list(reltol=0.00001,maxit=1000),method="BFGS")$par
      
      if(length(beta.new)==1){
        e = log(obs.time)-beta.new*z1
      }
      if(length(beta.new)>1){
        e = log(obs.time)-z1%*%beta.new 
      }
      index = seq(1,n,1)
      index1 = index[delta[order(e)]==1]
      Lambda = Lam.f(e,delta,w,h.opt)
      s.fit = exp(-Lambda)
      s.new = s_ext(s.fit,e,delta)
      
      dx = max(abs(gamma.new-gamma.ini),abs(beta.new-beta.ini))
      gamma.ini = gamma.new
      beta.ini = beta.new
      gamma.temp = gamma.new[-i]
      s.ini = s.new
      iter = iter + 1  
      #cat(iter,dx,beta.ini,gamma.new,"\n")
    }
    w = as.vector(eta_prob(gamma.new,z2,delta,s.new))
    pl.lik(beta.new,gamma.new,obs.time,z1,z2,delta,w,h.opt)
  })
  
  smat = cbind(sapply(1:length(beta),function(i){((pl.s1[,i]-pl.s2[,i])/(2*h.n))}),sapply(1:length(gamma),function(i){((pl.s3[,i]-pl.s4[,i])/(2*h.n))}))
  est.var = solve(t(smat)%*%(smat))
  err=as.matrix(sqrt(diag(est.var)))
  coef=as.matrix(c(beta,gamma))
  pvalue=sapply(1:length(coef),function(i){2*(1-pnorm(abs(coef[i]/err[i])))})
  return(rbind(t(coef),t(err),t(pvalue)))
}


LuMixFull<-function(obs.time,delta,bz,bzd,gz,gzd,sample,out=1,snp=0){
  all_content = readLines(sample)
  skip_second = all_content[c(-2)]
  a1 = read.csv(textConnection(skip_second), header = TRUE, stringsAsFactors = FALSE, sep=" ")
  
  #Sort out NA and 0 times
  a=a1[complete.cases(a1[,c(obs.time,delta,bz,gz)])&a1[,obs.time]!=0&a1[,16]!=1,]
  
  #cat("Start Alt...",sep="\n")
  Alt=LuMixFullNew(a[,obs.time],a[,delta],a[,bz],a[,gz],1)
  
  paras=Alt$para#,c(Alt$pl.lik.fix,SurvNull$pl.lik.fix,p.surv.fix),c(Alt$pl.lik.fix,SusNull$pl.lik.fix,p.sus.fix),c(Alt$pl.lik.fix,Null$pl.lik,p.both.fix)
  colnames(paras)<-c(colnames(a)[bz],"Sus_Intercept",colnames(a)[gz])
  if(snp!=0){
  write.table(paras,file =paste("EpiPGX_",colnames(a)[snp],".txt",sep=""),sep=" ",quote=F,row.names = F,col.names = T)}
  if(snp==0){
    write.table(paras,file ="EpiPGX_CLINS.txt",sep=" ",quote=F,row.names = F,col.names = T)}
  
  betas=Alt$beta#,c(SurvNull.Fix$beta,NA),c(SusNull.Fix$beta),c(Null$beta,NA),c(Alt.Fix$beta)
  gammas=Alt$gamma#,c(SurvNull.Fix$gamma),c(SusNull.Fix$gamma,NA),c(Null$gamma,NA),c(Alt.Fix$gamma)
  if(out==1){
    maxbzd=max(bzd)
    gz_int<-c(0,gzd)
    maxgzd=max(gzd)
    
    sigbeta=as.matrix(Alt$para[1,1:(dim(a[,bz])[2])]*as.numeric(Alt$para[3,1:(dim(a[,bz])[2])]<0.05))
    for(i in seq(1,maxbzd)){sigbeta[bzd==i]<-as.matrix(Alt$para[1,1:(dim(a[,bz])[2])][bzd==i]*as.numeric(min(Alt$para[3,1:(dim(a[,bz])[2])][bzd==i])<0.05))}
    siggamma=as.matrix(Alt$para[1,(dim(a[,bz])[2]+1):((dim(a[,bz])[2]+dim(a[,gz])[2])+1)]*as.numeric(Alt$para[3,(dim(a[,bz])[2]+1):((dim(a[,bz])[2]+dim(a[,gz])[2])+1)]<0.05))
    for(i in seq(1,maxgzd)){siggamma[gz_int==i]<-as.matrix(Alt$para[1,(dim(a[,gz])[2]+1):((dim(a[,gz])[2]+dim(a[,gz])[2])+1)][gz_int==i]*as.numeric(min(Alt$para[3,(dim(a[,gz])[2]+1):((dim(a[,gz])[2]+dim(a[,gz])[2])+1)][gz_int==i])<0.05))}
    #gamma[pvalue[2:3]>0.05]<-0
    neta = pb.f(dim(a[,bz])[1],siggamma,as.matrix(a[,gz]))
    w = as.vector(eta_prob(siggamma,as.matrix(a[,gz]),delta,Alt$s))
    sus.res=rep(NA,dim(a1[,bz])[1])
    sus.res[complete.cases(a1[,c(obs.time,delta,bz,gz)])&a1[,obs.time]!=0]<- w-neta
    surv.res=rep(NA,dim(a1[,bz])[1])
    surv.res[complete.cases(a1[,c(obs.time,delta,bz,gz)])&a1[,obs.time]!=0]<- Alt$s
    
    #beta[pvalue[1]>0.05]<-0
    #rc=-log(s+0.0001)
    #mart2=delta-w*rc
    #surv.dev=sign(mart2)*sqrt(-2*(mart2+delta*log(delta-mart2)))
    #e = log(a[,obs.time])-as.matrix(a[,z])%*%sigbeta
    
    sample.res<-cbind(a1,sus.res,surv.res)
    second = all_content[c(2)]
    header2 = read.csv(textConnection(second), header = F, stringsAsFactors = FALSE, sep=" ",check.names = F)
    nheader2=as.data.frame(c(header2,c("P","P")))
    
    names(nheader2) <- names(sample.res) 
    rbind(nheader2,sample.res)
    
    write.table(nheader2,file = paste(sample,"wres",sep=""),sep=" ",quote=F,row.names = F,col.names = T)
    write.table(sample.res,file = paste(sample,"wres",sep=""),sep=" ",quote=F,row.names = F,col.names = F,append = T)
  }
  ParaList<-list("para"=paras,"gammas"=gammas,"betas"=betas)
  return(ParaList)
} 