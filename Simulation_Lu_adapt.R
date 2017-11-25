##Simulation for Mixture Model

##pull in and define variables

rm(list=ls())
args <- commandArgs(trailingOnly = TRUE) 
n=as.numeric(args[1])
maf=as.numeric(args[2])
#n=100;maf=0.01
rm(list=setdiff(ls(), c("n","maf")))

##pull in functions and set up parallel network

#library(statmod)
#source("/home/ben/Dropbox/Epilepsy/Simulation/SIMfunctions.R")
source("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SIMfunctions.R")
nodes = scan(Sys.getenv("OAR_NODE_FILE"), what=character(0))
oar_job_id = as.numeric(Sys.getenv("OAR_JOB_ID"))
connector = paste0("OAR_JOB_ID=", oar_job_id)
connector = paste0(connector, " ~jemeras/bin/roarsh")
comm_type = "PSOCK"
cl<-makeCluster(nodes, type = comm_type, rshcmd = connector)
registerDoParallel(cl)

##Simulation parameters

nsim = 500
max.iter = 1000

##Disperse functions to individual cores

clusterEvalQ(cl,library("MASS"))
clusterEvalQ(cl,library("survival"))
clusterExport(cl, varlist=c("n","max.iter","maf"))

##Begin simulations

l=0
Null_Sim=NULL
while(((nsim*9)-l)>0){
  Null_Sima=foreach(i=1:(nsim-(l/9)), .combine='rbind', .packages=c("MASS","survival"), .errorhandling="remove") %dopar% {
    data=SimData(n,alpha0 = -0.5,beta0 = c(0, 0, 0),gamma0 = c(0,0,0,0),tau = 8,maf=maf,it=i)
    dataSing=LuMixSing(obs.time = data[,5],delta = data[,4],z = data[,c(1:3,6:11)],it=i) 
    dataFull=LuMixFull(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    dataFix = LuMixFix(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    rbind(dataSing,dataFull,dataFix)
  }
  Null_Sima=Null_Sima[complete.cases(Null_Sima),]
  if(!is.null(Null_Sima)){Null_Sim=rbind(Null_Sim,Null_Sima)}
  l=dim(Null_Sim)[1]
  if(is.null(Null_Sim)){l=0}
}

write.table(Null_Sim,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/NullSimResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F)

l=0
Lu_Sim=NULL
while(((nsim*9)-l)>0){
  Lu_Sima=foreach(i=1:(nsim-(l/9)), .combine='rbind', .packages=c("MASS","survival"), .errorhandling="remove") %dopar% {
    data=SimData(n,alpha0 = -0.5,beta0 = c(1.0, 0, 0),gamma0 = c(0.5,-0.5,0,0),tau = 8,maf=maf,it=i)
    dataSing=LuMixSing(obs.time = data[,5],delta = data[,4],z = data[,c(1:3,6:11)],it=i) 
    dataFull=LuMixFull(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    dataFix = LuMixFix(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    rbind(dataSing,dataFull,dataFix)
  }
  Lu_Sima=Lu_Sima[complete.cases(Lu_Sima),]
  if(!is.null(Lu_Sima)){Lu_Sim=rbind(Lu_Sim,Lu_Sima)}
  l=dim(Lu_Sim)[1]
  if(is.null(Lu_Sim)){l=0}
}

write.table(Lu_Sim,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/LuSimResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F)

l=0
Lu_Surv_Sim=NULL
while(((nsim*9)-l)>0){
  Lu_Surv_Sima=foreach(i=1:(nsim-(l/9)), .combine='rbind', .packages=c("MASS","survival"), .errorhandling="remove") %dopar% {
    data=SimData(n,alpha0 = -0.5,beta0 = c(1.0, 0.5, 0),gamma0 = c(0.5,-0.5,0,0),tau = 8,maf=maf,it=i)
    dataSing=LuMixSing(obs.time = data[,5],delta = data[,4],z = data[,c(1:3,6:11)],it=i) 
    dataFull=LuMixFull(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    dataFix = LuMixFix(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    rbind(dataSing,dataFull,dataFix)
  }
  Lu_Surv_Sima=Lu_Surv_Sima[complete.cases(Lu_Surv_Sima),]
  if(!is.null(Lu_Surv_Sima)){Lu_Surv_Sim=rbind(Lu_Surv_Sim,Lu_Surv_Sima)}
  l=dim(Lu_Surv_Sim)[1]
  if(is.null(Lu_Surv_Sim)){l=0}
}


write.table(Lu_Surv_Sim,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Lu_SurvResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F)

l=0
Lu_Sus_Sim=NULL
while(((nsim*9)-l)>0){
  Lu_Sus_Sima=foreach(i=1:(nsim-(l/9)), .combine='rbind', .packages=c("MASS","survival"), .errorhandling="remove") %dopar% {
    data=SimData(n,alpha0 = -0.5,beta0 = c(1.0, 0, 0),gamma0 = c(0.5,-0.5,0.5,0),tau = 8,maf=maf,it=i)
    dataSing=LuMixSing(obs.time = data[,5],delta = data[,4],z = data[,c(1:3,6:11)],it=i) 
    dataFull=LuMixFull(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    dataFix = LuMixFix(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    rbind(dataSing,dataFull,dataFix)
  }
  Lu_Sus_Sima=Lu_Sus_Sima[complete.cases(Lu_Sus_Sima),]
  if(!is.null(Lu_Sus_Sima)){Lu_Sus_Sim=rbind(Lu_Sus_Sim,Lu_Sus_Sima)}
  l=dim(Lu_Sus_Sim)[1]
  if(is.null(Lu_Sus_Sim)){l=0}
}

write.table(Lu_Sus_Sim,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Lu_SusResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F)

l=0
Lu_Both_Sim=NULL
while(((nsim*9)-l)>0){
  Lu_Both_Sima=foreach(i=1:(nsim-(l/9)), .combine='rbind', .packages=c("MASS","survival"), .errorhandling="remove") %dopar% {
    data=SimData(n,alpha0 = -0.5,beta0 = c(1.0, 0.5, 0),gamma0 = c(0.5,-0.5,0.5,0),tau = 8,maf=maf,it=i)
    dataSing=LuMixSing(obs.time = data[,5],delta = data[,4],z = data[,c(1:3,6:11)],it=i) 
    dataFull=LuMixFull(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    dataFix = LuMixFix(obs.time = data[,5],delta = data[,4],z = data[,c(1:3)])
    rbind(dataSing,dataFull,dataFix)
  }
  Lu_Both_Sima=Lu_Both_Sima[complete.cases(Lu_Both_Sima),]
  if(!is.null(Lu_Both_Sima)){Lu_Both_Sim=rbind(Lu_Both_Sim,Lu_Both_Sima)}
  l=dim(Lu_Both_Sim)[1]
  if(is.null(Lu_Both_Sim)){l=0}
}

write.table(Lu_Both_Sim,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/Lu_BothResults_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F)

res_coll(n,maf)
