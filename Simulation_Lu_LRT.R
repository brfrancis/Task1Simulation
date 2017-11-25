##efficient Kernel estimation for AFT cure model
##define functions
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE) 
n=as.numeric(args[1])
maf=as.numeric(args[2])
survb=as.numeric(args[3])
susg=as.numeric(args[4])
#i=as.numeric(args[5])
#system("scp -r -P 8022 /home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/Simulation_Lu_LRT.R gaia:/work/projects/epipgx/users/LIV/Task1Surv/Simulation/")
#n=500;maf=0.4;survb=1;susg=-1;i=9
rm(list=setdiff(ls(), c("n","maf","survb","susg")))#,"i"
###############main function for simulation
#library(statmod)
#source("/home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/Simulation_Lu_Sing.R")
#source("/home/ben/Dropbox/Epilepsy/Task1Analysis/Task1Simulation/SIMfunctions.R")
source("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/SIMfunctions.R")
 nodes = scan(Sys.getenv("OAR_NODE_FILE"), what=character(0))
 oar_job_id = as.numeric(Sys.getenv("OAR_JOB_ID"))
 connector = paste0("OAR_JOB_ID=", oar_job_id)
 connector = paste0(connector, " ~jemeras/bin/roarsh")
 comm_type = "PSOCK"
 cl<-makeCluster(nodes, type = comm_type, rshcmd = connector)
 registerDoParallel(cl)
# 
nsim = 1
max.iter = 1000
# 
# clusterEvalQ(cl,library("MASS"))
# clusterEvalQ(cl,library("survival"))
# clusterExport(cl, varlist=c("n","max.iter","maf","survb","susg"))

Sim=NULL
dataSing=NULL
dataFull=NULL
dataFix=NULL
Sim=foreach(i=1:nsim, .combine='rbind', .packages=c("MASS","survival"), .errorhandling="remove") %dopar% {

data=SimData(n,alpha0 = -0.5,beta0 = c(1, 0, survb),gamma0 = c(2,1,0,susg),tau = 20,maf=maf,it=i)
MixLRT=LuMixFullLRT(obs.time=data[,5],delta=data[,4],z=data[,c(1,3)])

ptm<-proc.time()
dataSing=LuMixSing(obs.time = data[,5],delta = data[,4],z = data[,c(1,2,3,6:11)],clin=MixLRT$nullpara,s=MixLRT$s,it=i)
proc.time()-ptm

# MixLRT=LuMixFullLRT(obs.time=data[,5],delta=data[,4],z=data[,c(1,2)])
# dataSing=LuMixSing(obs.time = data[,5],delta = data[,4],z = data[,c(1,3,2,9:11,6:8)],clin=MixLRT$nullpara,s=MixLRT$s,it=i)
#Sim=rbind(dataSing$like,MixLRT$likelihoods)

# MixLRT$para
# cat(" ",sep="\n")
# MixLRT$nullpara
# cat(" ",sep="\n")
# MixLRT$likelihoods
# cat(" ",sep="\n")
# MixLRT$betas
# cat(" ",sep="\n")
# MixLRT$gammas
# cat(" ",sep="\n")
# dataSing$like
# cat(" ",sep="\n")
cbind(MixLRT$likelihoods,dataSing$like)}
Sim
write.table(Sim,file=paste("/work/projects/epipgx/users/LIV/Task1Surv/Simulation/temp/Results",survb,"_",susg,"_",n,"_",maf,".txt",sep=""),quote = F,row.names = F,col.names = F,append = T)
