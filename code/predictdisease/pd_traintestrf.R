rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')
source('code/misc/kfold.R')

#################
### Load data ###
#################

load(file = paste(savedir,'dzpredict_data',extralab,'.RData',sep=''))

############################
### Train and test model ###
############################

dz.res <- kfold.RF(x=df,y=dx,k.folds=5,nreps=100)
cluster.res <- kfold.RF(x=df,y=clusters,k.folds=5,nreps=100)
 
########################
### Save performance ###
########################

save(dz.res, file = paste(savedir,'predictdz_RFperf',extralab,'.RData',sep=''))
save(cluster.res, file = paste(savedir,'predictcluster_RFperf',extralab,'.RData',sep=''))

