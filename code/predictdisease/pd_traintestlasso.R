rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')

#################
### Load data ###
#################

load(file = paste(savedir,'dzpredict_data',extralab,'.RData',sep=''))

############################
### Train and test model ###
############################

dz.res <- train.test.Lasso(x=df,y=dx,nreps=500)
cluster.res <- train.test.Lasso(x=df,y=clusters,nreps=500)

########################
### Save performance ###
########################

save(dz.res, file = paste(savedir,'predictdz_Lassoperf',extralab,'.RData',sep=''))
save(cluster.res, file = paste(savedir,'predictcluster_Lassoperf',extralab,'.RData',sep=''))

