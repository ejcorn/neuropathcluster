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

dz.res <- train.test.new.GLM(x=df,y=dx,nreps=100)
cluster.res <- train.test.new.GLM(x=df,y=clusters,nreps=100)

########################
### Save performance ###
########################

save(dz.res,df, file = paste(savedir,'predictdz_GLMperf',extralab,'.RData',sep=''))
save(cluster.res,df, file = paste(savedir,'predictcluster_GLMperf',extralab,'.RData',sep=''))

