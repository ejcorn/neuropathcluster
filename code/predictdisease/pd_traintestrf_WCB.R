rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir.in <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir.in,recursive=T)
savedir.out <- paste(params$resultsdir,'predictdisease/weighted_downsample/',sep='')
dir.create(savedir.out,recursive = T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')
source('code/misc/kfold.R')

#################
### Load data ###
#################

load(file = paste(savedir.in,'dzpredict_data',extralab,'.RData',sep=''))

############################
### Train and test model ###
############################
# implement weighted class balancing 
WCB.classes <- c('AD','LBD','FTLD','Normal')
df$WCB_LABELS <- row.Which.Max.tie(dx[,WCB.classes])
dz.res <- kfold.RF(x=df,y=dx,k.folds=5,nreps=100)
save(dz.res, file = paste(savedir.out,'predictdz_WCB_RFperf',extralab,'.RData',sep=''))

df$WCB_LABELS <- row.Which.Max.tie(clusters)
cluster.res <- kfold.RF(x=df,y=clusters,k.folds=5,nreps=100)
save(cluster.res, file = paste(savedir.out,'predictcluster_WCB_RFperf',extralab,'.RData',sep=''))

