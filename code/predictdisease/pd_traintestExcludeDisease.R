rm(list = setdiff(ls(), c("params","extralab",'dz.exc','n.dx','exc.cl')))
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

#################################
### Perform disease exclusion ###
#################################

# Mask that will remove all pts of certain diseases
Rem.Mask <- exclude.dz(patientSubsample,dz.exc,n.dx)

partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))
# Only remove ABC int-high AD from specified Cluster(s)
Rem.Mask <- Rem.Mask | !partitionSample %in% paste('Cluster',exc.cl)
INDDIDs.exclude <- names(partitionSample)[!Rem.Mask] # get INDDIDs of subjects that are to be excluded

list[patientSubsample,partitionSample] <- 
  lapply(list(patientSubsample,partitionSample), function(X) flexdim.rowmask(X,Rem.Mask))

# index data frames by INDDID
df <- df[!df$INDDID %in% INDDIDs.exclude,]
clusters <- clusters[!rownames(clusters) %in% INDDIDs.exclude,]
dx <- dx[!rownames(dx) %in% INDDIDs.exclude,]

print(paste('Train/Test Exclude:',dz.exc))
cluster.counts <- colSums(clusters) # for dummy data frame, colSums == cluster counts output

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample,n=3)
INDDIDs.exclude <- names(partitionSample)[!Rem.Cl] # get INDDIDs of subjects that are to be excluded

list[patientSubsample,partitionSample] <- 
  lapply(list(patientSubsample,partitionSample), function(X) flexdim.rowmask(X,Rem.Cl))

# index data frames by INDDID
df <- df[!df$INDDID %in% INDDIDs.exclude,]
clusters <- clusters[!rownames(clusters) %in% INDDIDs.exclude,]
dx <- dx[!rownames(dx) %in% INDDIDs.exclude,]

# display counts after excluding small clusters
cluster.counts <- cluster.count(partitionSample,k)
print(paste('CSF: Exclude',dz.exc))
print(cluster.counts)
clusters <- clusters[,cluster.counts != 0] # remove columns if Rem.Cl removed that cluster

############################
### Train and test model ###
############################

#dz.res <- kfold.GLM(x=df,y=dx,k.folds=5,nreps=100)
#cluster.res <- kfold.GLM(x=df,y=clusters,nreps=100,k.folds=5)

########################
### Save performance ###
########################

#save(cluster.res, file = paste(savedir,'predictcluster_GLMperf',extralab,'Exclude',dz.exc,'.RData',sep=''))
#save(dz.res,df, file = paste(savedir,'predictdz_GLMperf',extralab,'Exclude',dz.exc,'.RData',sep=''))

#####################
### Random forest ###
#####################

#dz.res <- kfold.RF(x=df,y=dx,nreps=100,k.folds=5)
cluster.res <- kfold.RF(x=df,y=clusters,nreps=100,k.folds=5)

save(cluster.res, file = paste(savedir,'predictcluster_RFperf',extralab,'Exclude',dz.exc,'.RData',sep=''))
#save(dz.res,df, file = paste(savedir,'predictdz_RFperf',extralab,'Exclude',dz.exc,'.RData',sep=''))
