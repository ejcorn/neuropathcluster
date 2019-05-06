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
# Only remove Braak-CERAD AD from Cluster 2
Rem.Mask <- Rem.Mask | partitionSample != paste('Cluster',exc.cl)

list[patientSubsample,df,clusters,partitionSample] <- 
  lapply(list(patientSubsample,df,clusters,partitionSample), function(X) flexdim.rowmask(X,Rem.Mask))

print(paste('Train/Test Exclude:',dz.exc))
cluster.counts <- colSums(clusters) # for dummy data frame, colSums == cluster counts output

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample)
list[patientSubsample,df,clusters,partitionSample] <- 
  lapply(list(patientSubsample,df,clusters,partitionSample), function(X) flexdim.rowmask(X,Rem.Cl))

# display counts after excluding small clusters
cluster.counts <- cluster.count(partitionSample,k)
print(paste('CSF: Exclude',dz.exc))
print(cluster.counts)

############################
### Train and test model ###
############################

cluster.res <- kfold.GLM(x=df,y=clusters,nreps=100,k.folds=5)

########################
### Save performance ###
########################

save(cluster.res, file = paste(savedir,'predictcluster_GLMperf',extralab,'Exclude',dz.exc,'.RData',sep=''))
