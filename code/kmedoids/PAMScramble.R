rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'micekmeans/subsample/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

#########################
### load imputed data ###
#########################

library(parallel)

microSample.imp <- complete(microSample.imp,1)
colnames(microSample.imp) <- colnames(microSample)

list[pathItems.type,pathRegions.name] <- get.pathscore.names(vers = 'short')
list[pathItems.index,pathItems.labels] <- get.feature.labels(pathItems.type,colnames(microSample.imp))
# order microSample columns to appear in Fig. 4d displaying cluster centroids
microSample <- microSample.imp[,pathItems.index]

# 1. generate 10 pathology matrices with independently scrambled columns. 
# 2. Take 10 80% subsamples of each of these matrices, and cluster at the optimal threshold for a range of k values.
# 3.	Plot null distribution of silhouette values as a function of k for optimal threshold.

nobs <- nrow(microSample)
n.scrambles <- 100
microSample.Scramble <- W.scramble <- list()
for(S in 1:n.scrambles){
  microSample.Scramble[[S]] <- sapply(microSample, function(X) X[sample(x = 1:length(X),replace = F)])
  rownames(microSample.Scramble[[S]]) <- 1:nrow(microSample.Scramble[[S]]) # replace rownames
}

# use parallel apply functions to spread polychoric calculation across cores

numCores <- 4 # for some reason detectCores() says I have 8 cores when I know I have 4
parallel.polychoric <- function(X){ return( psych::polychoric(t(X)) ) } # function for polychoric correlation
cl <- makeCluster(numCores)
for(S in 1:n.scrambles){
  print(paste('Scramble',S))
  if(params$dist.met == 'polychoric'){W.scramble[[S]] <- parallel.polychoric(microSample.Scramble[[S]])}
  else if(params$dist.met == 'spearman'){W.scramble[[S]] <- cor(microSample.Scramble[[S]],method='spearman')}
  
}

#W.scramble <- parLapply(cl,microSample.Scramble,parallel.polychoric) # apply microSample.Scramble elements through parallel.polychoric function across cores
save(W.scramble,microSample.Scramble,file=paste0(savedir,'MicroSampleScramble_',params$dist.met,'.RData'))

