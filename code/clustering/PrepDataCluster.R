rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
source('code/misc/fxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]
#microSample <- as.data.frame(scale(microSample,center = T))
if(sum((colSums(is.na(microSample)) == nrow(microSample))) > 0){
  break
}

writeMat(paste(params$opdir,'processed/pathDataForClustering.mat',sep=''),X=as.matrix(microSample-.1)) # subtract 0.1 to make R save as double
