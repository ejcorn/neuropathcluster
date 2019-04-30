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

###########################################
##### Neuropath correlation structure #####
###########################################

pathItems.type <- list("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody")
pathItems.index <- sapply(1:length(pathItems.type), function(i) grep(pathItems.type[[i]], colnames(microSample)))
pathItems.labels <- sapply(1:length(pathItems.type), function(i) c(matrix("",floor(0.5*length(pathItems.index[[i]]))),
                                                                   pathItems.type[[i]], c(matrix("",ceiling(0.5*length(pathItems.index[[i]])-1)))))
# initialize everyone's cluster by the type of pathology they have most of
cluster.init.names <- c('Tau','TDP43','Thio','Syn')
cluster.init.idx <- sapply(cluster.init.names, function(i) which(pathItems.type == i))
cluster.init.amount <- sapply(cluster.init.idx, function(i) rowMeans(microSample[,pathItems.index[[i]]],na.rm=T))
cluster.init.assign <- row.Which.Max(cluster.init.amount)

writeMat(paste(params$opdir,'processed/pathDataForClustering.mat',sep=''),M0=cluster.init.assign,X=as.matrix(microSample-.1)) # subtract 0.1 to make R save as double
