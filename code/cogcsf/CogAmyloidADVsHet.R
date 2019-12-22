rm(list = setdiff(ls(), c("params","dz.exc","n.dx")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)

###########
### Cog ###
###########

cog <- read.csv(paste(params$opdir,'processed/MoCA_processed.csv',sep=''))
cog <- cog[cog$INDDID %in% INDDIDs,]

# look by MOCA subscores
COI <- grep('Total',colnames(cog))
cog.moca <- lapply(unique(cog$INDDID), function(ID) 
  colMeans(cog[cog$INDDID == ID,COI],na.rm= T)) # average over multiple tests
cog.moca <- do.call(what='rbind',args=cog.moca)
cog.moca <- as.data.frame(cog.moca[,colSums(is.na(cog.moca)) == 0]) # remove the qualitative features
exclude <- c('MoCATotal')
cog.moca <- cog.moca[,!(names(cog.moca) %in% exclude)]

partitionSample <- as.character(partition[INDDIDs %in% unique(cog$INDDID)])
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))

clusterColors <- getClusterColors(k)

# Select only pts of certain diseases to see how cluster parse heterogeneity within a disease
Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(cog$INDDID),],dz.exc,n.dx)

# only 2 patients left... not worth doing
lapply(as.data.frame(cog.moca), function(X) 
	wilcox.test(X[Rem.Mask & partitionSample=='Cluster 2'],
	X[!Rem.Mask & partitionSample=='Cluster 2'],conf.int=TRUE))