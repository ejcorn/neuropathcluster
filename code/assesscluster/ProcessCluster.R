rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]
#microSample <- as.data.frame(scale(microSample,center = T))
if(sum((colSums(is.na(microSample)) == nrow(microSample))) > 0){
  break
}


#################
### Load data ###
#################

pathItems.type <- c("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody")
list[pathItems.index,pathItems.labels] <- get.feature.labels(pathItems.type,colnames(microSample))
# order microSample columns to appear in Fig. 4d displaying cluster centroids
microSample <- microSample[,pathItems.index]

#cl.data <- readMat('data/subjectClusterLouvain.mat')
cl.data <- readMat(paste(params$opdir,'optimcluster/subjectClusterLouvainPartitionsByGammaNReps1000.mat',sep=''))
partitions.by.gamma <- cl.data$partitions.by.gamma
partitions.by.gamma <- round(partitions.by.gamma)
gamma.opt = params$gamma.opt
gamma.rng <- round(cl.data$gamma.rng*10) /10 # doubles got saved with 1e-16 error
partition <- partitions.by.gamma[,which(gamma.rng == gamma.opt)]

DisconnectedSubjects <- cl.data$DisconnectedSubjects

############################
### Annex small clusters ###
############################

thrsh <- 0.02 # if cluster is < 1% of sample, annex it to another cluster

cluster.counts.by.gamma <- lapply(1:ncol(partitions.by.gamma), function(i)
  sapply(1:max(partitions.by.gamma[,i]), function(k.i) sum(partitions.by.gamma[,i]==k.i)))
names(cluster.counts.by.gamma) <- as.character(gamma.rng)
centroids <- compute.centroids(microSample[-DisconnectedSubjects,],partition)

partition <- annex.small.clusters(partition,X = microSample[-DisconnectedSubjects,],centroids = centroids,thrsh = thrsh)
DisconnectedSubjects <- c(DisconnectedSubjects,which(is.na(partition)))
partition <- partition[!is.na(partition)] # remove disconnected subjects
#list[partition,DisconnectedSubjects] <- disconnect.small.clusters(partition)

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
partition <- remove.Partition.Gap(partition)
centroids <- compute.centroids(microSample,partition)

########################
### Reorder clusters ###
########################

# order clusters the same way 
# regardless of the arbitrary ordering of louvain
# only if 4 clusters

cluster.init.reorder <- order.cluster.by.feature.old(microSample,centroids,c('CB_Tau','Thio','TDP43','Syn'))
# vs.
#cluster.init.reorder <- order.cluster.by.feature(microSample,centroids)

# only reorder if you can make a unique match for each cluster
if(length(unique(cluster.init.reorder)) == length(cluster.init.reorder) & length(cluster.init.reorder) == length(unique(partition))){
	partition <- reorder.partition(partition,cluster.init.reorder)
}
# recompute centroids with reordered partition
centroids <- compute.centroids(microSample,partition)

k <- max(partition)
# save new partition
save(partition,centroids,k,DisconnectedSubjects,pathItems.labels,
	file = paste(savedir,'subjLouvainPartitionReordered.RData',sep=''))
