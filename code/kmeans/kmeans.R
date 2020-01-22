rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'micekmeans/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(savedir,'microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

##################
### Clustering ###
##################

library(cluster)

#subject.cormat <- polychoric(t(microSample.imp))
#imagesc(subject.cormat$rho)
#clust <- hclust(dist(microSample.imp/max(microSample.imp)))
microSample.imp <- complete(microSample.imp,1)
X <- microSample.imp/max(microSample.imp) # perform rank normalization before use of euclidean distance metric
k.rng <- 2:15
clust <- lapply(k.rng, function(k) kmeans(x = X,centers = k,nstart = 100))
sil <- lapply(clust, function(K) silhouette(x=K$cluster,dist=dist(X)))
sil.mean <- sapply(sil,function(S) mean(S[,3]))
p <- ggplot() + geom_line(aes(x=k.rng,y=sil.mean)) + #theme_classic() +
  scale_x_continuous(breaks = k.rng) + xlab('k') + ylab('Mean Silhouette')
ggsave(filename = paste0(savedir,'KMeansMeanSilhouetteByK.pdf'),plot = p,
       height = 18,width=18,units='cm')
clust <- clust[which(k.rng <= 6)]
p.list <- lapply(clust,function(X) imagesc(t(fliplr(micro.order.by(X$centers)))) + 
                   theme(axis.text.y = element_text(size=4,hjust=1,vjust=0.5)))
p <- plot_grid(plotlist = p.list)
ggsave(filename = paste0(savedir,'KMeansCentroids.pdf'),plot = p,
       height = 18,width=18,units='cm')

pathItems.type <- c("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody")
list[pathItems.index,pathItems.labels] <- get.feature.labels(pathItems.type,colnames(microSample))
# order microSample columns to appear in Fig. 4d displaying cluster centroids
microSample <- microSample[,pathItems.index]

for(k in c(4,6)){
  partition <- clust[[which(k.rng==k)]]$cluster
  DisconnectedSubjects <- c() # no disconnected subjects with k-means
  centroids <- compute.centroids(microSample,partition)
  
  # reorder clusters
  cluster.init.reorder <- order.cluster.by.feature.old(microSample,centroids,c('CBTau','Thio','TDP43','Syn'))
  
  # only reorder if you can make a unique match for each cluster
  if(length(unique(cluster.init.reorder)) == length(cluster.init.reorder) & length(cluster.init.reorder) == length(unique(partition))){
    partition <- reorder.partition(partition,cluster.init.reorder)
  }
  centroids <- compute.centroids(microSample,partition)
  # save such that gamma.opt parameter can index k
  savedir <- paste0(params$opdir,'results_G',k,'/analyzecluster/')
  dir.create(savedir,recursive = T)
  save(partition,centroids,k,DisconnectedSubjects,pathItems.labels,
       file = paste0(savedir,'subjLouvainPartitionReordered.RData'))
}


