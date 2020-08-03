rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'micekmeans/bootstrap/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

####################################################################################
### Process resampled clustering data from PAMBootstrapReliabilityScrambleNorm.R ###
####################################################################################

sampfrac <- 1
#k.rng <- 2:10
#thresh.rng <- seq(0,0.5,0.1)
load(file = paste0(savedir,'ClusterBootstrapping_SF',sampfrac,'.RData'))
load(file = paste0(savedir,'ScrambledClusterBootstrapping_SF',sampfrac,'.RData'))

# loop through all clustering on bootstrapped data, pull out silhouette values for each k and threshold combo
df.sub <- lapply(thresh.rng,function(thresh) do.call('rbind',lapply(k.rng, function(k) do.call('rbind',lapply(1:nreps, function(R)
  data.frame(sil=clust.thresh[[R]][[which(thresh.rng==thresh)]][[which(k.rng==k)]]$silinfo$avg.width,
             k=k,thresh=thresh,rep=R,null.or.data='Data'))))))
names(df.sub) <- as.character(thresh.rng)
df.sub.all <- do.call('rbind',df.sub)
silhouette.by.thresh  <- sapply(1:nreps, function(R) 
  sapply(thresh.rng, function(thresh) mean(df.sub.all[df.sub.all$rep==R & df.sub.all$thresh==thresh,'sil'])))

# pick threshold with highest mean silhouette
thresh.bestsil <- thresh.rng[which.max(rowMeans(silhouette.by.thresh))]

# 3. for each Bootstrap, compute similarity of partitions wrt overlappinig elements
# http://psb.stanford.edu/psb-online/proceedings/psb02/benhur.pdf
# need a similarity metric: matching coefficient, jaccard, rand ... should be similar
# alternatively can predict cluster for remaining samples
# https://people.eecs.berkeley.edu/~jordan/sail/readings/luxburg_ftml.pdf
#
# bootstrapping:
# http://www.indigo.lib.uic.edu:8080/bitstream/handle/10027/8612/bootstrap_clustering_revision.pdf?sequence=1&isAllowed=y
# k-medoids: Kaufman and Rousseeuw, 1990

partition.sim.mat <- list()
for(k in k.rng){
  print(paste('k =',k))
  partition.sim.mat[[k]] <- matrix(NA,nreps,nreps)
  for(R1 in 1:nreps){
    for(R2 in 1:nreps){
      partition1 <- clust.thresh[[R1]][[which(thresh.rng==thresh.bestsil)]][[which(k.rng==k)]]$clustering
      partition2 <- clust.thresh[[R2]][[which(thresh.rng==thresh.bestsil)]][[which(k.rng==k)]]$clustering
      
      # extract partition on joint sample
      partition1.union <- partition1[names(partition1) %in% names(partition2)]
      partition2.union <- partition2[names(partition2) %in% names(partition1)]
      partition2.union <- partition2.union[names(partition1.union)]
      partition.sim.mat[[k]][R1,R2] <- cluster_similarity(partition1.union,partition2.union,similarity = 'jaccard')
    }
  }
}

partition.sim.mat.scramble <- list()
for(k in k.rng){
  print(paste('k =',k))
  partition.sim.mat.scramble[[k]] <- matrix(NA,nreps,nreps)
  for(R1 in 1:nreps){
    for(R2 in 1:nreps){
      partition1 <- clust.thresh.scramble[[R1]][[which(thresh.rng==thresh.bestsil)]][[which(k.rng==k)]]$clustering
      partition2 <- clust.thresh.scramble[[R2]][[which(thresh.rng==thresh.bestsil)]][[which(k.rng==k)]]$clustering
      
      # extract partition on joint sample
      partition1.union <- partition1[names(partition1) %in% names(partition2)]
      partition2.union <- partition2[names(partition2) %in% names(partition1)]
      partition2.union <- partition2.union[names(partition1.union)]
      partition.sim.mat.scramble[[k]][R1,R2] <- cluster_similarity(partition1.union,partition2.union,similarity = 'jaccard')
    }
  }
}

save(partition.sim.mat,partition.sim.mat.scramble,thresh.bestsil,file=paste0(savedir,'PartitionSimilarity_ActualAndScramble_ClusterBootstrapping_SF',sampfrac,'.RData'))

