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

##############################
### load polychoric matrix ###
##############################

library(cluster)

microSample.imp <- complete(microSample.imp,1)
colnames(microSample.imp) <- colnames(microSample)

if(params$dist.met == 'spearman'){
  DisconnectedSubjects <- which(rowSD(microSample.imp)==0)
  W <- cor(t(microSample.imp[-DisconnectedSubjects,]),method = 'spearman')
} else if(params$dist.met == 'polychoric'){
  load(file = paste0(params$opdir,'processed/SubjectPolychoricMatrix.RData'))
  W <- subject.cormat$rho
}

#####################################################
### see if you get same centroids in split halves ###
#####################################################

# Take 60% subsamples. Cluster at range of k and thresholds.
# 1.	Plot distribution of silhouette values as a function of k and threshold. 
# 2.	For each subsample, figure out best threshold (by peak silhouette.
# 3.	For each subsample, attempt to match centroids to full sample, and then plot distribution of agreement (for each cluster or aggregated).

list[pathItems.type,pathRegions.name] <- get.pathscore.names(vers = 'short')
list[pathItems.index,pathItems.labels] <- get.feature.labels(pathItems.type,colnames(microSample.imp))
# order microSample columns to appear in Fig. 4d displaying cluster centroids
microSample <- microSample.imp[,pathItems.index]

nreps <- 100

nobs <- length(colnames(W))
sampfrac <- 0.8

clust.thresh <- list() # save clustering of subsamples
k.rng <- 2:10
thresh.rng <- seq(0,0.5,0.1)

# Take 60% subsamples. Cluster at range of k and thresholds.
for(R in 1:nreps){
  print(paste0('Subsample #',R))
  # define split sample
  samp <- sample(x = colnames(W),replace = F,size = floor(sampfrac*nobs))
  # Similarity matrix,subsampled, thresholded, converted to dissimilarity matrix
  D.samp <- lapply(thresh.rng, function(thresh) 1-thresh.mat(W[samp,samp],'>',thresh)) 
  start_time <- Sys.time()
  clust.thresh[[R]] <- lapply(D.samp, function(D) lapply(k.rng, function(k)
    pam(x = D,diss = TRUE,k=k))) # cluster subsampled matrix for a range of threshold and k values
  end_time <- Sys.time()
  print(paste(end_time - start_time,'seconds'))
}

save(clust.thresh, file = paste0(savedir,'ClusterSubsampling_SF',sampfrac,'.RData'))

# 1.	Plot distribution of silhouette values as a function of k and threshold. 

# loop through all clustering on subsampled data, pull out silhouette values for each k and threshold combo
df.sub <- lapply(thresh.rng,function(thresh) do.call('rbind',lapply(k.rng, function(k) do.call('rbind',lapply(1:nreps, function(R)
  data.frame(sil=clust.thresh[[R]][[which(thresh.rng==thresh)]][[which(k.rng==k)]]$silinfo$avg.width,
             k=k,thresh=thresh,rep=R))))))
names(df.sub) <- as.character(thresh.rng)

max.sil <- max(do.call('rbind',df.sub)$sil) # scale y axis of all plots the same
p.list <- lapply(as.character(thresh.rng),function(t) ggplot() + #geom_line(data=df.full[[t]],aes(x=k,y=sil)) + #theme_classic() +
                   geom_boxplot(data=df.sub[[t]],aes(x=as.character(k),y=sil))+
                   scale_x_discrete(limits=as.character(k.rng),breaks = as.character(k.rng)) + 
                   scale_y_continuous(limits=c(-0.1,max.sil)) +
                   xlab('k') + ylab('Mean Silhouette') + ggtitle(paste('Threshold =',t))+ theme_bw()+
                   theme(text=element_text(size=8))+
                   theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=6),text=element_text(color='black')))
p <- plot_grid(plotlist = p.list)
ggsave(filename = paste0(savedir,'PAM_SubsampleMeanSilhouetteByKByThresh_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 12,width=18,units='cm')

# 2.	For each subsample, figure out best threshold (by peak silhouette.

df.sub.all <- do.call('rbind',df.sub)
silhouette.by.thresh  <- sapply(1:nreps, function(R) 
  sapply(thresh.rng, function(thresh) mean(df.sub.all[df.sub.all$rep==R & df.sub.all$thresh==thresh,'sil'])))
thresh.indicator <- sapply(1:nreps, function(R) sapply(thresh.rng, function(thresh) thresh)) # for plot, store threshold values
df.plot <- data.frame(sil=as.vector(silhouette.by.thresh),thresh=as.vector(thresh.indicator))
p <- ggplot(df.plot) + geom_boxplot(aes(x=as.character(thresh),y=sil)) + theme_classic() +
  theme(text=element_text(size=8))+
  ylab('Mean Silhouette') + xlab('Threshold')
ggsave(filename = paste0(savedir,'PAM_SubsampleMeanSilhouetteByThresh_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 4,width=4,units='cm')

# pick threshold with highest mean silhouette
thresh.bestsil <- thresh.rng[which.max(rowMeans(silhouette.by.thresh))]

# 4. for each subsample, compute similarity of partitions wrt overlappinig elements
# http://psb.stanford.edu/psb-online/proceedings/psb02/benhur.pdf
# need a similarity metric: matching coefficient, jaccard, rand ... should be similar
# alternatively can predict cluster for remaining samples
# https://people.eecs.berkeley.edu/~jordan/sail/readings/luxburg_ftml.pdf
#
# bootstrapping:
# http://www.indigo.lib.uic.edu:8080/bitstream/handle/10027/8612/bootstrap_clustering_revision.pdf?sequence=1&isAllowed=y
# k-medoids: Kaufman and Rousseeuw, 1990
library(clusteval)

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
mean.jac <- sapply(partition.sim.mat, function(X) mean(X*as.numeric(!diag(nreps)))) # ignore diagonals
mean.jac <- mean.jac[-1] # remove k =1

k.cdf <- lapply(partition.sim.mat[-1], function(X) ecdf(X)(seq(0,1,length.out = 100)))
k.auc <- sapply(k.cdf, mean)

# compute stability metrics in null data

load(file=paste0(savedir,'MicroSampleScramble_',params$dist.met,'.RData'))
n.scramble <- length(W.scramble)
nreps <- 10 # use 10 subsamples for each scrambled matrix
clust.thresh.scramble <- list() # save clustering of subsamples

# Take X% subsamples. Cluster at range of k and thresholds.
for(S in 1:n.scramble){
  for(R in 1:nreps){
    print(paste0('Subsample #',10*(S-1)+R))
    # define split sample
    samp <- sample(x = colnames(W.scramble[[S]]$rho),replace = F,size = floor(sampfrac*nobs))
    # Similarity matrix,subsampled, thresholded, converted to dissimilarity matrix
    D.samp <- lapply(thresh.rng, function(thresh) 1-thresh.mat(W.scramble[[S]]$rho[samp,samp],'>',thresh)) 
    start_time <- Sys.time()
    clust.thresh.scramble[[10*(S-1)+R]] <- lapply(D.samp, function(D) lapply(k.rng, function(k)
      pam(x = D,diss = TRUE,k=k))) # cluster subsampled matrix for a range of threshold and k values
    end_time <- Sys.time()
    print(paste(end_time - start_time,'seconds'))
  }
}

df.sub.scramble <- lapply(thresh.rng,function(thresh) do.call('rbind',lapply(k.rng, function(k) do.call('rbind',lapply(1:nreps, function(R)
  data.frame(sil=clust.thresh.scramble[[R]][[which(thresh.rng==thresh)]][[which(k.rng==k)]]$silinfo$avg.width,
             k=k,thresh=thresh,rep=R))))))
names(df.sub.scramble) <- as.character(thresh.rng)

max.sil <- max(do.call('rbind',df.sub)$sil) # scale y axis of all plots the same -- based on data
p.list <- lapply(as.character(thresh.rng),function(t) ggplot() + #geom_line(data=df.full[[t]],aes(x=k,y=sil)) + #theme_classic() +
                   geom_boxplot(data=df.sub.scramble[[t]],aes(x=as.character(k),y=sil))+
                   scale_x_discrete(limits=as.character(k.rng),breaks = as.character(k.rng)) + 
                   scale_y_continuous(limits=c(-0.1,max.sil)) +
                   xlab('k') + ylab('Mean Silhouette') + ggtitle(paste('Threshold =',t))+ theme_bw()+
                   theme(text=element_text(size=8))+
                   theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=6),text=element_text(color='black')))
p <- plot_grid(plotlist = p.list)
ggsave(filename = paste0(savedir,'PAM_ScrambleSubsampleMeanSilhouetteByKByThresh_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 12,width=18,units='cm')

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

mean.jac.scramble <- sapply(partition.sim.mat.scramble[-1], function(X) mean(X*as.numeric(!diag(nreps)))) # ignore diagonals, remove k=1

p <- ggplot() + geom_line(aes(x=k.rng,y=mean.jac/mean.jac.scramble)) + theme_bw() + scale_y_continuous(limits=c(0,5)) + theme_bw() +
  ylab('Mean Jaccard Index') + xlab('k') + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(savedir,'PAM_ScrambledNormalizedSubsampleMeanJaccard_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 4,width=4,units='cm')

k.cdf.scramble <- lapply(partition.sim.mat.scramble[-1], function(X) ecdf(X)(seq(0,1,length.out = 100)))
k.auc.scramble <- sapply(k.cdf.scramble, mean)
p <- ggplot() + geom_line(aes(x=k.rng,y=(1-k.auc)/(1-k.auc.scramble))) + theme_bw() + scale_y_continuous(limits=c(0,5)) + theme_bw() +
  ylab('1- AUC(Jaccard CDF)') + xlab('k') + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(savedir,'PAM_ScrambledNormalizedSubsampleAUC_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 4,width=4,units='cm')
