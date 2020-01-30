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

#####################################################
### see if you get same centroids in split halves ###
#####################################################

# Take 60% subsamples. Cluster at range of k and thresholds.
# 1.	Plot distribution of silhouette values as a function of k and threshold. 
# 2.	For each subsample, figure out best threshold (by peak silhouette.
# 3.	For each subsample, attempt to match centroids to full sample, and then plot distribution of agreement (for each cluster or aggregated).

sampfrac <- 1
#k.rng <- 2:10
#thresh.rng <- seq(0,0.5,0.1)
load(file = paste0(savedir,'ClusterBootstrapping_SF',sampfrac,'.RData'))
load(file = paste0(savedir,'ScrambledClusterBootstrapping_SF',sampfrac,'.RData'))

# 1.	Plot distribution of silhouette values as a function of k and threshold. 

# loop through all clustering on bootstrapped data, pull out silhouette values for each k and threshold combo
df.sub <- lapply(thresh.rng,function(thresh) do.call('rbind',lapply(k.rng, function(k) do.call('rbind',lapply(1:nreps, function(R)
  data.frame(sil=clust.thresh[[R]][[which(thresh.rng==thresh)]][[which(k.rng==k)]]$silinfo$avg.width,
             k=k,thresh=thresh,rep=R,null.or.data='Data'))))))
names(df.sub) <- as.character(thresh.rng)

df.sub.scramble <- lapply(thresh.rng,function(thresh) do.call('rbind',lapply(k.rng, function(k) do.call('rbind',lapply(1:nreps, function(R)
  data.frame(sil=clust.thresh.scramble[[R]][[which(thresh.rng==thresh)]][[which(k.rng==k)]]$silinfo$avg.width,
             k=k,thresh=thresh,rep=R,null.or.data='Null'))))))
names(df.sub.scramble) <- as.character(thresh.rng)

max.sil <- max(do.call('rbind',df.sub)$sil) # scale y axis of all plots the same -- based on data
p.list.silhouette <- lapply(as.character(thresh.rng),function(t) ggplot() + #geom_line(data=df.full[[t]],aes(x=k,y=sil)) + #theme_classic() +
                   geom_boxplot(data=rbind(df.sub[[t]],df.sub.scramble[[t]]),aes(x=as.character(k),y=sil,fill=null.or.data),size=0.25,outlier.size = 0.5)+
                   scale_x_discrete(limits=as.character(k.rng),breaks = as.character(k.rng)) + 
                   scale_y_continuous(limits=c(-0.05,max.sil)) + scale_fill_manual(values=c("#C6CDF7",'grey70'))+
                   xlab('k') + ylab('Mean Silhouette') + ggtitle(paste('Threshold =',t))+ theme_bw()+
                   theme(text=element_text(size=8),legend.position = 'none')+
                   theme(plot.title = element_text(hjust=0.5,size=8),axis.text.x = element_text(size=6),text=element_text(color='black')))
p.silhouette <- plot_grid(plotlist = p.list.silhouette)
ggsave(filename = paste0(savedir,'PAM_ScramblevsDataBootstrapMeanSilhouetteByKByThresh_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p.silhouette,
       height = 8,width=14,units='cm')

# 2.	For each Bootstrap, figure out best threshold (by peak silhouette.

df.sub.all <- do.call('rbind',df.sub)
silhouette.by.thresh  <- sapply(1:nreps, function(R) 
  sapply(thresh.rng, function(thresh) mean(df.sub.all[df.sub.all$rep==R & df.sub.all$thresh==thresh,'sil'])))
thresh.indicator <- sapply(1:nreps, function(R) sapply(thresh.rng, function(thresh) thresh)) # for plot, store threshold values
df.plot <- data.frame(sil=as.vector(silhouette.by.thresh),thresh=as.vector(thresh.indicator))
p.mean.sil <- ggplot(df.plot) + geom_boxplot(aes(x=as.character(thresh),y=sil),fill="#7294D4") + theme_classic() +
  ggtitle('')+theme(text=element_text(size=8))+ scale_y_continuous(limits=c(-0.05,max(silhouette.by.thresh)))+
  ylab('Mean Silhouette') + xlab('Threshold')
ggsave(filename = paste0(savedir,'PAM_BootstrapMeanSilhouetteByThresh_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p.mean.sil,
       height = 4,width=4,units='cm')

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

mean.jac <- sapply(partition.sim.mat[-1], function(X) mean(X*as.numeric(!diag(nreps)))) # ignore diagonals, remove k =1
p <- ggplot() + geom_line(aes(x=k.rng,y=mean.jac),color="#7294D4") + theme_bw() + scale_y_continuous(limits=c(0,5)) +
  ylab('Mean Jaccard Index') + xlab('k') + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(savedir,'PAM_BootstrapMeanJaccard_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 4,width=4,units='cm')

k.cdf <- lapply(partition.sim.mat[-1], function(X) ecdf(X)(seq(0,1,length.out = 100)))
k.auc <- sapply(k.cdf, mean)
p <- ggplot() + geom_line(aes(x=k.rng,y=(1-k.auc)),color="#7294D4") + theme_bw() + scale_y_continuous(limits=c(0,5)) +
  ylab('1 - AUC(Jaccard CDF)') + xlab('k') + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(savedir,'PAM_BootstrapAUC_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 4,width=4,units='cm')

mean.jac.scramble <- sapply(partition.sim.mat.scramble[-1], function(X) mean(X*as.numeric(!diag(nreps)))) # ignore diagonals, remove k=1

p.meanjac.norm <- ggplot() + geom_line(aes(x=k.rng,y=mean.jac/mean.jac.scramble),color="#7294D4") + theme_bw() + scale_y_continuous(limits=c(0,6)) +
  ggtitle(paste('Threshold =',thresh.bestsil)) + ylab('Mean Jaccard Index') + xlab('k') + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5,size=8))
ggsave(filename = paste0(savedir,'PAM_ScrambledNormalizedBootstrapMeanJaccard_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 4,width=4,units='cm')

k.cdf.scramble <- lapply(partition.sim.mat.scramble[-1], function(X) ecdf(X)(seq(0,1,length.out = 100)))
k.auc.scramble <- sapply(k.cdf.scramble, mean)
p.auc.norm <- ggplot() + geom_line(aes(x=k.rng,y=(1-k.auc)/(1-k.auc.scramble)),color="#7294D4") + theme_bw() + scale_y_continuous(limits=c(0,6)) +
  ggtitle(paste('Threshold =',thresh.bestsil)) + ylab('1 - AUC(Jaccard CDF)') + xlab('k') + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(savedir,'PAM_ScrambledNormalizedBootstrapAUC_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p,
       height = 4,width=4,units='cm')

# align all plots into one figure
p.right <- plot_grid(plotlist=list(p.mean.sil,p.auc.norm),ncol=1)
p.all <- plot_grid(plotlist=list(p.silhouette,p.right),rel_widths = c(3,1),nrow=1)
ggsave(filename = paste0(savedir,'FigS2_',params$dist.met,'SF',sampfrac,'.pdf'),plot = p.all,
       height = 8,width=18,units='cm')
