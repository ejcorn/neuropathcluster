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

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

##################
### Clustering ###
##################

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

k.rng <- 2:15
thresh.rng <- seq(0,0.5,0.1)
clust.thresh <- lapply(thresh.rng, function(thresh) lapply(k.rng, function(k)  # use 1 - r to make into dissimilarity
  pam(x = 1-thresh.mat(W,'>',thresh),diss = TRUE,k=k)))
sil.mean <- lapply(clust.thresh, function(clust) sapply(clust,function(K) K$silinfo$avg.width))
names(sil.mean)<- thresh.rng
p.list <- lapply(as.character(thresh.rng),function(t) ggplot() + geom_line(aes(x=k.rng,y=sil.mean[[t]])) + #theme_classic() +
  scale_x_continuous(breaks = k.rng) + scale_y_continuous(limits=c(-0.1,max(unlist(sil.mean)))) +
    xlab('k') + ylab('Mean Silhouette') + ggtitle(paste('Thresh.',t))+
    theme(plot.title = element_text(hjust=0.5)))
p <- plot_grid(plotlist = p.list)
ggsave(filename = paste0(savedir,'PAM_MeanSilhouetteByKByThresh_',params$dist.met,'.pdf'),plot = p,
       height = 18,width=18,units='cm')

# threshold value for which highest silhouette value is achieved at any value of k --
# this is a parameter that ought to be fit and cross-validated

thresh.bestsil.idx <- which.max(sapply(sil.mean,function(sil) max(sil)))
thresh.bestsil <- thresh.rng[thresh.bestsil.idx]
clust <- clust.thresh[[thresh.bestsil.idx]]

# compute medoids
for(k in 1:length(clust)){
  clust[[k]]$medoid.matrix <- microSample.imp[-DisconnectedSubjects,][as.numeric(clust[[k]]$medoids),]
}

p.list <- lapply(clust[which(k.rng <= 6)],function(X) imagesc(t(fliplr(micro.order.by(X$medoid.matrix)))) + 
                   theme(axis.text.y = element_text(size=4,hjust=1,vjust=0.5)))
p <- plot_grid(plotlist = p.list)
ggsave(filename = paste0(savedir,'PAMMedoids_',params$dist.met,'.pdf'),plot = p,
       height = 18,width=18,units='cm')

list[pathItems.type,pathRegions.name] <- get.pathscore.names(vers = 'short')
list[pathItems.index,pathItems.labels] <- get.feature.labels(pathItems.type,colnames(microSample.imp))
# order microSample columns to appear in Fig. 4d displaying cluster centroids
microSample <- microSample.imp[,pathItems.index]

# lists of features to automatically guide cluster reordering
reordering.features <- list('6'=list(c('Tau','CP_Thio','GP_Thio','TS_Thio'),
                                     c('Amyg_Tau','CS_Tau','Angiopathy'),c('TDP43'),
                                     c('SN_Gliosis','SN_NeuronLoss'),
                                     c('Syn','Tau'),unlist(pathItems.type)),
                                     '4'=c('CB_Tau','Thio','TDP43','Syn')) 

for(k in c(4,6)){
  partition <- unname(clust[[which(k.rng==k)]]$cluster)
  centroids <- compute.centroids(microSample[-DisconnectedSubjects,],partition)
  
  cluster.init.reorder <- order.cluster.by.feature.old(microSample,centroids,
                                                       reordering.features[[as.character(k)]],
                                                       max.min = c(rep('max',5),'min'))
  
  # only reorder if you can make a unique match for each cluster
  if(length(unique(cluster.init.reorder)) == length(cluster.init.reorder) & length(cluster.init.reorder) == length(unique(partition))){
    partition <- reorder.partition(partition,cluster.init.reorder)
  }
  
  centroids <- compute.centroids(microSample[-DisconnectedSubjects,],partition)
  
  # using medoids and calling them centroids so variable names match up for subsequent scripts
  #centroids <- as.matrix(microSample[-DisconnectedSubjects,][as.numeric(clust[[which(k.rng==k)]]$medoids),]) 
  rownames(centroids) <- sapply(1:k,function(k.i) paste('Cluster',k.i))
  # save such that gamma.opt parameter can index k
  savedir <- paste0(params$opdir,'results_G',k,'/analyzecluster/')
  dir.create(savedir,recursive = T)
  save(partition,centroids,k,DisconnectedSubjects,pathItems.labels,thresh.bestsil,centroids,
       file = paste0(savedir,'subjLouvainPartitionReordered.RData'))
}


