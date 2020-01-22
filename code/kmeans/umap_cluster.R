# cluster subjects based on UMAP component loadings

rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/pathspace/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(savedir,'microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

#############################################################
### cluster subjects basesd on weights on UMAP components ###
#############################################################

clusterColors <- getClusterColors(k)
dz.names <- sort(unique(patientSample$NPDx1)) # get all disease names

pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.names)) # assign a color to each disease
names(dz.colors) <- dz.names

library(umap)
library(cluster)
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
microSample.comp <- complete(microSample.imp,1) # temporary
colnames(microSample.comp) <- colnames(microSample) # temporary
microSample.comp <- remove.Disconnected.Subjects(microSample.comp,DisconnectedSubjects)
path.umap <- umap(microSample.comp,metric='pearson',n_components=3)

k.rng <- 2:15
X <- path.umap$layout
clust <- lapply(k.rng, function(k) kmeans(x = X,centers = k,nstart = 100))
sil <- lapply(clust, function(K) silhouette(x=K$cluster,dist=dist(X)))
sil.mean <- sapply(sil,function(S) mean(S[,3]))
k.idx <- which.max(sil.mean)
k <- k.rng[k.idx]
partition.umap <- clust[[k.idx]]$cluster
clusterColors <- getClusterColors(k)

p <- ggplot() + geom_line(aes(x=k.rng,y=sil.mean)) + #theme_classic() +
  scale_x_continuous(breaks = k.rng) + xlab('k') + ylab('Mean Silhouette')
ggsave(filename = paste0(savedir,'UMAP_KMeansSilhouetteByK.pdf'),plot = p,
       height = 8,width=8,units='cm')

df <- data.frame(x=X[,1],y=X[,2],z=X[,3],g=partition.umap)
pdf(paste0(savedir,'UMAP_KmeansInUMAPSpace3D.pdf'),width=3,height = 3)
scatterplot3d(x=df$x,y=df$y,z=df$z, pch = 20, color=clusterColors[partition.umap],cex.symbols = 0.2,cex.axis = 0.25,cex.lab = 0.5)
legend(x = max(df$x),y=max(df$y), legend = levels(factor(unique(as.character(partition)))),
       col = clusterColors, pch = 20,cex=0.5,pt.cex = 1)
dev.off()

list[patientSample,dz.short]<- other.dz(patientSample)

list[p.k.dz,df.plt] <- plot.dz.by.cluster(patientSample$NPDx1,partition.umap,dz.short,clusterColors,'Primary\nHistopathologic Diagnosis')
ggsave(filename = paste(savedir,'UMAP_ClustersByPrimaryDisease.pdf',sep=''),plot = p.k.dz,
       height = 2,width=3,units='in')
