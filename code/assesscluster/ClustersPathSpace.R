# visualize clusters in pathology space using exploratory factor analysis of
# polychoric correlation matrix + UMAP

rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/pathspace/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

#####################################
### plot clusters in factor space ###
#####################################

load(file = paste0(params$opdir,'micekmeans/FactorAnalysis.RData')) # load factor analysis
scores <- scores[-DisconnectedSubjects,]
clusterColors <- getClusterColors(k)

list[patientSample,dz.short]<- other.dz(patientSample)
dz.names <- sort(unique(patientSample$NPDx1)) # get all disease names
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.names)) # assign a color to each disease
names(dz.colors) <- dz.names

df <- data.frame(x=scores[,1],y=scores[,2],z=scores[,3],g=partition,d=patientSample$NPDx1)
asp.ratio <- 0.75
p <- ggplot(df) + geom_point(aes(x=x,y=y,color=as.character(g)),stroke=0,alpha=0.5) + theme_classic() +
  xlab('Component 1') + ylab('Component 2') +
  scale_color_manual(values=clusterColors,name='Cluster')+ coord_fixed(asp.ratio)
ggsave(plot = p,filename = paste0(savedir,'ClustersInFactorSpace2D.pdf'),
       width = 8, height = 8, units = "cm")

p <- ggplot(df) + geom_point(aes(x=x,y=y,color=as.character(d)),stroke=0,alpha=0.5) + theme_classic() +
  xlab('Component 1') + ylab('Component 2') +
  scale_color_manual(values=dz.colors,name='Disease') + theme(legend.key.size = unit(0.5,'cm'),legend.text = element_text(size=4))+ 
  coord_fixed(asp.ratio)
ggsave(plot = p,filename = paste0(savedir,'DiseasesInFactorSpace2D.pdf'),
       width = 18, height = 8, units = "cm")

pdf(paste0(savedir,'ClustersInFactorSpace3D.pdf'),width=3,height = 3)
scatterplot3d(x=df$x,y=df$y,z=df$z, pch = 20, color=clusterColors[partition],cex.symbols = 0.2,cex.axis = 0.25,cex.lab = 0.5)
legend(x = max(df$x),y=max(df$y), legend = levels(factor(unique(as.character(partition)))),
       col = clusterColors, pch = 20,cex=0.5,pt.cex = 1)
dev.off()

pdf(paste0(savedir,'DiseasesInFactorSpace3D.pdf'),width=3,height = 3)
scatterplot3d(x=df$x,y=df$y,z=df$z, pch = 20, color=dz.colors[patientSample$NPDx1],cex.symbols = 0.2,cex.axis = 0.25,cex.lab = 0.5)
legend(x = max(df$x),y=max(df$y), legend = dz.names,
       col = dz.colors, pch = 20,cex=0.5,pt.cex = 1)
dev.off()

#############################################
### plot clusters in UMAP component space ###
#############################################

library(umap)
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
microSample.comp <- complete(microSample.imp,1) # temporary
colnames(microSample.comp) <- colnames(microSample) # temporary
microSample.comp <- remove.Disconnected.Subjects(microSample.comp,DisconnectedSubjects)
path.umap <- umap(microSample.comp,metric='pearson',n_components=3)

t.umap <- umap(as.data.frame(matrix(rnorm(1000),nrow=100,ncol=10)),metric='pearson',n_components=3)
plot(t.umap$layout[,1],t.umap$layout[,2])
df <- data.frame(x=path.umap$layout[,1],y=path.umap$layout[,2],z=path.umap$layout[,3],g=partition,d=patientSample$NPDx1)

p <- ggplot(df) + geom_point(aes(x=x,y=y,color=as.character(g)),stroke=0,alpha=0.5) + theme_classic() +
  xlab('Component 1') + ylab('Component 2') +
  scale_color_manual(values=clusterColors,name='Cluster')+ 
  coord_fixed(asp.ratio)
ggsave(plot = p,filename = paste0(savedir,'ClustersInUMAPSpace2D.pdf'),
       width = 8, height = 8, units = "cm")

p <- ggplot(df) + geom_point(aes(x=x,y=y,color=as.character(d)),stroke=0,alpha=0.5) + theme_classic() +
  xlab('Component 1') + ylab('Component 2') +
  scale_color_manual(values=dz.colors,name='Disease') + theme(legend.key.size = unit(0.5,'cm'),legend.text = element_text(size=4))+ 
  coord_fixed(asp.ratio)
ggsave(plot = p,filename = paste0(savedir,'DiseasesInUMAPSpace2D.pdf'),
       width = 18, height = 8, units = "cm")

pdf(paste0(savedir,'ClustersInUMAPSpace3D.pdf'),width=3,height = 3,useDingbats = F)
scatterplot3d(x=df$x,y=df$y,z=df$z, pch = 20, color=clusterColors[partition],cex.symbols = 0.2,cex.axis = 0.25,cex.lab = 0.5)
legend(x = max(df$x),y=max(df$y), legend = paste('Cluster',1:k),
       col = clusterColors, pch = 20,cex=0.5,pt.cex = 1)
dev.off()

pdf(paste0(savedir,'DiseasesInUMAPSpace3D.pdf'),width=3,height = 3,useDingbats = F)
scatterplot3d(x=df$x,y=df$y,z=df$z, pch = 20, color=dz.colors[patientSample$NPDx1],cex.symbols = 0.2,cex.axis = 0.25,cex.lab = 0.5)
legend(x = max(df$x),y=max(df$y), legend = dz.names,
       col = dz.colors, pch = 20,cex=0.5,pt.cex = 1)
dev.off()

