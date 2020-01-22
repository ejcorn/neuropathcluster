rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleBraakCERAD.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,1]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

############################
### Plot distance matrix ###
############################

cl.data <- readMat(paste(params$opdir,'optimcluster/subjectCorrMat.mat',sep=''))
load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
W <- cl.data$W
if(params$dist.met == 'polychoric'){ 
  # if using polychoric correlations -> pam clustering, then threshold correlation matrix
  # best threshold is saved
  #W <- thresh.mat(W,'>',thresh.bestsil)
}

s <- sample(nrow(W),replace = F)
melted_cormat <- melt(W[s,s])
melted_cormat$Var1 <- melted_cormat$Var1[length(melted_cormat$Var1):1]

pal <- colorRampPalette(c('dark red','#c23b22','white','#779ecb','dark blue'))
pal <- pal(100)
p1 <- ggplot() + 
  geom_tile(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
  scale_y_discrete(labels = NULL,expand=c(0,0)) + 
  scale_x_discrete(labels = NULL,expand=c(0,0)) + coord_fixed() +
  scale_fill_gradientn(name = "r",limits = c(-1,1),breaks=c(-1,0,1),colours = pal) + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size = 6,color='black'),
        axis.text.y = element_text(size = 6,color='black'), text = element_text(size = 6,color='black'),
        axis.ticks = element_blank(),axis.line = element_blank(),
        legend.position = 'none',
        #legend.key.height = unit(0.15,'in'),legend.key.width = unit(0.05,'in')
        )
#ggsave(p1,filename = 'results/subjectCorrShuffled.png',units = 'in',height = 3,width = 3)

########################
### Order by cluster ###
########################

k <- max(partition)
k.Index <- lapply(1:k, function(k.i) which(partition == k.i)) 
k.Labels <- lapply(1:k, function(i) c(matrix("",floor(0.5*length(k.Index[[i]]))), # make axis labels
              paste('Cluster',i), c(matrix("",ceiling(0.5*length(k.Index[[i]])-1)))))
k.Labels <- Reduce(c,k.Labels)

k.rect.idx <- fliplr(sapply(1:k, function(k.i) sum(partition == k.i)))
cs.rect <- cumsum(k.rect.idx)
df <- data.frame(xmx = 1+cs.rect, xmi = c(1,1+cs.rect[-length(cs.rect)]))

pal.k <- colorRampPalette(brewer.pal(name = 'Dark2',n=8))
clusterColors <- pal.k(8)[(1:k)+2]

melted_cormat <- melt(W[order(partition),order(partition)])
melted_cormat$Var1 <- melted_cormat$Var1[length(melted_cormat$Var1):1]
p2 <- ggplot() + 
  geom_tile(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
  # boxes on diagonal only
  geom_rect(data=df,aes(xmin = xmi, xmax = xmx, ymin = nrow(W)-xmi+1, ymax = nrow(W)-xmx +1),fill=NA,color='black',size=0.4) +
  scale_y_discrete(limits=1:length(k.Labels),labels = k.Labels,expand=c(0,0)) + 
  scale_x_discrete(limits=1:length(k.Labels),labels = fliplr(k.Labels),expand=c(0,0)) + coord_fixed() +
  scale_fill_gradientn(name = "r",limits = c(-1,1),breaks=c(-1,0,1),colours = pal) + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size = 6,color='black'),
        axis.text.y = element_text(size = 6,color='black'), text = element_text(size = 6,color='black'),
        axis.ticks = element_blank(),axis.line = element_blank(),
        legend.position = 'none'
        #legend.key.height = unit(0.15,'in'),legend.key.width = unit(0.05,'in')
        )
#ggsave(p2,filename = 'results/subjectCorrLouvain.png',units = 'in',height = 3,width = 3)

#######################################
### Order by disease for comparison ###
#######################################

list[patientSample,dz.short]<- other.dz(patientSample)

dz.Names <- fliplr(sort(unique(as.character(patientSample$NPDx1))))
dz.short <- fliplr(dz.short)
list[dz.Index,dz.Labels]<- get.feature.labels(dz.Names,patientSample$NPDx1)
for(i in 1:length(dz.Names)){dz.Labels[which(dz.Labels == dz.Names[i])] <- dz.short[i]} # replace long disease names with short labels
# 
# dz.Index <- lapply(dz.Names, function(name.i) grep(name.i, patientSample$NPDx1))
# dz.Labels <- lapply(1:length(dz.Names), function(i) c(matrix("",floor(0.5*length(dz.Index[[i]]))),
#               dz.short[[i]], c(matrix("",ceiling(0.5*length(dz.Index[[i]])-1)))))
# dz.Index <- Reduce(c,dz.Index)
# dz.Labels <- Reduce(c,dz.Labels)

dz.rect.idx <- unlist(lapply(dz.Names, function(x) 
  length(grep(x,patientSample$NPDx1[dz.Index]))))
#dz.rect.idx <- fliplr(dz.rect.idx)
cs.rect <- cumsum(dz.rect.idx)
df <- data.frame(xmx = 1+cs.rect, xmi = c(1,1+cs.rect[-length(cs.rect)]))

pal.dz <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal.dz(length(dz.Names))

melted_cormat <- melt(W[dz.Index,dz.Index])
melted_cormat$Var1 <- melted_cormat$Var1[length(melted_cormat$Var1):1]
p3 <- ggplot() + 
  geom_tile(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
  # boxes on diagonal only
  geom_rect(data=df,aes(xmin = xmi, xmax = xmx, ymin = nrow(W)-xmi+1, ymax = nrow(W)-xmx +1),fill=NA,color='black',size=0.4) +
  scale_y_discrete(limits=1:length(dz.Labels),labels = dz.Labels,expand=c(0,0)) + 
  scale_x_discrete(limits=1:length(dz.Labels),labels = fliplr(dz.Labels),expand=c(0,0)) + coord_fixed() +
  scale_fill_gradientn(name = "r",limits = c(-1,1),breaks=c(-1,0,1),colours = pal) + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size = 6,color='black'),
        axis.text.y = element_text(size = 6,color='black'), text = element_text(size = 6,color='black'),
        axis.ticks = element_blank(),axis.line = element_blank(),
        legend.position = 'none'
        #legend.key.height = unit(0.15,'in'),legend.key.width = unit(0.05,'in')
        )
#ggsave(p3,filename = 'results/subjectCorrDisease.png',units = 'in',height = 3,width = 3)

p1p2p3 <- plot_grid(plotlist=list(p1,p3,p2),align = 'hv',nrow=1)
ggsave(p1p2p3,filename = paste(savedir,'subjectCorrMats.png',sep=''),
       units = 'in',height = 8/3,width = 8)

save(W,dz.Names,partition,file = paste(savedir,'Figure2a-c_SourceData.RData',sep=''))
