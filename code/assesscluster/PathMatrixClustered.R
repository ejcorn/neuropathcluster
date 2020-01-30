rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleABC.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
microSample.comp <- complete(microSample.imp,1) # temporary
colnames(microSample.comp) <- colnames(microSample) # temporary
microSample <- microSample.comp

INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

##############################################
### plot pathology matrix after clustering ###
##############################################

microSample <- micro.order.by(microSample,vers='short') # order by feature
idx <- order(partition)
partition.order <- partition[idx]
midpoints <- sapply(1:k,function(k.i) floor(median(which(partition.order==k.i)))) # get midpoints of each stretch of subjects per cluster
rnames <- rownames(microSample[idx,]) # store rownames of reordered matrix
p <- imagesc(microSample[idx,],cmap = 'Blues',clim=c(1,5),caxis_labels = c('0','Rare','1+','2+','3+')) + 
  scale_y_discrete(limits= rev(rnames),breaks=rnames[midpoints],labels=paste('Cluster',1:k))+
  theme(axis.text.x = element_text(size=6,angle=90,hjust=1,vjust=0.5),axis.text.y = element_text(hjust=0.5,angle=90,color = getClusterColors(k)),
        legend.key.width = unit(0.1,'cm'),legend.background = element_blank(),legend.margin = ggplot2::margin(0,0,0,0))

ggsave(p,filename = paste0(savedir,'PathologyMatrixClustered.pdf'),
       units= 'cm',height = 12,width=18.5)
