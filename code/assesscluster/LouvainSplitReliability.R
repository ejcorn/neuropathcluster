rm(list = setdiff(ls(),c('params','sampfrac')))
savedir <- paste(params$opdir,'optimcluster/',sep='')
source('code/misc/fxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

########################################
### Get plotting order for centroids ###
########################################

pathItems.type <- fliplr(list("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody"))
pathItems.index <- sapply(1:length(pathItems.type), function(i) grep(pathItems.type[[i]], colnames(microSample)))
pathItems.labels <- sapply(1:length(pathItems.type), function(i) c(matrix("",floor(0.5*length(pathItems.index[[i]]))),
                                                                   pathItems.type[[i]], c(matrix("",ceiling(0.5*length(pathItems.index[[i]])-1)))))
pathItems.index <- Reduce(c,pathItems.index)
pathItems.labels <- Reduce(c,pathItems.labels)

microSample <- microSample[,pathItems.index]

# load full centroids

e1 <- new.env()
load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
centroids.full <- get('centroids', e1)

##################################
### Load subsampled partitions ###
##################################

cl.data <- readMat(paste(savedir,'subsamplePartitionsSF',sampfrac,'.mat',sep=''))

partitions <- cl.data$partitions
samples <- cl.data$samples
n_samps <- dim(partitions)[2]

############################################
### Compute centroids for each partition ###
############################################

centroids <- centroid.orig.cor <- list()
cluster.order <- c('Tau','Thio','TDP43','Syn')
k.all <- rep(NA,n_samps)
for(S in 1:n_samps){ # iterate through each subsample
  samp <- samples[,S] # get indices of subjects
  microSample.subsample <- microSample[samp,] # make subsampled data matrix
  partition <- partitions[,S] # get partition for given subsample
  # sometimes clusters are 1,2,3,5 etc. bc multislice pair w/ variability so get new k
  # compute centroids
  centroids[[S]] <- compute.centroids(microSample.subsample,partition)
  # annex small clusters to nearest neighbor (measured by spearman cor)
  partition <- annex.small.clusters(partition,microSample.subsample,centroids[[S]])
  partition <- remove.Partition.Gap(partition)
  # recompute centroids after annexing small clusters
  centroids[[S]] <- compute.centroids(microSample.subsample,partition)  
  # reorder centroids based on correlations with original centroids
  centroid.orig.cor[[S]] <- cor(t(centroids.full),t(centroids[[S]]),use='pairwise.complete.obs')
  # only attempt to reorder if unique indices & at least 4 cluster
  shuffIdx <- row.Which.Max(centroid.orig.cor[[S]])
  if(length(unique(shuffIdx)) == length(shuffIdx) & length(shuffIdx) == nrow(centroids.full)){    
    k.all[S] <- max(partition)
    centroids[[S]] <- centroids[[S]][shuffIdx,] 
  }
}

##################################################################
### Compute similarity matrices between subsample and original ###
##################################################################

sub.vs.orig.cors <- lapply(centroids, function(C) cor(t(centroids.full),t(C),use='pairwise.complete.obs'))
sub.vs.orig.cors <- lapply(sub.vs.orig.cors, function(x) x[diag(TRUE,nrow=nrow(x),ncol=ncol(x))])
sub.vs.orig.cors <- unlist(sub.vs.orig.cors)

############################################
### Plot histogram of correlation values ###
############################################

h <- hist(sub.vs.orig.cors)
binned.counts <- unlist(mapply(function(mid,count) rep(mid,count),mid=h$mids,count=h$counts))
df <- data.frame(x=rep(1,length(sub.vs.orig.cors)),
                 y=factor(binned.counts))

# set up breaks at 0,25,50,75,100%
brks <- round(sum(h$counts) * c(0,0.25,0.5,0.75,1))
save(brks,df,file=paste(savedir,'FigS2b_SourceData.RData',sep=''))
p <- ggplot(df) + 
  #geom_bar(aes(x=y,y=(..count..)/sum(..count..)),fill='red',alpha=0.6) +
  geom_histogram(aes(x=sub.vs.orig.cors),fill='red',alpha=0.6)+ 
  scale_y_continuous(limits=c(min(h$counts),sum(h$counts)),breaks=brks, labels = paste(signif(100*brks/sum(h$counts),2),'%',sep=''))+
  #scale_y_continuous(labels=percent) +
  xlab('Pearson\'s r') + ylab('Percent of Total') + ggtitle('Similarity to Full Sample Centroids') +
  theme_classic() + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5,size=8,face='bold'))
p

ggsave(plot = p,filename = paste(savedir,'SubsampleSimilarityToOriginal_SF',sampfrac,'.pdf',sep=''),width = 2.5, height = 2, units = "in")

#####################################################
### Compute similarity matrices between centroids ###
#####################################################

c.2 <- combn(x=1:n_samps,m=2) # unique combinations of samples
intersample.centroid.cors <- 
  lapply(as.data.frame(c.2), function(i) 
    cor(t(centroids[[i[1]]]),t(centroids[[i[2]]]),use='pairwise.complete.obs'))

intersample.centroid.cors <- 
  lapply(intersample.centroid.cors, function(x) x[diag(TRUE,nrow=nrow(x),ncol=ncol(x))])
intersample.centroid.cors <- unlist(unlist(intersample.centroid.cors))
 

