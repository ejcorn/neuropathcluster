rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/patientcharacteristics/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleABC.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics122219.csv',sep=''))
onset <- read.csv(paste0('data/INDD_GlobalOnset122219.csv'),stringsAsFactors = F)

## subset onset data to our sample, in same order
rownames(onset) <- as.character(onset$INDDID)
onset <- onset[as.character(INDDIDs),]

## process patientSample to group patients with > low ADNC
AD.status.oi <- c('Intermediate','High')
EOAD.cutoff <- 65
patientSample <- cbind(patientSample,GlobalAgeOnset=onset$GlobalAgeOnset,Cluster=paste('Cluster',partition),EOAD=onset$GlobalAgeOnset<EOAD.cutoff)
patientSample <- patientSample[patientSample$ADStatus %in% AD.status.oi,]

boot.prop <- function(x1,x2,nperms=10000){
  # INPUTS:
  # x1: binary vector
  # x2: binary vector
  #
  # OUTPUTS:
  # bootstrapped p value for difference in proportion
  x1.prop <- sapply(1:nperms, function(j) mean(sample(x=x1,size=length(x1),replace = T),na.rm=T))
  x2.prop <- sapply(1:nperms, function(j) mean(sample(x=x2,size=length(x2),replace = T),na.rm=T))
  pval <- pval.2tail.np(0,x1.prop-x2.prop) # test null hypothesis that difference in means = 0
  return(pval)
}

clusterNames <- paste('Cluster',1:k)
sapply(clusterNames, function(k1) mean(patientSample$EOAD[patientSample$Cluster == k1],na.rm=T))
p.vals <- sapply(clusterNames, function(k1) 
  sapply(clusterNames, function(k2) boot.prop(patientSample$EOAD[patientSample$Cluster == k1],patientSample$EOAD[patientSample$Cluster == k2])))
pvals <- p.adjust(p.vals,method = 'fdr')
diffs <- sapply(clusterNames, function(k1) 
  sapply(clusterNames, function(k2) mean(patientSample$EOAD[patientSample$Cluster == k1],na.rm=T) - mean(patientSample$EOAD[patientSample$Cluster == k2],na.rm=T)))

clusterColors <- getClusterColors(k)
p <- plot.allele.beta.matrix(diffs,p.vals,'Early Onset AD','',min(diffs),max(diffs),clusterColors)
ggsave(filename = paste(savedir,'EOADByCluster.pdf',sep=''),plot = p,
       height = 5,width=4.5,units='cm')
save(diffs,p.vals,file=paste0(params$sourcedatadir,'FigS6b_SourceData_EOADProportion.RData'))