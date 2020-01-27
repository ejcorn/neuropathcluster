rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1)] # Get rid of index column only

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)

###########
### Cog ###
###########

cog <- read.csv(paste(params$opdir,'processed/MoCA_processed.csv',sep=''),stringsAsFactors = F)[,-1] # remove index column
cog <- cog[cog$INDDID %in% INDDIDs,]
cog$TestDate <- as.Date(cog$TestDate,format='%m/%d/%Y')
cog.by.subj <- lapply(unique(cog$INDDID), function(ID) 
  cog[cog$INDDID == ID,]) # extract all tests for patient
cog.by.subj <- lapply(cog.by.subj, function(X) X[order(X$TestDate),]) # sort so more recent tests come last
cog.by.subj <- lapply(cog.by.subj, function(X) X[nrow(X),]) # only keep most recent test (closest to autopsy)
cog.by.subj <- do.call('rbind',cog.by.subj)

names(partition) <- as.character(INDDIDs)
df.all <- cbind(cog.by.subj, data.frame(Cluster=paste('Cluster',partition[as.character(cog.by.subj$INDDID)]),stringsAsFactors = F))

clusterNames <- sort(unique(df.all$Cluster))
comps <- lapply(as.data.frame(combn(clusterNames,m = 2)),function(x) as.character(x))

pairwise.wtests <- function(df.all){return(sapply(clusterNames, function(k1) sapply(clusterNames, function(k2) wilcox.test(df.all$VisuospatialTotal[df.all$Cluster == k1],df.all$VisuospatialTotal[df.all$Cluster == k2],conf.int=TRUE)$estimate)))}
wtests <- pairwise.wtests(df.all)
nreps <- 100
wtests.null <- lapply(1:nreps, function(R) pairwise.wtests(data.frame(VisuospatialTotal=df.all$VisuospatialTotal,Cluster = sample(df.all$Cluster,replace=F))))

df.plot1 <- data.frame(diff=as.vector(wtests),grp=rep('Actual',length(as.vector(wtests))))
df.plot2 <- data.frame(diff=unlist(wtests.null),grp=rep('Null',length(unlist(wtests.null))))
ggplot(rbind(df.plot1,df.plot2)) + geom_boxplot(aes(x=grp,y=abs(diff)))
