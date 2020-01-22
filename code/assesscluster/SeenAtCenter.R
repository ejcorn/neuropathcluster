## make plots for whether clusters are biased by patients in different centers that subserve the CDNR INDD
# idea being that maybe patients from certain centers were autopsied in a characteristic way
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleBraakCERAD.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics122219.csv',sep=''))
demo <- demo[demo$INDDID %in% INDDIDs,]

centers <- c('PTID','PDCID','FTDID','ALSID','CNDRID')
id.exist <- sapply(centers,function(s) !is.na(demo[,s]))
colnames(id.exist) <- centers
hist(rowSums(id.exist))

cluster.labels <- paste('Cluster',partition)
p.list <- lapply(centers, function(s) ggplot() + geom_bar(aes(x=cluster.labels[which(id.exist[,s])]))+
                   ggtitle(s) + theme(plot.title = element_text(hjust=0.5),axis.text.x=element_text(angle=90))+xlab(''))
plot_grid(plotlist = p.list)
