rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'diseasenetwork/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

########################################################
### Construct a network of overlap between diagnoses ###
########################################################

list[patientSample.othered,dz.short] <- other.dz(patientSample)

dzs.othered <- sort(unique(patientSample.othered$NPDx1))
dzs.othered <- dzs.othered[!dzs.othered %in% c('Other','Unremarkable adult','PART')]
codx.network <- matrix(NA,nrow=length(dzs.othered),ncol=length(dzs.othered),dimnames=list(dzs.othered,dzs.othered))

for(dz1 in dzs.othered){
  dz1.mask <- !exclude.dz(patientSample.othered,dz1,5)
  for(dz2 in dzs.othered){
    if(dz1 != dz2){
      dz2.mask <- !exclude.dz(patientSample.othered,dz2,5)
      codx.network[dz1,dz2] <- sum(dz1.mask & dz2.mask)/ sum(dz1.mask) # proportion of dz1 that is + for dz2
    }
  }
}

p <- imagesc(codx.network,caxis_name = 'Pr(Dz2 | Dz1)', xlabel = 'Disease 2',ylabel = 'Disease 1') + 
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))
p
library(igraph)
plot(graph_from_adjacency_matrix(adjmatrix = codx.network,weighted = TRUE,mode = 'directed',diag = FALSE))

