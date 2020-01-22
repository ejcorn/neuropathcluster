rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')
source('code/misc/plottingfxns.R')
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs

##################
### Make plots ###
##################

load(file = paste(savedir,'predictcluster_RFperf',extralab,'.RData',sep=''))
clusterColors <- getClusterColors(length(cluster.res))

p1 <- plot.featureweights.rf(cluster.res,clusterColors)
ggsave(filename = paste(savedir,'RFClusterFeatureWeights',extralab,'.pdf',sep=''),plot = p1,
       height = 5,width=7,units='cm')

# traditional disease label
load(file = paste(savedir,'predictdz_RFperf',extralab,'.RData',sep=''))
list[patientSample,dz.short]<- other.dz(patientSample)
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.short))

p2 <- plot.featureweights.rf(dz.res,dz.colors[dz.short %in% names(dz.res)])
ggsave(filename = paste(savedir,'RFTraditionalDiseaseFeatureWeights',extralab,'.pdf',sep=''),plot = p2,
       height = 5,width=7,units='cm')