rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')
source('code/misc/plottingfxns.R')
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

##################
### Make plots ###
##################

load(file = paste(savedir,'predictcluster_GLMperf',extralab,'.RData',sep=''))
clusterColors <- getClusterColors(length(cluster.res))

p1 <- plot.featureweights.lm.contintsplit(cluster.res,df,clusterColors)

w.multiplier <- 1
if(extralab == 'CSFGene'){w.multiplier <- 2} # make plot wider if both csf and gene betas 
ggsave(filename = paste(savedir,'GLMClusterFeatureWeights',extralab,'.pdf',sep=''),plot = p1,
       height = 5,width=5*w.multiplier,units='cm')

# traditional disease label
load(file = paste(savedir,'predictdz_GLMperf',extralab,'.RData',sep=''))
list[patientSample,dz.short]<- other.dz(patientSample)
if(extralab == 'CSFOnly' | extralab == 'CSFGene'){
	pal <- colorRampPalette(brewer.pal(name = 'Blues',n=9))
	dz.colors <- pal(9)[c(4,6,8)] # pick good range of blues
}
if(extralab == 'GeneOnly'){
	pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
	dz.colors <- pal(length(dz.short)) # assign a color to each disease
	# select only diseases included in analysis
	dz.colors <- dz.colors[dz.short %in% names(dz.res)]
}

w.multiplier <- 1
if(extralab == 'CSFGene'){w.multiplier <- 2} # make plot wider if both csf and gene betas 
p2 <- plot.featureweights.lm.contintsplit(dz.res,df,dz.colors)
ggsave(filename = paste(savedir,'GLMTraditionalDiseaseFeatureWeights',extralab,'.pdf',sep=''),plot = p2,
       height = 5,width=5*w.multiplier,units='cm')