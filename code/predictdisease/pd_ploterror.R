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

# traditional disease label
list[patientSample,dz.short]<- other.dz(patientSample)
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.short))

load(file = paste(savedir,'predictdz_RFperf',extralab,'.RData',sep=''))

for(c.res in cluster.res){
	for(c.rep in c.res){
		
	}
}