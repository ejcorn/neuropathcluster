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

# traditional disease label
list[patientSample,dz.short]<- other.dz(patientSample)
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.short))

load(file = paste(savedir,'predictdz_Lassoperf',extralab,'.RData',sep=''))

met <- lapply(dz.res, function(a) sapply(a, function(c) c$Sensitivity))
p.sen <- plot.model.perf.met(met=met,perf.met='Sensitivity',ttl='Histopathologic Dx',colors = dz.colors[dz.short %in% names(met)])
met <- lapply(dz.res, function(a) sapply(a, function(c) c$Specificity))
p.spec <- plot.model.perf.met(met=met,perf.met='Specificity',ttl='Histopathologic Dx',colors = dz.colors[dz.short %in% names(met)])
met <- lapply(dz.res, function(a) sapply(a, function(c) c$AUC))
p.auc <- plot.model.perf.met(met=met,perf.met='AUC',ttl='Histopathologic Dx',colors = dz.colors[dz.short %in% names(met)])
met <- lapply(dz.res, function(a) sapply(a, function(c) c$OverallAccuracy))
p.acc <- plot.model.perf.met(met=met,perf.met='Accuracy',ttl='Histopathologic Dx',colors = dz.colors[dz.short %in% names(met)])
p.roc <- plot.model.roc(met=dz.res,ttl='Histopathologic Dx',colors=dz.colors[dz.short %in% names(dz.res)])
p.all <- plot_grid(plotlist = 
	list(p.sen,remove.y.ticklabels(p.spec),remove.y.ticklabels(p.auc),p.roc),
		align = 'hv',nrow=1,axis = 'b',rel_widths = c(1,1,1,1))

ggsave(filename = paste(savedir,'LassoPerformanceTraditionalDisease',extralab,'.pdf',sep=''),plot = p.all,
       height = 5,width=18,units='cm')

# cluster label
load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
clusterColors <- getClusterColors(k)

load(file = paste(savedir,'predictcluster_Lassoperf',extralab,'.RData',sep=''))

met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Sensitivity))
p.sen <- plot.model.perf.met(met=met,perf.met='Sensitivity',ttl='Data-driven',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Specificity))
p.spec <- plot.model.perf.met(met=met,perf.met='Specificity',ttl='Data-driven',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$AUC))
p.auc <- plot.model.perf.met(met=met,perf.met='AUC',ttl='Data-driven',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$OverallAccuracy))
p.acc <- plot.model.perf.met(met=met,perf.met='Accuracy',ttl='Data-driven',colors = clusterColors)
p.roc <- plot.model.roc(met=cluster.res,ttl='Data-driven',colors=clusterColors)
p.all <- plot_grid(plotlist = 
	list(p.sen,remove.y.ticklabels(p.spec),remove.y.ticklabels(p.auc),p.roc),
		align = 'hv',nrow=1,axis = 'b',rel_widths = c(1,1,1,1))


ggsave(filename = paste(savedir,'LassoPerformanceClusters',extralab,'.pdf',sep=''),plot = p.all,
       height = 5,width=18,units='cm')
