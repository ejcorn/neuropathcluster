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

met <- lapply(dz.res, function(a) sapply(a, function(c) c$Sensitivity))
p.sen.d <- plot.model.perf.met(met=met,perf.met='Test Sensitivity',ttl='',colors = dz.colors[dz.short %in% names(met)])
met <- lapply(dz.res, function(a) sapply(a, function(c) c$Specificity))
p.spec.d <- plot.model.perf.met(met=met,perf.met='Test Specificity',ttl='',colors = dz.colors[dz.short %in% names(met)])
met <- lapply(dz.res, function(a) sapply(a, function(c) c$AUC))
p.auc.d <- plot.model.perf.met(met=met,perf.met='Test AUC',ttl='',colors = dz.colors[dz.short %in% names(met)])
met <- lapply(dz.res, function(a) sapply(a, function(c) c$OverallAccuracy))
p.acc.d <- plot.model.perf.met(met=met,perf.met='Test Accuracy',ttl='',colors = dz.colors[dz.short %in% names(met)])
p.roc.d <- plot.model.roc(met=dz.res,ttl='',colors=dz.colors[dz.short %in% names(dz.res)])


# cluster label
load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
clusterColors <- getClusterColors(k)

load(file = paste(savedir,'predictcluster_RFperf',extralab,'.RData',sep=''))

met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Sensitivity))
p.sen.c <- plot.model.perf.met(met=met,perf.met='Test Sensitivity',ttl='',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Specificity))
p.spec.c <- plot.model.perf.met(met=met,perf.met='Test Specificity',ttl='',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$AUC))
p.auc.c <- plot.model.perf.met(met=met,perf.met='Test AUC',ttl='',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$OverallAccuracy))
p.acc.c <- plot.model.perf.met(met=met,perf.met='Test Accuracy',ttl='',colors = clusterColors)
p.roc.c <- plot.model.roc(met=cluster.res,ttl='',colors=clusterColors)

# align all plots
p.all <- plot_grid(plotlist = 
	list(p.sen.d,remove.y.ticklabels(p.spec.d),remove.y.ticklabels(p.acc.d),remove.y.ticklabels(p.auc.d),p.roc.d,
		p.sen.c,remove.y.ticklabels(p.spec.c),remove.y.ticklabels(p.acc.c),remove.y.ticklabels(p.auc.c),p.roc.c),
		align = 'hv',nrow=2,axis = 'b')


ggsave(filename = paste(savedir,'RFPerformanceClustersDisease',extralab,'.pdf',sep=''),plot = p.all,
       height = 10,width=18,units='cm')