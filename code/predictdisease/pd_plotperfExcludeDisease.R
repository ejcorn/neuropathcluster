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

# cluster label
load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
clusterColors <- getClusterColors(k)

load(file = paste(savedir,'predictcluster_GLMperf',extralab,'.RData',sep=''))

met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Sensitivity))
p.sen <- plot.model.perf.met(met=met,perf.met='Test Sensitivity',ttl='',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Specificity))
p.spec <- plot.model.perf.met(met=met,perf.met='Test Specificity',ttl='',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$AUC))
p.auc <- plot.model.perf.met(met=met,perf.met='Test AUC',ttl='',colors = clusterColors)
p.roc <- plot.model.roc(met=cluster.res,ttl='',colors=clusterColors)

p.all <- plot_grid(plotlist = 
	list(p.sen,remove.y.ticklabels(p.spec),remove.y.ticklabels(p.auc),p.roc),
	align = 'h',nrow=1,axis = 'b',
		rel_widths = c(1.2,1,1,1))

if(extralab == "CSFOnlyExcludeAlzheimer's disease"){
	save(cluster.res,file = paste(savedir,'FigS11a-b_',extralab,'SourceData.RData',sep=''))
}

ggsave(filename = paste(savedir,'GLMPerformanceClusters',extralab,'.pdf',sep=''),plot = p.all,
       height = 5,width=18,units='cm')
