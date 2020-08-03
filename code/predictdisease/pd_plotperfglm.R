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

load(file = paste(savedir,'predictdz_GLMperf',extralab,'.RData',sep=''))

if(grepl('CSFOnly',extralab) | grepl('CSFGene',extralab)){
	pal <- colorRampPalette(brewer.pal(name = 'Blues',n=9))
	dz.colors <- pal(9)[3:(3+length(dz.res)-1)] # pick good range of blues
}
if(grepl('GeneOnly',extralab)){
	pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
	dz.colors <- pal(length(dz.short)) # assign a color to each disease
	# select only diseases included in analysis
	dz.colors <- dz.colors[dz.short %in% names(dz.res)]
}

met <- lapply(dz.res, function(a) sapply(a, function(c) c$Sensitivity))
p.sen.d <- plot.model.perf.met(met=met,perf.met='Test Sensitivity',ttl='',colors = dz.colors)
met <- lapply(dz.res, function(a) sapply(a, function(c) c$Specificity))
p.spec.d <- plot.model.perf.met(met=met,perf.met='Test Specificity',ttl='',colors = dz.colors)
met <- lapply(dz.res, function(a) sapply(a, function(c) c$AUC))
p.auc.d <- plot.model.perf.met(met=met,perf.met='Test AUC',ttl='',colors = dz.colors)
p.roc.d <- plot.model.roc(met=dz.res,ttl='',colors=dz.colors)
p.fw.d <- plot.featureweights.lm.contintsplit(dz.res,df,dz.colors)

# cluster label
load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
clusterColors <- getClusterColors(k)

load(file = paste(savedir,'predictcluster_GLMperf',extralab,'.RData',sep=''))

met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Sensitivity))
p.sen.c <- plot.model.perf.met(met=met,perf.met='Test Sensitivity',ttl='',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$Specificity))
p.spec.c <- plot.model.perf.met(met=met,perf.met='Test Specificity',ttl='',colors = clusterColors)
met <- lapply(cluster.res, function(a) sapply(a, function(c) c$AUC))
p.auc.c <- plot.model.perf.met(met=met,perf.met='Test AUC',ttl='',colors = clusterColors)
p.roc.c <- plot.model.roc(met=cluster.res,ttl='',colors=clusterColors)
p.fw.c <- plot.featureweights.lm.contintsplit(cluster.res,df,clusterColors)

# align all plots
p.d.list <- list(p.sen.d,remove.y.ticklabels(p.spec.d),remove.y.ticklabels(p.auc.d),p.roc.d,p.fw.d)
p.d <- plot_grid(plotlist = p.d.list,align = 'h',nrow=1,axis = 'b',rel_widths = c(1.2,1,1,1,1.2))
p.c.list <- list(p.sen.c,remove.y.ticklabels(p.spec.c),remove.y.ticklabels(p.auc.c),p.roc.c,p.fw.c)
p.c <- plot_grid(plotlist = p.c.list,align = 'h',nrow=1,axis = 'b',rel_widths = c(1.2,1,1,1,1.2),
                 rel_heights = c(rep(4,4),0.1))
rel_heights <- c(1.2,1) #rel_heights <- c(1,1.05)
if(grepl('GeneOnly',extralab)){rel_heights <- c(1.2,1)}
p.all <- plot_grid(plotlist= list(p.d,p.c), align = 'hv',nrow = 2,axis='b',
		rel_heights = rel_heights)

w.multiplier <- 1; h.multiplier <- 1
if(grepl('CSFGene',extralab) | grepl('GeneOnly',extralab)){w.multiplier <- 1.22; h.multiplier <- 0.75}
ggsave(filename = paste(savedir,'GLMPerformanceClustersDisease',extralab,'.pdf',sep=''),plot = p.all,
       height = 14*h.multiplier,width=19*w.multiplier,units='cm')

if(extralab == 'CSFOnlyAddNormalMMSE'){
  Fig6ad <- lapply(c(p.d.list,p.c.list),function(X) X$data)
  save(Fig6ad,file = paste(params$sourcedata.dir,'Fig6a-d_',extralab,'SourceData.RData',sep=''))
}
if(extralab == 'CSFGeneAddNormalMMSE'){
  FigS16ad <- lapply(c(p.d.list,p.c.list),function(X) X$data)
	save(FigS16ad,file = paste(params$sourcedata.dir,'FigS16a-d_',extralab,'SourceData.RData',sep=''))
}
if(extralab == 'GeneOnlyAddNormalMMSE'){
	FigS18ad <- lapply(c(p.d.list,p.c.list),function(X) X$data)
  save(FigS18ad,file = paste(params$sourcedata.dir,'FigS18ad_',extralab,'SourceData.RData',sep=''))
}