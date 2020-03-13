rm(list = setdiff(ls(), c("params","extralab",'dz.exc','n.dx','exc.cl')))
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

# cluster label
load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
clusterColors <- getClusterColors(k)

load(file = paste(savedir,'predictcluster_GLMperf',extralab,'.RData',sep=''))
load(file = paste(savedir,'predictdz_GLMperf',extralab,'.RData',sep=''))
pal <- colorRampPalette(brewer.pal(name = 'Blues',n=9))
dz.colors <- pal(9)[3:(3+length(dz.res)-1)] # pick good range of blues

setups <- list(Clusters=list(colors=clusterColors,res=cluster.res),
     Diseases=list(colors=dz.colors,res=dz.res))
for(setup.name in names(setups)){
  setup <- setups[[setup.name]] # run through set ups
  colors.setup <- setup$colors # extract colors
  res <- setup$res  # extract data
  
  met <- lapply(res, function(a) sapply(a, function(c) c$Sensitivity))
  p.sen <- plot.model.perf.met(met=met,perf.met='Test Sensitivity',ttl='',colors = colors.setup)
  met <- lapply(res, function(a) sapply(a, function(c) c$Specificity))
  p.spec <- plot.model.perf.met(met=met,perf.met='Test Specificity',ttl='',colors = colors.setup)
  met <- lapply(res, function(a) sapply(a, function(c) c$AUC))
  p.auc <- plot.model.perf.met(met=met,perf.met='Test AUC',ttl='',colors = colors.setup)
  p.roc <- plot.model.roc(met=res,ttl='',colors=colors.setup)
  
  p.list <- list(p.sen,remove.y.ticklabels(p.spec),remove.y.ticklabels(p.auc),p.roc)
  p.all <- plot_grid(plotlist = p.list, align = 'h',nrow=1,axis = 'b',rel_widths = c(1.2,1,1,1))
  
  if(extralab == "CSFOnlyAddNormalMMSEExcludeAlzheimer's disease"){
    if(setup.name == 'Clusters'){
      FigS17ab <- lapply(p.list,function(X) X$data)
  	  save(FigS17ab,file = paste(params$sourcedata.dir,'FigS17a-b_',extralab,'SourceData.RData',sep=''))
    }
  }
  
  ggsave(filename = paste(savedir,'GLMPerformance',setup.name,extralab,'FromC',paste0(exc.cl,collapse=','),'NDx',n.dx,'.pdf',sep=''),plot = p.all,
         height = 5,width=18,units='cm')
}


load(file = paste(savedir,'predictcluster_RFperf',extralab,'.RData',sep=''))
load(file = paste(savedir,'predictdz_RFperf',extralab,'.RData',sep=''))
setups <- list(Clusters=list(colors=clusterColors,res=cluster.res),
               Diseases=list(colors=dz.colors,res=dz.res))

for(setup.name in names(setups)){
  setup <- setups[[setup.name]] # run through set ups
  colors.setup <- setup$colors # extract colors
  res <- setup$res  # extract data
  
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
  
  if(extralab == "CSFGeneAddNormalMMSEExcludeAlzheimer's disease"){
    if(setup.name == 'Clusters'){
      FigS17cd <- lapply(p.list,function(X) X$data)
      save(FigS17cd,file = paste(params$sourcedata.dir,'FigS17c-d_',extralab,'SourceData.RData',sep=''))
    }
  }
  
  ggsave(filename = paste(savedir,'RFPerformance',setup.name,extralab,'FromC',paste0(exc.cl,collapse=','),'NDx',n.dx,'.pdf',sep=''),plot = p.all,
         height = 5,width=18,units='cm')
}