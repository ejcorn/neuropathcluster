rm(list = setdiff(ls(), c("params","dz.exc","n.dx","exc.cl")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'genecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
if(!is_empty(DisconnectedSubjects)){
  microSample <- microSample[-DisconnectedSubjects,]
  patientSample <- patientSample[-DisconnectedSubjects,]
}

clusterColors <- getClusterColors(k)

load(file=paste(savedir,'AlleleTablesCluster.RData',sep=''))

#################################################
### Plot allele proportions excluding disease ###
#################################################

Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(Genes$INDDID),],dz.exc,n.dx)
Rem.Mask <- Rem.Mask | partitionSample != paste('Cluster',exc.cl)

for(g.i in names(Allele.Tables)){
  Allele.Tables[[g.i]] <- Allele.Tables[[g.i]][Rem.Mask,]
}
partitionSample <- partitionSample[Rem.Mask]

# display number of patients left in each cluster
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Genes: Exclude',dz.exc))
print(cluster.counts)

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample,n=20)
# apply mask to relevant variables
list[clusterColors,clusterNames] <- 
  lapply(list(clusterColors,clusterNames), function(X) X.n.exclude(X,cluster.counts))

partitionSample <- flexdim.rowmask(partitionSample,Rem.Cl)

G.color.inds <- list(APOE=c(1,3,7),MAPTHaplotype=c(4:5))
AlleleProportion.byCluster <- list()
for(g.i in names(Allele.Tables)){
  AlleleProportion.byCluster[[g.i]] <- sapply(clusterNames, function(k.i) colSums(Allele.Tables[[g.i]][partitionSample == k.i,]) / (2*sum(partitionSample == k.i)))
  Allele.Proportions <- AlleleProportion.byCluster[[g.i]]
  df <- data.frame(y=as.vector(Allele.Proportions),
                   g=as.vector(sapply(clusterNames, function(i) matrix(rownames(Allele.Proportions),ncol=1))),
                   x=rep(clusterNames,each=nrow(Allele.Proportions)))
  save(df,file=paste(savedir,'FigS5a-b_',g.i,'SourceData.RData',sep=''))
  pal.g <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
  pal.g <- pal.g(12)[G.color.inds[[g.i]]]
  p <- ggplot(data=df,aes(y=y,x=x,fill=g)) + geom_col(position=position_dodge(width = 0.9)) + theme_classic()+ 
    scale_fill_manual(values=pal.g,name='') + ylab('Proportion of Cluster') + xlab('') + ggtitle(g.i) +
    theme(
      #legend.position = c(0.32,0.8),
      legend.key.size = unit(0.1,'in'),
      plot.title = element_text(face='bold',size=8,hjust=0.5)) +
    theme(text= element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color=clusterColors))
  p
  
  ggsave(filename = paste(savedir,g.i,'AlleleProportionsByClusterExclude',dz.exc,'.pdf',sep=''),plot = p,
         height = 2,width=3,units='in')
}