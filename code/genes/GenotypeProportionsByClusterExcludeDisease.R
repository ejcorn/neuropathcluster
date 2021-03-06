rm(list = setdiff(ls(), c("params",'exc.cl','dz.exc','n.dx')))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'genecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

clusterColors <- getClusterColors(k)

GOIs <- c('APOE','MAPTHaplotype')#,'C9orf72','LRRK2')
Genes <- read.csv(paste(params$opdir,'processed/Genetics_processed.csv',sep=''),stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% INDDIDs,]
missing.mask <- rowSums(Genes[,GOIs] != '') == length(GOIs)
Genes <- Genes[missing.mask,]  # remove missing
Genes.df <- Genes[,GOIs] # select genes of interest

clusterColors <- getClusterColors(k)

load(file=paste(savedir,'AlleleTablesCluster.RData',sep=''))

############################################
### Plot Genotype proportions by cluster ###
############################################

Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(Genes$INDDID),],dz.exc,n.dx)
Rem.Mask <- Rem.Mask | !partitionSample %in% paste('Cluster',exc.cl)
Genes.df <- Genes.df[Rem.Mask,]
partitionSample <- partitionSample[Rem.Mask]

# display number of patients left in each cluster
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Genes: Exclude',dz.exc))
print(cluster.counts)

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample,n=9)
# apply mask to relevant variables
list[clusterColors,clusterNames] <- 
  lapply(list(clusterColors,clusterNames), function(X) X.n.exclude(X,cluster.counts,n=9))

partitionSample <- flexdim.rowmask(partitionSample,Rem.Cl)

GenotypeProportion.byCluster <- list()
p.names <- list('Blues','Reds')
names(p.names) <- names(Genes.df)
for(g.i in names(Genes.df)){
  
  GenotypeProportion.byCluster[[g.i]] <- sapply(clusterNames, function(k.i) count.ejc(Genes.df[partitionSample == k.i,g.i],items = unique(Genes.df[,g.i])) / sum(partitionSample==k.i))
  Genotype.Proportions <- GenotypeProportion.byCluster[[g.i]]
  df <- data.frame(y=as.vector(Genotype.Proportions),
                   g=as.vector(sapply(clusterNames, function(i) matrix(rownames(Genotype.Proportions),ncol=1))),
                   x=rep(clusterNames,each=nrow(Genotype.Proportions)))
  save(df,file=paste(params$sourcedata.dir,'FigS14a-b_SourceData',g.i,'Exclude',dz.exc,'From',paste(exc.cl,collapse=','),'.RData',sep=''))
  
  p <- ggplot(data=df,aes(y=y,x=x,fill=g)) + geom_col(position=position_dodge(width = 0.9)) + theme_classic()+ 
    scale_fill_brewer(palette=p.names[[g.i]],name='') + scale_y_continuous(limits=c(0,1)) +
    ylab('Proportion of Cluster') + xlab('') + ggtitle(g.i) +
    theme(
      #legend.position = c(0.32,0.8),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.key.size = unit(0.1,'in'),
      plot.title = element_text(face='bold',size=8,hjust=0.5)) +
    theme(text= element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color=clusterColors))
  p
  
  ggsave(filename = paste(savedir,g.i,'GenotypeProportionsByClusterExclude',dz.exc,'FromC',paste0(exc.cl,collapse=','),'NDX',n.dx,'.pdf',sep=''),plot = p,
         height = 4.5,width=9,units='cm')
}