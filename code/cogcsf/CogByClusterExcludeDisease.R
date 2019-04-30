rm(list = setdiff(ls(), c("params","dz.exc","n.dx","exc.cl")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

if(!is_empty(DisconnectedSubjects)){
  microSample <- microSample[-DisconnectedSubjects,]
  patientSample <- patientSample[-DisconnectedSubjects,]
}
if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)

###########
### Cog ###
###########

cog <- read.csv(paste(params$opdir,'processed/MoCA_processed.csv',sep=''))
cog <- cog[cog$INDDID %in% INDDIDs,]

# look by MOCA subscores
COI <- grep('Total',colnames(cog))
cog.moca <- lapply(unique(cog$INDDID), function(ID) 
  colMeans(cog[cog$INDDID == ID,COI],na.rm= T)) # average over multiple tests
cog.moca <- do.call(what='rbind',args=cog.moca)
cog.moca <- as.data.frame(cog.moca[,colSums(is.na(cog.moca)) == 0]) # remove the qualitative features
exclude <- c('MoCATotal')
cog.moca <- cog.moca[,!(names(cog.moca) %in% exclude)]

partitionSample <- as.character(partition[INDDIDs %in% unique(cog$INDDID)])
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))

clusterColors <- getClusterColors(k)

# Select only pts of certain diseases to see how cluster parse heterogeneity within a disease
Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(cog$INDDID),],dz.exc,n.dx)
# Only remove Braak-CERAD AD from Cluster 2
Rem.Mask <- Rem.Mask | partitionSample != paste('Cluster',exc.cl)

partitionSample <- partitionSample[Rem.Mask]
cog.moca <- cog.moca[Rem.Mask,]

# display number of patients left in each cluster
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Cog: Exclude',dz.exc))
print(cluster.counts)
print(paste('Exclude',dz.exc,'n =',length(partitionSample)))

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample)
# apply mask to relevant variables
clusterColors <- X.n.exclude(clusterColors,cluster.counts)
list[partitionSample,cog.moca] <- 
	lapply(list(partitionSample,cog.moca), function(X) flexdim.rowmask(X,Rem.Cl))

# display counts after excluding small clusters
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Cog: Exclude',dz.exc))
print(cluster.counts)

df <- data.frame(x=rep(partitionSample,ncol(cog.moca)),y=as.vector(as.matrix(cog.moca)),g=rep(colnames(cog.moca),each=nrow(cog.moca)))
comps <- lapply(as.data.frame(combn(levels(df$x),m = 2)),function(x) as.character(x))
p <- ggplot(data=df,aes(x=x,y=y,fill=x)) + geom_boxplot(outlier.size=0.5) + theme_classic() +
  facet_wrap(~g) +
  scale_fill_manual(values=clusterColors,name='') + ylab('Score') +
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + # leg size c(0.1,0.75)
  theme(text= element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color=clusterColors)) +
  stat_compare_means(size=2,comparisons = comps) + xlab('')
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/MOCASubscoreBoxplotsbyClusterLouvainExclude',dz.exc,'NDx',n.dx,'.pdf',sep=''),height = unit(4,'in'),width=unit(3.5,'in'),useDingbats = F)
plot(ggplot_gtable(p))
dev.off()
