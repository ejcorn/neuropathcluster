rm(list = setdiff(ls(), c("params","dz.exc","n.dx","exc.cl")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'csfcluster/',sep='')
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

###########
### CSF ###
###########

CSF.name <- 'Luminex'
CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% INDDIDs,]
CSF.vars <- c('LuminexTTau','LuminexPTau','LuminexAbeta42')

CSF.by.pt <- lapply(sort(unique(CSF$INDDID)), function(id) # order tests by date w/in subjs
  CSF[CSF$INDDID==id,c('CSFDate',CSF.vars)])
CSF.by.pt <- lapply(CSF.by.pt, function(X) X[order(X$CSFDate),-1])
#CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) colMeans(X)))
CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) X[1,])) # just use first sample
# for some patients, could compute a feature based on change in CSF proteins
#CSF.diff <- do.call('rbind',lapply(CSF.by.pt, function(X) X[nrow(X),] - X[1,]))
CSF.mean <- cbind(data.frame(INDDID=sort(unique(CSF$INDDID))),CSF.mean)
#CSF.mean <- CSF.mean[-(which(CSF.mean$LuminexTTau >1000)),] # ttau outliers

partitionSample <- as.character(partition[INDDIDs %in% unique(CSF.mean$INDDID)])
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))

# Remove all pts of certain diseases to see what is driving cluster diffs
Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(CSF.mean$INDDID),],dz.exc,n.dx)
# Find AD within cluster 2
Rem.Mask <- Rem.Mask & partitionSample == paste('Cluster',exc.cl)
# Give them a separate label
partitionSample[Rem.Mask] <- 'tmp'
Rem.Mask <- partitionSample == 'tmp' | partitionSample == paste('Cluster',exc.cl)
# Isolate all cluster 2 people
partitionSample <- partitionSample[Rem.Mask]
# Relabel them as a cluster for temporary implementation ease
partitionSample[partitionSample == 'tmp'] <- paste('Cluster',(1:k)[-exc.cl][1])

CSF.mean <- CSF.mean[Rem.Mask,]
CSF.mean <- CSF.mean[,-1]

# display number of patients left in each cluster
print(paste('CSF: Exclude',dz.exc))
cluster.counts <- cluster.count(partitionSample,k)

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample)
# apply mask to relevant variables
clusterColors <- X.n.exclude(clusterColors,cluster.counts)
list[partitionSample,CSF.mean] <- 
  lapply(list(partitionSample,CSF.mean), function(X) flexdim.rowmask(X,Rem.Cl))

# display counts after excluding small clusters
cluster.counts <- cluster.count(partitionSample,k)
print(paste('CSF: Exclude',dz.exc))
print(cluster.counts)
print(paste('Exclude',dz.exc,'n =',length(partitionSample)))

ymax <- 1000
n.breaks <- 4 # set the number of breaks on the original plot so the tick label scaling is done automatically
df <- data.frame(y=as.vector(as.matrix(CSF.mean)),
                 x=rep(colnames(CSF.mean),each=nrow(CSF.mean)),
                 g=rep(partitionSample,ncol(CSF.mean)))
comps <- lapply(as.data.frame(combn(levels(df$g),m = 2)),function(x) as.character(x))
p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5) + theme_classic()+
  facet_wrap(~x) + scale_y_continuous(breaks=pretty_breaks(n=n.breaks)(df$y)) +
  scale_fill_manual(values=clusterColors,name='') + ylab('CSF Protein (pg/ml)') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=2,comparisons = comps,position=5) 
#p <- ggplot_build(p)
p <- mult.comp.ggpubr(p)
p <- fix.yscale(p,ymax,n.breaks)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/CSF',CSF.name,'BoxplotsbyClusterLouvainCluster2ADvsNonAD',dz.exc,'NDx',n.dx,'.pdf',sep=''),height = unit(4,'in'),width=unit(3.5,'in'),useDingbats = F)
plot(ggplot_gtable(p))
dev.off()
