rm(list = setdiff(ls(), c("params",'dz.exc','exc.cl','n.dx')))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1)] # Get rid of index column only

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
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

cog <- read.csv(paste(params$opdir,'processed/MoCA_processed.csv',sep=''),stringsAsFactors = F)[,-1] # remove index column
cog <- cog[cog$INDDID %in% INDDIDs,]
cog$TestDate <- as.Date(cog$TestDate,format='%m/%d/%Y')
cog.by.subj <- lapply(unique(cog$INDDID), function(ID) 
  cog[cog$INDDID == ID,]) # extract all tests for patient
cog.by.subj <- lapply(cog.by.subj, function(X) X[order(X$TestDate),]) # sort so more recent tests come last
cog.by.subj <- lapply(cog.by.subj, function(X) X[nrow(X),]) # only keep most recent test (closest to autopsy)
cog.by.subj <- do.call('rbind',cog.by.subj)

names(partition) <- as.character(INDDIDs) # name partition by INDDIDs
partitionSample <- partition[as.character(cog.by.subj$INDDID)] # index by INDDIDs from cog dataframe so it matches
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))
clusterColors <- getClusterColors(k)
print(paste('Cog + Path n =',length(partitionSample)))

# Select only pts of certain diseases to see how cluster parse heterogeneity within a disease
Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(cog$INDDID),],dz.exc,n.dx)
# Only remove Braak-CERAD AD from Cluster 2
Rem.Mask <- Rem.Mask | !partitionSample %in% paste('Cluster',exc.cl)

partitionSample <- partitionSample[Rem.Mask]
cog.by.subj <- cog.by.subj[Rem.Mask,]

# display number of patients left in each cluster
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Cog: Exclude',dz.exc))
print(cluster.counts)
print(paste('Exclude',dz.exc,'n =',length(partitionSample)))

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample)
# apply mask to relevant variables
clusterColors <- X.n.exclude(clusterColors,cluster.counts)
list[partitionSample,cog.by.subj] <- 
  lapply(list(partitionSample,cog.by.subj), function(X) flexdim.rowmask(X,Rem.Cl))

# display counts after excluding small clusters
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Cog: Exclude',dz.exc))
print(cluster.counts)
print(paste('Cog + Path n =',length(partitionSample)))
# plot mean moca score
meanMOCA <- cog.by.subj$MoCATotal
df <- data.frame(x=partitionSample,y=meanMOCA)
comps <- lapply(as.data.frame(combn(levels(df$x),m = 2)),function(x) as.character(x))
p <- ggplot(data=df,aes(x=x,y=y,fill=x)) + geom_boxplot() + theme_classic() +
  ylab('Mean MOCA') + xlab('') + scale_fill_manual(values=clusterColors) +
  scale_y_continuous(breaks = c(0,10,20,30))+
  theme(legend.position = 'none',text=element_text(size=8),
        axis.text.x = element_text(color=clusterColors,angle=90,vjust=0.5)) + theme(axis.text.y = element_text(size=6))+
  stat_compare_means(size=2,comparisons=comps)

p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste0(savedir,'/OverallMoCABoxplotsbyClusterLouvainExclude',dz.exc,'FromC',paste(exc.cl,collapse = ','),'NDx',n.dx,'.pdf'),height = 4.5/2.54,width=4.5/2.54,useDingbats = F)
plot(ggplot_gtable(p))
dev.off()

# look by MOCA subscores
COI <- grep('Total',colnames(cog))
cog.moca <- cog.by.subj[,COI] # isolate scores with Total in name from most recent test for each patient
cog.moca <- as.data.frame(cog.moca[,colSums(is.na(cog.moca)) == 0]) # remove the qualitative features
exclude <- c('MoCATotal')
cog.moca <- cog.moca[,!(names(cog.moca) %in% exclude)]
cog.moca <- cog.moca[,order(colnames(cog.moca))] # reorder column names to alphabetical
colnames(cog.moca) <- c('Digit Attention','Repetition','Naming','Orientation','Recall','Visuospatial')

df <- data.frame(x=rep(partitionSample,ncol(cog.moca)),y=as.vector(as.matrix(cog.moca)),g=rep(colnames(cog.moca),each=nrow(cog.moca)))
save(df,file = paste(savedir,'Fig3a_SourceData.RData',sep=''))

p <- ggplot(data=df,aes(x=x,y=y,fill=x)) + geom_boxplot(outlier.size=0.5) + theme_classic() +
  facet_wrap(~g) + scale_y_continuous(breaks=c(0:max(cog.moca))) + # make ticks only for 1:max moca subscore score, which is 6
  scale_fill_manual(values=clusterColors,name='') + ylab('Score') +
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + # leg size c(0.1,0.75)
  theme(text= element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color=clusterColors)) +
  stat_compare_means(size=2,comparisons = comps) + xlab('')
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste0(savedir,'/MOCASubscoreBoxplotsbyClusterLouvainExclude',dz.exc,'FromC',paste(exc.cl,collapse = ','),'NDx',n.dx,'.pdf'),height = unit(4,'in'),width=unit(3.5,'in'),useDingbats = F)
plot(ggplot_gtable(p))
dev.off()

# to get df and effect size for each test:


clusterNames <- sort(unique(partitionSample))
k <- length(clusterNames)
results <- latex <- samp.size <- pvals <- list()
for(MoCA.Category in colnames(cog.moca)){
  results[[MoCA.Category]] <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
  samp.size[[MoCA.Category]] <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
  pvals[[MoCA.Category]] <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
  latex[[MoCA.Category]] <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
  for(k1 in clusterNames){
    for(k2 in clusterNames){
      if(k1 != k2){m <- wilcox.test(cog.moca[partitionSample == k1,MoCA.Category],cog.moca[partitionSample == k2,MoCA.Category],conf.int = TRUE)}
      else if(k1==k2){m <- list(estimate=NA,p.value=NA)}
      
      results[[MoCA.Category]][k1,k2] <- m$estimate
      samp.size[[MoCA.Category]][k1,k2] <- sum(partitionSample %in% c(k1,k2))
      pvals[[MoCA.Category]][k1,k2] <- m$p.value
    }
  }
}

pvals <- list.fdr.correct(pvals) # fdr correct

for(MoCA.Category in colnames(cog.moca)){
  for(k1 in clusterNames){
    for(k2 in clusterNames){
      latex[[MoCA.Category]][k1,k2] <- paste0('$M_{',substr(k1,nchar(k1),nchar(k1)),'-',substr(k2,nchar(k2),nchar(k2)),'}=',
                                              round(results[[MoCA.Category]][k1,k2],2),'$, $n=',samp.size[[MoCA.Category]][k1,k2],'$, $p_','\\','mathrm{FDR}=',signif(pvals[[MoCA.Category]][k1,k2],2),'$')
    }
  }
}
