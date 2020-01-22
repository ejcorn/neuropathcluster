rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

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

cog <- read.csv(paste(params$opdir,'processed/MoCA_processed.csv',sep=''))
cog <- cog[cog$INDDID %in% INDDIDs,]
meanMOCA <- sapply(unique(cog$INDDID), function(ID) mean(cog$MoCATotal[cog$INDDID == ID],na.rm= T))

partitionSample <- as.character(partition[INDDIDs %in% unique(cog$INDDID)])
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))
print(paste('Cog + Path n =',length(partitionSample)))

# plot mean moca score
clusterColors <- getClusterColors(k)
df <- data.frame(x=partitionSample,y=meanMOCA)
comps <- lapply(as.data.frame(combn(levels(df$x),m = 2)),function(x) as.character(x))
p <- ggplot(data=df,aes(x=x,y=y,fill=x)) + geom_boxplot() + theme_classic() +
  ylab('Mean MOCA') + xlab('') + scale_fill_manual(values=clusterColors) +
  theme(legend.position = 'none',text=element_text(size=8),
        axis.text.x = element_text(color=clusterColors)) +
  stat_compare_means(data=df,size=2,mapping=aes(label = ..p.adj..),
                     comparisons=comps,formula = y~x,p.adjust.method='holm')

p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = file.path(savedir,'MeanMOCAbyClusterLouvain.pdf'),height = unit(3,'in'),width=unit(3,'in'))
plot(ggplot_gtable(p))
dev.off()

# look by MOCA subscores
COI <- grep('Total',colnames(cog))
cog.moca <- lapply(unique(cog$INDDID), function(ID) 
  colMeans(cog[cog$INDDID == ID,COI],na.rm= T)) # average over multiple tests
cog.moca <- do.call(what='rbind',args=cog.moca)
cog.moca <- as.data.frame(cog.moca[,colSums(is.na(cog.moca)) == 0]) # remove the qualitative features
exclude <- c('MoCATotal')
cog.moca <- cog.moca[,!(names(cog.moca) %in% exclude)]

MOCA.mean.sub <- sapply(cog.moca, function(X) sapply(1:k, function(i) mean(X[partitionSample==paste('Cluster',i)])))
cls <- sapply(1:k,function(i) paste('Cluster',i))
df <- data.frame(g=rep(cls,ncol(MOCA.mean.sub)),y=as.vector(MOCA.mean.sub),x=rep(colnames(MOCA.mean.sub),each=k))
p <- ggplot(data=df,aes(x=x,y=y,fill=g)) + geom_col(position = position_dodge(0.9)) + theme_classic()+
  scale_fill_manual(values=clusterColors,name='') + ylab('Score') + xlab('') +
  theme(legend.position = c(0.1,0.75),legend.key.size = unit(0.1,'in')) +
  theme(text= element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color=clusterColors))
p

ggsave(filename = file.path(savedir,'MOCASubscoreMeansbyClusterLouvain.pdf'),plot = p,
       height = 3,width=5,units='in')

df <- data.frame(x=rep(partitionSample,ncol(cog.moca)),y=as.vector(as.matrix(cog.moca)),g=rep(colnames(cog.moca),each=nrow(cog.moca)))
save(df,file = paste(savedir,'Fig3a_SourceData.RData',sep=''))

p <- ggplot(data=df,aes(x=x,y=y,fill=x)) + geom_boxplot(outlier.size=0.5) + theme_classic() +
  facet_wrap(~g) +
  scale_fill_manual(values=clusterColors,name='') + ylab('Score') +
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + # leg size c(0.1,0.75)
  theme(text= element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color=clusterColors)) +
  stat_compare_means(size=2,comparisons = comps) + xlab('')
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = file.path(savedir,'MOCASubscoreBoxplotsbyClusterLouvain.pdf'),height = unit(4,'in'),width=unit(3.5,'in'),useDingbats = F)
plot(ggplot_gtable(p))
dev.off()


# to get df and effect size for each test:

clusterNames <- sapply(1:k, function(k.i) paste('Cluster',k.i))
results <- samp.size <- list()
for(MoCA.Category in colnames(cog.moca)){
  results[[MoCA.Category]] <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
  samp.size[[MoCA.Category]] <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
  for(k1 in clusterNames){
    for(k2 in clusterNames){
      m <- wilcox.test(cog.moca[partitionSample == k1,MoCA.Category],cog.moca[partitionSample == k2,MoCA.Category],conf.int = TRUE)
      results[[MoCA.Category]][k1,k2] <- m$estimate
      samp.size[[MoCA.Category]][k1,k2] <- sum(partitionSample %in% c(k1,k2))
    }
  }
}