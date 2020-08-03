rm(list = setdiff(ls(), "params"))
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
meanMOCA <- cog.by.subj$MoCATotal

names(partition) <- as.character(INDDIDs) # name partition by INDDIDs
partitionSample <- partition[as.character(cog.by.subj$INDDID)] # index by INDDIDs from cog dataframe so it matches
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))
print(paste('Cog + Path n =',length(partitionSample)))

patientSample$AutopsyDate <- as.Date(patientSample$AutopsyDate,format='%Y-%m-%d')
df.all <- merge(cog.by.subj,patientSample,by='INDDID')
data.date <- '122219'
demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics',data.date,'.csv',sep=''),stringsAsFactors = F)
df.all <- merge(df.all,demo,by='INDDID') # merge with demographics to control for age at death
df.all$AutopsyMinusMoCA <- df.all$AgeatTest - df.all$AgeatDeath  
df.all$Cluster <- partitionSample[as.character(df.all$INDDID)]

# p <- ggplot(df.all) + geom_boxplot(aes(x=Cluster,y=AutopsyMinusMoCA/ 365,fill=Cluster)) +
#   scale_fill_manual(values = getClusterColors(k))+ xlab('') + ylab('Last MoCA - Autopsy Date (years)')+ theme_classic()+
#   theme(text=element_text(size=8),axis.text.x=element_text(angle=90),legend.position = 'None')
# p

clusterColors <- getClusterColors(k)
df <- data.frame(x=df.all$Cluster,y=df.all$AutopsyMinusMoCA) # /365 to convert days to years
comps <- lapply(as.data.frame(combn(levels(df$x),m = 2)),function(x) as.character(x))
p <- ggplot(data=df,aes(x=x,y=y,fill=x)) + geom_boxplot() + theme_classic() + 
  ylab('Last MoCA - Death (years)') + xlab('') + scale_fill_manual(values=clusterColors) +
  theme(legend.position = 'none',text=element_text(size=8),
        axis.text.x = element_text(color=clusterColors,angle=90,vjust=0.5)) + theme(axis.text.y = element_text(size=6))+
  stat_compare_means(size=2,comparisons=comps)
p
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = file.path(savedir,'LastMoCAvsDeathDatebyClusterLouvain.pdf'),height = 4.5/2.54,width=4.5/2.54,useDingbats = F)
plot(ggplot_gtable(p))
dev.off()
save(df.all, file = paste0(params$sourcedata.dir,'FigS4b_SourceData_MoCAvsDateofDeath.RData'))

# plot mean moca score
clusterColors <- getClusterColors(k)
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
pdf(file = file.path(savedir,'MeanMOCAbyClusterLouvain.pdf'),height = 4.5/2.54,width=4.5/2.54,useDingbats = F)
plot(ggplot_gtable(p))
dev.off()
save(df, file = paste0(params$sourcedata.dir,'FigS4a_SourceData_MeanMoCA.RData'))

# look by MOCA subscores
COI <- grep('Total',colnames(cog))
cog.moca <- cog.by.subj[,COI] # isolate scores with Total in name from most recent test for each patient
cog.moca <- as.data.frame(cog.moca[,colSums(is.na(cog.moca)) == 0]) # remove the qualitative features
exclude <- c('MoCATotal')
cog.moca <- cog.moca[,!(names(cog.moca) %in% exclude)]
cog.moca <- cog.moca[,order(colnames(cog.moca))] # reorder column names to alphabetical
colnames(cog.moca) <- c('Digit Attention','Repetition','Naming','Orientation','Recall','Visuospatial')
  
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
save(df,file = paste(params$sourcedata.dir,'Fig4a_SourceData.RData',sep=''))

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

get.mat.inds <- function(mat,inds) sapply(1:ncol(inds), function(j) mat[inds[1,j],inds[2,j]]) # iterate through unique combinations of clusters and get p-values from p-value matrix
cluster.combs <- combn(clusterNames,2)
p.table <- as.data.frame(sapply(pvals, function(X) get.mat.inds(X,cluster.combs)))
rownames(p.table) <- sapply(1:ncol(cluster.combs), function(j) paste0(cluster.combs[1,j],' vs. ',cluster.combs[2,j]))
p.table.print <- signif(p.table,digits = 2)
p.table.print[p.table.print<0.001] <- 'p < 0.001'
xtable(p.table.print,caption = 'FDR-corrected $p$-values for Figure \\ref{fig:figure4}',label = 'table:figure6pvals')

which.sig.any <- which(rowSums(p.table<0.05)>1)
comps <- lapply(as.data.frame(cluster.combs[,which.sig.any]),function(x) as.character(x))

p <- ggplot(data=df,aes(x=x,y=y,fill=x)) + geom_boxplot(outlier.size=0.5) + theme_classic() +
  facet_wrap(~g) + scale_y_continuous(breaks=c(0:max(cog.moca))) + # make ticks only for 1:max moca subscore score, which is 6
  scale_fill_manual(values=clusterColors,name='') + ylab('Score') +
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + # leg size c(0.1,0.75)
  theme(text= element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color=clusterColors)) +
  stat_compare_means(size=2,comparisons = comps) + xlab('')
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = file.path(savedir,'MOCASubscoreBoxplotsbyClusterLouvain.pdf'),height = unit(4,'in'),width=unit(3.5,'in'),useDingbats = F)
plot(ggplot_gtable(p))
dev.off()

# to get df and effect size for overall score tests:

clusterNames <- sort(unique(partitionSample))
k <- length(clusterNames)
results <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
samp.size <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
pvals <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
latex <- matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames))
for(k1 in clusterNames){
  for(k2 in clusterNames){
    if(k1 != k2){m <- wilcox.test(meanMOCA[partitionSample == k1],meanMOCA[partitionSample == k2],conf.int = TRUE)}
    else if(k1==k2){m <- list(estimate=NA,p.value=NA)}
    
    results[k1,k2] <- m$estimate
    samp.size[k1,k2] <- sum(partitionSample %in% c(k1,k2))
    pvals[k1,k2] <- m$p.value
  }
}

pvals <- list.fdr.correct(pvals) # fdr correct
for(k1 in clusterNames){
  for(k2 in clusterNames){
    latex[k1,k2] <- paste0('$M_{',substr(k1,nchar(k1),nchar(k1)),'-',substr(k2,nchar(k2),nchar(k2)),'}=',
                                            round(results[k1,k2],2),'$, $n=',samp.size[k1,k2],'$, $p_','\\','mathrm{FDR}=',signif(pvals[k1,k2],2),'$')
  }
}
