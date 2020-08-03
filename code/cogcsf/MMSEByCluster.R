# relationship between path scores and cognition
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/trainfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1)] # Get rid of index column only
load(paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
microSample.comp <- complete(microSample.imp)
colnames(microSample.comp) <- colnames(microSample)

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
microSample.comp <- remove.Disconnected.Subjects(microSample.comp,DisconnectedSubjects)

if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)

# load MMSE data
data.date <- '122219'
demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics',data.date,'.csv',sep=''),stringsAsFactors = F)
mmse <- read.csv(paste0(params$homedir,'data/INDD_MMSE',data.date,'.csv'),stringsAsFactors=F)
subj.scores <- lapply(INDDIDs, function(id) mmse[mmse$INDDID == id,]) # get each subject's scores 
age.at.death <- sapply(INDDIDs, function(id) demo[demo$INDDID ==id,'AgeatDeath']) # get subjects age at death
mmse <- sapply(subj.scores, function(df) ifelse(nrow(df)>0,yes=df$MMSETotal[which.max(df$AgeatTest)],no=NaN)) # either NaN if no scores or take most recent
age.at.test <- sapply(subj.scores, function(df) ifelse(nrow(df)>0,yes=max(df$AgeatTest),no=NaN))
mmse[mmse>30] <-NA
test.minus.death <- age.at.test - age.at.death
test.minus.death[test.minus.death > 0] <- NaN # can't get an MMSE after you die
# plot MMSE scores by cluster
df <- data.frame(g=paste('Cluster',partition),y=mmse,stringsAsFactors = F)
clusterNames <- sort(unique(df$g))
clusterColors <- getClusterColors(k)
comps <- combn(clusterNames,m = 2,simplify = FALSE)

p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5,size=0.25) + theme_classic()+
  scale_y_continuous(breaks = c(0,10,20,30),limits=c(0,NA)) + ggtitle('MMSE')+ 
  scale_fill_manual(values=clusterColors,name='') + ylab('MMSE Total') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),plot.title=element_text(size=8,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=1.75,comparisons = comps,position=5,method = "wilcox.test") 
p <- mult.comp.ggpubr(p,method='fdr')
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/MMSEByClusterPairwiseWilcox.pdf',sep=''),height = 4.5/2.54,width=4.5/2.54,useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()

# plot age at test relative to death by cluster

df <- data.frame(g=paste('Cluster',partition),y=test.minus.death,stringsAsFactors = F)
clusterNames <- sort(unique(df$g))
clusterColors <- getClusterColors(k)
comps <- combn(clusterNames,m = 2,simplify = FALSE)

p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5,size=0.25) + theme_classic()+
  scale_y_continuous() + ggtitle('MMSE')+ 
  scale_fill_manual(values=clusterColors,name='') + ylab('Last MMSE - Death (years)') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),plot.title=element_text(size=8,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=1.75,comparisons = comps,position=5,method = "wilcox.test") 
p <- mult.comp.ggpubr(p,method='fdr')
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/LastMMSEMinusDeathByClusterPairwiseWilcox.pdf',sep=''),height = 4.5/2.54,width=4.5/2.54,useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()
