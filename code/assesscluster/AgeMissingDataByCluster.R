rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleBraakCERAD.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics122219.csv',sep=''))
onset <- read.csv(paste0('data/INDD_GlobalOnset122219.csv'),stringsAsFactors = F)

demo <- demo[demo$INDDID %in% INDDIDs,]	# only look at patients in our sample
age.death <- demo$AgeatDeath # extract age.death at death
age.death[age.death == 0] <- NA # nobody is 0 years old

age.onset <- onset$GlobalAgeOnset

# missing onset data
dx.mis <- patientSample$NPDx1[!INDDIDs %in% onset$INDDID]
ggplot() + geom_bar(aes(x=as.character(dx.mis))) + theme(axis.text.x=element_text(angle=90))
######################
### Age by Cluster ###
######################

clusterColors <- getClusterColors(k)
partition.names <- sapply(partition, function(i) paste('Cluster',i))
names(partition.names) <- INDDIDs

age.data <- list(AgeAtDeath=data.frame(INDDID=INDDIDs,y=age.death),
                 AgeOnset=data.frame(INDDID=onset$INDDID,y=age.onset))
for(age.lab in names(age.data)){
  df <- data.frame(INDDID=INDDIDs,g=partition.names,stringsAsFactors = F)
  clusterNames <- sort(unique(df$g))
  comps <- combn(clusterNames,m = 2,simplify = FALSE)
  df <- merge(df,age.data[[age.lab]],by='INDDID')

  p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5) + theme_classic()+
    scale_y_continuous(breaks=seq(40,160,by = 10)) +
    scale_fill_manual(values=clusterColors,name='') + ylab('Age (y)') + xlab('') +  
    theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
    theme(text= element_text(size=8),
          axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
    stat_compare_means(size=2,comparisons = comps,position=5,method = "wilcox.test") 
  p <- mult.comp.ggpubr(p)
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(file = paste(savedir,'/',age.lab,'ByClusterPairwiseWilcox.pdf',sep=''),height = unit(3,'in'),width=unit(3,'in'),useDingbats = F)
  plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
  dev.off()
  
  p <- ggplot(df) + geom_violin(aes(x=g,y=y,fill=g)) + 
    scale_y_continuous(breaks=seq(20,110,by = 10)) +
  	scale_fill_manual(breaks = clusterNames, 
  		limits = sapply(1:k, function(i) paste('Cluster',i)),
  		values=clusterColors) + theme_classic() + xlab('') + ylab('Age') +
  	theme(legend.position= 'none',text=element_text(size=8),axis.text.x=element_text(angle=90))
  
  ggsave(filename = paste(savedir,age.lab,'ByCluster.pdf',sep=''),plot = p,
         height = 3,width=3,units='in')
}
###############################
### Missing data by cluster ###
###############################

missingdata <- rowMeans(is.na(microSample))*100
df <- data.frame(g=sapply(partition, function(i) paste('Cluster',i)),
                 y = missingdata,stringsAsFactors = F)
clusterNames <- sort(unique(df$g))
comps <- combn(clusterNames,m = 2,simplify = FALSE)

p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5) + theme_classic()+
  scale_y_continuous(breaks=seq(0,100,by = 10)) +
  scale_fill_manual(values=clusterColors,name='') + ylab('Missing Features (%)') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=2,comparisons = comps,position=5,method = "wilcox.test") 
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/MissingFeaturesByClusterPairwiseWilcox.pdf',sep=''),height = unit(3,'in'),width=unit(3,'in'),useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()

p <- ggplot(df) + geom_violin(aes(x=g,y=y,fill=g)) + 
  scale_y_continuous(breaks=seq(0,100,by = 10)) +
  scale_fill_manual(breaks = clusterNames, 
                    limits = sapply(1:k, function(i) paste('Cluster',i)),
                    values=clusterColors) + theme_classic() + xlab('') + ylab('Missing Features (%)') +
  theme(legend.position= 'none',text=element_text(size=8),axis.text.x=element_text(angle=90))

ggsave(filename = paste(savedir,'MissingFeaturesByCluster.pdf',sep=''),plot = p,
       height = 3,width=3,units='in')

###############################################
### Missing data by cluster by feature type ###
###############################################

pathItems.type <- list("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody")
pathItems.prettylab <- list("Neuron Loss","Gliosis","Angiopathy","Ubiquitin","Neuritic Plaques","TDP-43","Tau","Synuclein","Amyloid-beta")
pathItems.index <- sapply(1:length(pathItems.type), function(i) grep(pathItems.type[[i]], colnames(microSample)))
missing.feat.mask <- lapply(pathItems.index,length) == 0 # find feature types that have been entirely removed

# remove missing feature classes
pathItems.type <- pathItems.type[!missing.feat.mask]
pathItems.prettylab <- pathItems.prettylab[!missing.feat.mask]
pathItems.index <- pathItems.index[!missing.feat.mask]
nfeats <- length(pathItems.index)

# compute percentage.death of missing features within each feature type
missingdata <- lapply(pathItems.index, function(P) rowMeans(is.na(microSample[,P,drop=FALSE]))*100)
names(missingdata) <- pathItems.prettylab
missingdata <- do.call('cbind',missingdata)

# construct data frame containing missing features, cluster labels, and feature types
df <- data.frame(g=rep(sapply(partition, function(i) paste('Cluster',i)),ncol(missingdata)),
                 y = as.vector(as.matrix(missingdata)),
                 x=rep(colnames(missingdata),each = length(partition)),
                 stringsAsFactors = F)

comps <- combn(clusterNames,m = 2,simplify = FALSE)

p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5) + theme_classic()+
  facet_wrap(~x,scales='free') +
  scale_y_continuous(breaks=seq(0,100,by = 10)) +
  scale_fill_manual(values=clusterColors,name='') + ylab('Missing Features (%)') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=2,comparisons = comps,position=5,method = "wilcox.test") 
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/MissingFeaturesByClusterByTypePairwiseWilcox.pdf',sep=''),height = unit(9,'in'),width=unit(9,'in'),useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()

p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_jitter(alpha=0.5,stroke=0) + theme_classic()+
  facet_wrap(~x,scales='free') +
  scale_y_continuous(breaks=seq(0,100,by = 10)) +
  scale_fill_manual(values=clusterColors,name='') + ylab('Missing Features (%)') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=2,comparisons = comps,position=5,method = "wilcox.test") 
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/MissingFeaturesByClusterByTypeJitterPairwiseWilcox.pdf',sep=''),height = unit(9,'in'),width=unit(9,'in'),useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()

