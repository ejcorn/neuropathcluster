rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/patientcharacteristics/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleABC.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

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
                 AgeOnset=data.frame(INDDID=onset$INDDID,y=age.onset),
                 # add disease duration
                 AutopsyDate=data.frame(INDDID=INDDIDs,y=as.Date(patientSample$AutopsyDate)))
date.breaks.p <- seq.Date(as.Date('1990/01/01'),as.Date('2020/01/01'),by = '5 years')

for(age.lab in names(age.data)){
  df <- data.frame(INDDID=INDDIDs,g=partition.names,stringsAsFactors = F)
  clusterNames <- sort(unique(df$g))
  comps <- combn(clusterNames,m = 2,simplify = FALSE)
  df <- merge(df,age.data[[age.lab]],by='INDDID')

  p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5) + theme_classic()+
    scale_y_continuous(breaks=seq(20,110,by = 20),limits=c(10,NA)) + ggtitle(age.lab)+ 
    scale_fill_manual(values=clusterColors,name='') + ylab('Age (y)') + xlab('') +  
    theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
    theme(text= element_text(size=8),plot.title=element_text(size=8,hjust=0.5),
          axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
    stat_compare_means(size=1.75,comparisons = comps,position=5,method = "wilcox.test") 
  if(is.Date(df$y)){p <- p + scale_y_date(breaks=date.breaks.p,labels=year(date.breaks.p))} # if autopsy date add special breaks for y axis
  p <- mult.comp.ggpubr(p,method='fdr')
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(file = paste(savedir,'/',age.lab,'ByClusterPairwiseWilcox.pdf',sep=''),height = unit(6/2.54,'in'),width=unit(4.5/2.54,'in'),useDingbats = F)
  plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
  dev.off()
  
  p <- ggplot(df) + geom_violin(aes(x=g,y=y,fill=g)) + 
    scale_y_continuous(breaks=seq(20,110,by = 20),limits=c(10,NA)) +
  	scale_fill_manual(breaks = clusterNames, 
  		limits = sapply(1:k, function(i) paste('Cluster',i)),
  		values=clusterColors) + theme_classic() + xlab('') + ylab('Age') + ggtitle(age.lab)+
  	theme(legend.position= 'none',text=element_text(size=8),axis.text.x=element_text(angle=90),plot.title=element_text(size=8,hjust=0.5))
  if(is.Date(df$y)){p <- p + scale_y_date(breaks=date.breaks.p,labels=year(date.breaks.p))} # if autopsy date add special breaks for y axis
  ggsave(filename = paste(savedir,age.lab,'ByCluster.pdf',sep=''),plot = p,
         height = 5,width=4.5,units='cm')
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
p <- mult.comp.ggpubr(p,method='none')
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
p <- mult.comp.ggpubr(p,method='none') # don't multiple comparisons correct
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/MissingFeaturesByClusterByTypePairwiseWilcox.pdf',sep=''),height = unit(9,'in'),width=unit(9,'in'),useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()

p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_jitter(alpha=0.5,stroke=0,size=0.5) + theme_classic()+
  facet_wrap(~x,scales='free',nrow=2) +
  scale_y_continuous(breaks=seq(0,100,by = 20)) +
  scale_fill_manual(values=clusterColors,name='') + ylab('Missing Features (%)') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=1.5,comparisons = comps,position=5,method = "wilcox.test",label.y=seq(100,275,length.out=length(comps)))
p <- mult.comp.ggpubr(p,method='none')
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/MissingFeaturesByClusterByTypeJitterPairwiseWilcox.pdf',sep=''),height = unit(10/2.54,'in'),width=unit(18/2.54,'in'),useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()

# plot sex ratio by cluster
sex <- as.character(demo$Sex)
female.overall <- 100*mean(sex=='Female')
female.by.cluster <- sapply(1:k, function(k.i) 100*mean(sex[partition==k.i]=='Female'))
p <- ggplot() + geom_col(aes(x=clusterNames,y=female.by.cluster,fill=clusterNames))+
  geom_hline(yintercept = female.overall,linetype='dashed',color='grey70')+
  theme_classic() + scale_fill_manual(values=getClusterColors(k)) + 
  scale_y_continuous(expand=c(0,0),limits=c(0,50),breaks=seq(0,50,by=10))+
  theme(text= element_text(size=8),legend.position = 'none',
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) +
  ylab('% Female') + xlab('')
ggsave(filename = paste0(savedir,'SexByCluster.pdf'),plot = p,
       height = 6,width=4.5,units='cm')

# make table and save as latex code in .txt file
female.by.cluster <- sapply(1:k, function(k.i) 100*mean(sex[partition==k.i]=='Female'))
names(female.by.cluster) <- clusterNames
female.by.cluster['Overall'] <- female.overall
sex.table <- t(as.data.frame(female.by.cluster))
rownames(sex.table) <- '% Female'
lt <- xtable(x=sex.table,caption = 'Table 1. Sex by cluster and in overall sample.',label='table:sexbycluster')
write.table(x=print(lt),file = paste0(savedir,'SexByClusterTable.txt'),row.names = F,col.names = F)

