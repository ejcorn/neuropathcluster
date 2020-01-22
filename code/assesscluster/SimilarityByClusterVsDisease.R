# purpose of this script is to compare similarity of pathology within clusters vs. within diseases
# similarity is measured by the given distance metric

rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/comparesimilarity_',params$dist.met,'/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)

#######################################################
### Get similarity matrix for given distance metric ###
#######################################################

if(params$dist.met == 'spearman'){
  W <- cor(t(microSample),method = 'spearman')
} else if(params$dist.met == 'polychoric'){
  load(file = paste0(params$opdir,'processed/SubjectPolychoricMatrix.RData'))
  W <- subject.cormat$rho
  #W <- thresh.mat(W,'>',thresh.bestsil)
}

################################################
### Define masks for each disease or cluster ###
################################################

list[patientSample,dz.short]<- other.dz(patientSample)
dz.masks <- lapply(sort(unique(patientSample$NPDx1)), function(dx.i) 
  !exclude.dz(patientSample = patientSample,dz.exc = dx.i,n.dx = 5)) # get masks for each main disease
dz.masks <- c(dz.masks, lapply(1:k, function(k.i) partition==k.i)) # add cluster masks
dz.short <- c(dz.short,sapply(1:k, function(k.i) paste('Cluster',k.i))) # add cluster names

###################################################
### Compute and compare within group similarity ###
###################################################

Within.Group.Similarity <- lapply(dz.masks, function(M) rowMedians(W[M,M])) # extract similarity values within each group (on-diagonal blocks)
names(Within.Group.Similarity) <- dz.short
Within.Group.Labels <- lapply(dz.short, function(dx.i) rep(dx.i,length(Within.Group.Similarity[[dx.i]])))
names(Within.Group.Labels) <- dz.short
df.plt <- data.frame(Similarity=unlist(Within.Group.Similarity),Group=unlist(Within.Group.Labels))
p <- ggplot(df.plt) + geom_boxplot(aes(x=Group,y=Similarity)) + theme_classic() +
  scale_x_discrete(limits=dz.short) +
  theme(text=element_text(size=8),axis.text.x = element_text(angle=90)) + xlab('')

# check that you extracted the right on-diagonal blocks for each disease
all(unname(sapply(Within.Group.Similarity,function(x) sqrt(length(x))))==sapply(dz.masks,sum))

# perform targeted tests comparing similarity within each cluster to its constituent diseases
tests <- list(list(Cluster='Cluster 2',Diseases='AD'),
              list(Cluster='Cluster 5',Diseases=c('PiD','PSP','CBD')),
              list(Cluster='Cluster 3',Diseases=c('PD','MSA','LBD')),
              list(Cluster='Cluster 4',Diseases=c('ALS','FTLD')),
              list(Cluster='Cluster 6',Diseases=c('LBD','AD')))

for(test in tests){
  dz.test <- unlist(test) # compile cluster names and disease names to be tested into one vector
  df <- data.frame(g=unname(unlist(Within.Group.Labels[dz.test])),
                   y = unname(unlist(Within.Group.Similarity[dz.test])),stringsAsFactors = F)
  groupNames <- sort(unique(df$g))
  comps <- combn(groupNames,m = 2,simplify = FALSE)
  groupColors <- brewer.pal(9,'Blues')[seq(3,9,length.out=length(groupNames))]
  p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5) + theme_classic()+
    scale_y_continuous(limits=c(-0.25,2),breaks=seq(-0.25,0.5,length.out=4)) +
    scale_x_discrete(limits=dz.test) +
    scale_fill_manual(values=groupColors,name='') + ylab('Similarity') + xlab('') +
    ggtitle('Within Group') + theme(plot.title = element_text(hjust=0.5))+
    theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
    theme(text= element_text(size=8),
          axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = groupColors)) + 
    stat_compare_means(size=2,comparisons = comps,position=5,method = "wilcox.test") 
  p <- mult.comp.ggpubr(p)
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(file = paste0(savedir,'/WithinGroupSimilarity',paste(dz.test,collapse = ','),'.pdf'),height = unit(3,'in'),width=unit(3,'in'),useDingbats = F)
  plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
  plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
  dev.off()
}

####################################################
### Compute and compare between group similarity ###
####################################################

# extract similarity values between each group and all other subjects (bands excluding diagonal block)
Between.Group.Similarity <- lapply(dz.masks, function(M) rowMedians(W[which(M),-which(M)]))
names(Between.Group.Similarity) <- dz.short
Between.Group.Labels <- lapply(dz.short, function(dx.i) rep(dx.i,length(Between.Group.Similarity[[dx.i]])))
names(Between.Group.Labels) <- dz.short

# perform targeted tests comparing similarity within each cluster to its constituent diseases
tests <- list(list(Cluster='Cluster 2',Diseases='AD'),
              list(Cluster='Cluster 5',Diseases=c('PiD','PSP','CBD')),
              list(Cluster='Cluster 3',Diseases=c('PD','MSA','LBD')),
              list(Cluster='Cluster 4',Diseases=c('ALS','FTLD')),
              list(Cluster='Cluster 6',Diseases=c('LBD','AD')))

for(test in tests){
  dz.test <- unlist(test) # compile cluster names and disease names to be tested into one vector
  df <- data.frame(g=unname(unlist(Between.Group.Labels[dz.test])),
                   y = unname(unlist(Between.Group.Similarity[dz.test])),stringsAsFactors = F)
  groupNames <- sort(unique(df$g))
  comps <- combn(groupNames,m = 2,simplify = FALSE)
  groupColors <- brewer.pal(9,'Blues')[seq(3,9,length.out=length(groupNames))]
  p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5) + theme_classic()+
    scale_y_continuous(limits=c(-0.25,2),breaks=seq(-0.25,0.5,length.out=4)) +
    scale_x_discrete(limits=dz.test) +
    scale_fill_manual(values=groupColors,name='') + ylab('Similarity') + xlab('') +  
    ggtitle('Between Group') + theme(plot.title = element_text(hjust=0.5))+
    theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
    theme(text= element_text(size=8),
          axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = groupColors)) + 
    stat_compare_means(size=2,comparisons = comps,position=5,method = "wilcox.test") 
  p <- mult.comp.ggpubr(p)
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(file = paste0(savedir,'/BetweenGroupSimilarity',paste(dz.test,collapse = ','),'.pdf'),height = unit(3,'in'),width=unit(3,'in'),useDingbats = F)
  plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
  plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
  dev.off()
}

