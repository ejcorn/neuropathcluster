rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

#######################################################################
### Plot the full distribution of pathology scores for each cluster ###
#######################################################################

# extract scores for all patients in each cluster
microSample.byCluster <- lapply(1:k, function(k.i) microSample[partition==k.i,])
pathScores <- c('NA','1','2','3','4','5')
pathScores.pretty <- c('NA','0','Rare','1+','2+','3+')
# count occurrence of each score within each cluster
microSample.byCluster.counts <- lapply(microSample.byCluster,function(Y) 
  sapply(Y,function(X) count.ejc(X,items=pathScores)))
clusterCounts <- count.ejc(partition)
p <- lapply(1:k,function(k.i) 
  imagesc(100*microSample.byCluster.counts[[k.i]]/clusterCounts[[k.i]],
          caxis_name='%',clim=c(0,100),yticklabels=pathScores.pretty,cmap = 'Blues') + 
    ggtitle(paste('Cluster',k.i))+
    theme(axis.text.x = element_text(angle=90,size=4,hjust=0,vjust=0.5),
          plot.title = element_text(size=8,hjust=0.5))+
    nice_cbar(pos='bottom'))
p <- plot_grid(plotlist = p,nrow=2)       
ggsave(filename = paste0(savedir,'PathologyScoreDistributionsByClusters.pdf'),plot = p,
       height = 12,width=18,units='cm')

#####################################################################################
### Measure deviation from uniformity for each pathology measure for each cluster ###
#####################################################################################

# get chi-square stat for goodness of fit test for path scores for each feature for each cluster
# for each feature this quantifies the extent to which the distribution of 
# discrete scores deviates from uniformity

CV.cluster <- lapply(microSample.byCluster.counts, 
                     function(X) sapply(colnames(X),function(P) CramerV.GOF(chisq.test(X[,P]))))
clusterNames <- sapply(1:k,function(k.i) rep(paste('Cluster',k.i)))
names(CV.cluster) <- clusterNames
clusterAssignments <- lapply(1:k, function(k.i) rep(paste('Cluster',k.i),ncol(microSample)))

# get chi-square and cramer V for disease diagnoses
list[patientSample.othered,dz.short] <- other.dz(patientSample)
dzs.othered <- sort(unique(patientSample.othered$NPDx1)) # list of diseases, lumping together several diseases into other
dz.masks <- lapply(dzs.othered, function(dz) !exclude.dz(patientSample.othered,dz,5)) # find diseases in NPDx1-5
# extract scores for all patients with each disease label
microSample.byDisease <- lapply(dz.masks, function(dz.mask)
  microSample[dz.mask,])
# count occurrence of each score within each disease label
microSample.byDisease.counts <- lapply(microSample.byDisease,function(Y) 
  sapply(Y,function(X) count.ejc(X,items=pathScores)))
# quantify
CV.disease <- lapply(microSample.byDisease.counts, 
                     function(X) sapply(colnames(X),function(P) CramerV.GOF(chisq.test(X[,P]))))
names(CV.disease) <- dzs.othered
diseaseAssignments <- lapply(dzs.othered, function(d.i) rep(d.i,ncol(microSample)))

df.cluster <- data.frame(Diagnosis=unlist(clusterAssignments),CV=unlist(CV.cluster))
df.disease <- data.frame(Diagnosis=unlist(diseaseAssignments),CV=unlist(CV.disease))

p <- ggplot(rbind(df.cluster,df.disease)) + geom_boxplot(aes(x=Diagnosis,y=CV)) +
  theme_classic() + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
p

# compile all cramer V stats into one matrix
X <- rbind(do.call('rbind',CV.disease),do.call('rbind',CV.cluster))
colnames(X) <- substr(colnames(X),start=1,stop=nchar(colnames(X))-10)
p <- imagesc(X) + theme(axis.text.x= element_text(angle=90,hjust=1,vjust=0.5))

# perform targeted tests comparing uniformity of scores within each cluster to its constituent diseases
dzs <- c('Pick\'s disease','PSP','Alzheimer\'s disease','Amyotrophic Lateral Sclerosis','FTLD-TDP','Parkinson\'s disease','Multiple System Atrophy','LBD')
tests <- list(list(Cluster='Cluster 2',Diseases='Alzheimer\'s Disease'),
     list(Cluster='Cluster 1',Diseases=c('Pick\'s disease','PSP')),
     list(Cluster='Cluster 4',Diseases=c('Parkinson\'s disease','Multiple System Atrophy','LBD')),
     list(Cluster='Cluster 3',Diseases=c('Amyotrophic Lateral Sclerosis','FTLD-TDP')))

CV.diseasevscluster <- do.call('cbind',lapply(CV.cluster, function(X) 
  sapply(dzs,function(dz) wilcox.test(X-CV.disease[[dz]],conf.int = TRUE)$estimate)))
rownames(CV.diseasevscluster) <- dzs

CV.diseasevscluster.p <- do.call('cbind',lapply(CV.cluster, function(X) 
  sapply(dzs,function(dz) wilcox.test(X,CV.disease[[dz]],conf.int = TRUE)$p.value)))

melted_p <- melt(t(CV.diseasevscluster.p))
melted_p$value <- p.signif(melted_p$value)

clim <- c(-0.2,0.2)
p <- imagesc(CV.diseasevscluster,clim = clim,caxis_name = 'Cluster > Dz',cmap='redblue') +
  geom_text(data=melted_p,aes(x=Var1,y=Var2,label=value),size=2.5,color='white')
ggsave(filename = paste0(savedir,'PathologyScoreCramerVDiseaseVsCluster.pdf'),plot = p,
       height = 12,width=18,units='cm')
