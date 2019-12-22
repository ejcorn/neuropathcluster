rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleBraakCERAD.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

######################
### Plot centroids ###
######################
clusterColors <- getClusterColors(k)
melted_cormat <- melt(centroids)
melted_cormat$Var2 <- melted_cormat$Var2[length(melted_cormat$Var1):1]

p2 <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + xlab("") + ylab("") +
  scale_y_discrete(labels = pathItems.labels[length(pathItems.labels):1]) + 
  scale_x_discrete(labels = rownames(centroids),expand=c(0,0)) +
  #scale_fill_viridis(option = 'plasma',limits=c(1,5),breaks=c(1:5)) + 
  scale_fill_gradientn(colours = c('white','#779ecb','#283480',"#0D0887"),
    limits=c(1,5),breaks=c(1:5)) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,size=6,colour = clusterColors),
        axis.text.y = element_text(size = 6), axis.ticks = element_blank(),
        legend.key.size = unit(0.3,'cm'),legend.text = element_text(size=6),
        legend.title = element_text(size=6),
        plot.margin = unit(c(0, 0, 0, 0), "null")) 
p2
ggsave(plot = p2,filename = paste(savedir,"CentroidColorbar.pdf",sep=''),width = 1.5, height = 2.25, units = "in")

p2 <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + xlab("") + ylab("") +
  scale_y_discrete(labels = pathItems.labels[length(pathItems.labels):1]) + 
  scale_x_discrete(labels = rownames(centroids),expand=c(0,0)) +
  #scale_fill_viridis(option = 'plasma') + 
  scale_fill_gradientn(colours = c('white','#779ecb','#283480','#0D0887')) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,size=6,colour = clusterColors),
        axis.text.y = element_text(size = 6), axis.ticks = element_blank(),
        legend.position = 'none')
        # legend.key.size = unit(0.1,'cm'),legend.text = element_text(size=6),
        # legend.title = element_text(size=6),legend.margin=ggplot2::margin(-5,-5,-5,-5),
        # #legend.box.margin=ggplot2::margin(-1,-1,-1,-1),
        # plot.margin = unit(c(0, 0, 0, 0), "null")) 
p2
ggsave(plot = p2,filename = paste(savedir,"LouvainClusterCentroids.pdf",sep=''),width = 1.5, height = 2, units = "in")

save(centroids,file = paste(savedir,'Fig2d_SourceData.RData',sep=''))

############################
### Cluster by diagnoses ###
############################

list[patientSample,dz.short]<- other.dz(patientSample)

list[p.k.dz,df.plt] <- plot.dz.by.cluster(patientSample$NPDx1,partition,dz.short,clusterColors,'Primary\nHistopathologic Diagnosis')
ggsave(filename = paste(savedir,'ClustersByPrimaryDiseaseLouvain.pdf',sep=''),plot = p.k.dz,
       height = 2,width=3,units='in')
save(df.plt,file = paste(savedir,'Fig2e_SourceData.RData',sep=''))

# now isolate AD patients and plot by secondary diagnoses
AD.mask <- patientSample$NPDx1 == 'Alzheimer\'s disease'
partition.AD <- partition[AD.mask]
patientSample.AD <- patientSample[AD.mask,]
patientSample.AD$NPDx2 <- as.character(patientSample.AD$NPDx2)
patientSample.AD$NPDx2[patientSample.AD$NPDx2 == ''] <- 'None' # label empty NPDx2 as none
patientSample.AD$NPDx2[grep('Lewy',patientSample.AD$NPDx2)] <- 'LBD'
list[patientSample.AD,dz.short.AD]<- other.dz(patientSample.AD,NPDx='NPDx2')
list[p.k.dz2,df.plt] <- plot.dz.by.cluster(patientSample.AD$NPDx2,partition.AD,dz.short.AD,clusterColors,'Secondary\nHistopathologic Diagnosis')

save(df.plt,file = paste(savedir,'Fig2f_SourceData.RData',sep=''))

############################
### Diagnoses by cluster ###
############################

dx <- data.frame(NPDx1 = patientSample$NPDx1,stringsAsFactors = F)
cl <- as.data.frame(as.factor(partition))
cl <- dummy.data.frame(data=cl)
colnames(cl) <- sapply(1:k, function(i) paste('Cluster',i))
dz.names <- as.character(sort(unique(dx$NPDx1)))
dz.k <- lapply(dz.names, function(d.i) colMeans(cl[dx==d.i,]))
dz.k <- do.call('cbind',dz.k)

df.plt <- data.frame(pr = as.vector(dz.k),
                     Cl = rep(rownames(dz.k),length(dz.names)),
                     Dz = as.vector(sapply(dz.short, function(d.i) rep(d.i,k))))

n.by.dz <- count.ejc(dx$NPDx1)

p.dz.k <- ggplot(data=df.plt) + geom_col(aes(x=Dz,y=pr,fill=Cl),position=position_stack(1)) +
  annotate(geom='text',x=dz.short,y=rep(1.1,length(dz.short)),label=paste('',n.by.dz,sep='\n'),size=1.5)+
  scale_y_continuous(expand=c(0,0)) + xlab('') + ylab('Proportion') +
  scale_fill_manual(values = clusterColors,name='Disease Cluster') + theme_classic() + 
  theme(axis.text.x = element_text(size=6,angle=90,vjust=0.5),
        legend.key.size = unit(0.01,'in'),
        text=element_text(size=6))
p.dz.k

ggsave(filename = paste(savedir,'DiseasesByClusterLouvain.pdf',sep=''),plot = p.dz.k,
       height = 4,width=6,units='in')

p.all <- plot_grid(plotlist = list(p2,p.k.dz,p.k.dz2),align = 'hv',nrow=1,axis = 'b',rel_widths = c(1,2,2))
ggsave(filename = paste(savedir,'LouvainCentroidsDiseaseNPDx1ADNPDx2.pdf',sep=''),plot = p.all,
       height = 5,width=18,units='cm')