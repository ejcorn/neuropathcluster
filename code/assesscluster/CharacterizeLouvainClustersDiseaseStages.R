# same as CharacterizeLouvainClusters.R but stratifies AD and LBD by disease stages

rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleABC.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

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
                       limits=c(1,5),breaks=c(1:5),labels=c('0','Rare','1+','2+','3+'),name='Score') +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,size=6,colour = clusterColors),
        axis.text.y = element_text(size = 6), axis.ticks = element_blank(),
        legend.key.size = unit(0.2,'cm'),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6),
        plot.margin = unit(c(0, 0, 0, 0), "null")) 
p2
ggsave(plot = p2,filename = paste(savedir,"CentroidColorbar.pdf",sep=''),width = 4.5, height = 5, units = "cm")

# add lines that separate each type of pathological feature
pathItems.type <- pathItems.labels[pathItems.labels!='']
sep.intercepts <- sapply(pathItems.type, function(ptype) 0.5 + max(which(grepl(ptype,rownames(t(fliplr(centroids)))))))
p2 <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + xlab("") + ylab("") + 
  geom_hline(yintercept = sep.intercepts,size=0.1)+
  scale_y_discrete(labels = pathItems.labels[length(pathItems.labels):1]) + 
  scale_x_discrete(labels = rownames(centroids),expand=c(0,0)) +
  #scale_fill_viridis(option = 'plasma') + 
  scale_fill_gradientn(colours = c('white','#779ecb','#283480','#0D0887')) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,size=6,colour = clusterColors),
        axis.text.y = element_text(size = 6), axis.ticks = element_blank(),
        legend.position = 'none',panel.border =element_rect(color = 'black',fill=NA))
# legend.key.size = unit(0.1,'cm'),legend.text = element_text(size=6),
# legend.title = element_text(size=6),legend.margin=ggplot2::margin(-5,-5,-5,-5),
# #legend.box.margin=ggplot2::margin(-1,-1,-1,-1),
# plot.margin = unit(c(0, 0, 0, 0), "null")) 
p2
ggsave(plot = p2,filename = paste(savedir,"LouvainClusterCentroids.pdf",sep=''),width = 1.5, height = 2, units = "in")

save(centroids,file = paste(params$sourcedata.dir,'Fig2d_SourceData.RData',sep=''))

############################
### Cluster by diagnoses ###
############################

list[patientSample,dz.short]<- other.dz(patientSample)

# add in graded colors for AD stages -- first get colors for each disease
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.short))
AD.Color <- dz.colors[which(dz.short == 'AD')] # then get AD color
LBD.Color <- dz.colors[which(dz.short == 'LBD')] # and get LBD color
# then make 3 colors for int high that get progressively lighter (because high, int is alphabetical and dz.short must match sort(unique(patientSample$NPDx1)))
AD.Stage.Colors <- c(darken(col = AD.Color,amount = 0.1,method='relative'),AD.Color,lighten(col = AD.Color,amount = 0.3,method='relative'))
#AD.Stage.Colors <- AD.Stage.Colors[1:length(unique(patie))]
# insert new colors
dz.colors <- c(dz.colors[0:(which(dz.short == 'AD')-1)],AD.Stage.Colors,dz.colors[(which(dz.short == 'AD')+1):length(dz.colors)])
# replace AD diagnoses with ADNC levels
patientSample$NPDx1[patientSample$NPDx1 %in% c('Alzheimer\'s disease')] <- paste('ADNC -',patientSample$ADStatus[patientSample$NPDx1 == 'Alzheimer\'s disease'])
# if ADNC - None shows up as primary dx... that means Alzheimer's disease listed but ABC criteria not met ... so AD listing probably a mistake AND probably no other pertinent diagnoses
# this only occurs in one patient and they are essentially normal except with very mild tau
patientSample$NPDx1[patientSample$NPDx1 == 'ADNC - None'] <- 'Unremarkable adult' 
list[patientSample,dz.short]<- other.dz(patientSample)
# make 3 colors for brainstem limbic neocortical that get darker (b/c alphabetical order, see above)
LBD.Stage.Colors <- c(lighten(col = LBD.Color,amount = 0.3,method='relative'),lighten(col = LBD.Color,amount = 0.2,method='relative'),LBD.Color,darken(col=LBD.Color,amount=0.1,method='relative'))
# insert new colors
dz.colors <- c(dz.colors[0:(which(dz.short == 'LBD')-1)],LBD.Stage.Colors,dz.colors[(which(dz.short == 'LBD')+1):length(dz.colors)])
# replace LBD diagnoses with LBD patterns
patientSample$NPDx1[patientSample$NPDx1 == 'LBD'] <- paste('LBD -',patientSample$DLBType[patientSample$NPDx1 == 'LBD'])

list[patientSample,dz.short]<- other.dz(patientSample)
list[p.k.dz,df.plt2e] <- plot.dz.by.cluster(patientSample$NPDx1,partition,dz.short,clusterColors,'Primary\nHistopathologic Diagnosis',dz.colors)
ggsave(filename = paste(savedir,'ClustersByPrimaryDiseaseStages.pdf',sep=''),plot = p.k.dz,
       height = 2,width=3,units='in')
save(df.plt2e,file = paste(params$sourcedata.dir,'Fig2e_SourceData.RData',sep=''))

# now isolate patients with high or intermediate ADNC and plot by secondary diagnoses
#AD.mask <- patientSample$ADStatus %in% c('High','Intermediate') 
AD.mask <- patientSample$NPDx1 %in% c('ADNC - Intermediate','ADNC - High')
partition.AD <- partition[AD.mask]
patientSample.AD <- patientSample[AD.mask,]
patientSample.AD$NPDx2[patientSample.AD$NPDx2 == ''] <- 'None' # label empty NPDx2 as none
patientSample.AD$NPDx2[grep('Lewy',patientSample.AD$NPDx2)] <- 'LBD'
patientSample.AD$NPDx2[patientSample.AD$NPDx2 == 'LBD'] <- paste('LBD -',patientSample.AD$DLBType[patientSample.AD$NPDx2 == 'LBD'])
list[patientSample.AD,dz.short.AD]<- other.dz(patientSample.AD,NPDx='NPDx2')
list[p.k.dz2,df.plt2f] <- plot.dz.by.cluster(patientSample.AD$NPDx2,partition.AD,dz.short.AD,clusterColors,'Secondary\nHistopathologic Diagnosis')

# reorder legend labels and colors to match Fig. 2e
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.order <- as.character(unique(df.plt2e$Dz)) # use disease labels from fig 2e
#dz.order[which(dz.order == 'Unremarkable adult')] <- 'None' # replace unremarkable adult with none so they have the same color
#dz.order <- c(dz.order,'None') # make it the last entry in legend
#dz.colors <- pal(length(dz.order)+1)
#dz.colors <- pal(length(dz.order))

dz.colors <- dz.colors[-which(grepl('ADNC',dz.order))]
dz.order <- dz.order[-which(grepl('ADNC',dz.order))]
# insert amygdala LBD
dz.colors <- c(dz.colors[1:which(grepl('LBD - NOS',dz.order))],lighten(col = LBD.Color,amount = 0.3,method='relative'),dz.colors[(which(grepl('LBD - NOS',dz.order))+1):length(dz.colors)])
dz.order <- c(dz.order[1:which(grepl('LBD - NOS',dz.order))],'LBD - Amygdala',dz.order[(which(grepl('LBD - NOS',dz.order))+1):length(dz.order)])
dz.order <- c(dz.order,'None') # add none as last entry
dz.colors <- c(dz.colors,'grey80') # add blank color
p.k.dz2 <- p.k.dz2 + scale_fill_manual(limits=dz.order, values = dz.colors,name='Secondary\nHistopathologic Diagnosis')

save(df.plt2f,file = paste(params$sourcedata.dir,'Fig2f_SourceData.RData',sep=''))
p.all <- plot_grid(plotlist = list(p2,p.k.dz,p.k.dz2),align = 'hv',nrow=1,axis = 'b',rel_widths = c(1,2,2))
ggsave(filename = paste(savedir,'ClusterCentroidsDiseaseStageNPDx1ADNPDx2.pdf',sep=''),plot = p.all,
       height = 5,width=18,units='cm')
