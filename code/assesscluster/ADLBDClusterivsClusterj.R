rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'target_diagnoses/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)

##############################################
### Compare cluster i and cluster j AD-LBD ###
##############################################

# choose two clusters to compare the patients within them that have both AD and LBD

ci.cj.list <- list(list(C.i = 2,C.j=4),list(C.i = 2, C.j = 5),list(C.i = 4, C.j = 5))
for(ci.cj in ci.cj.list){
  C.i <- ci.cj$C.i
  C.j <- ci.cj$C.j
  # select patients with both AD and LBD in any diagnosis, in either order
  dzs <- c('LBD','Alzheimer\'s disease')
  LBD <- !exclude.dz(patientSample = patientSample,dz.exc = 'LBD',n.dx = 5)
  AD <- !exclude.dz(patientSample = patientSample,dz.exc = 'Alzheimer\'s disease',n.dx = 5)
  dz.mask <- AD & LBD
  
  # now hone in on two AD-LBD patients in two clusters
  # impute missing data
  load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
  microSample.comp <- complete(microSample.imp,1) # temporary
  colnames(microSample.comp) <- colnames(microSample) # temporary
  microSample <- microSample.comp
  
  microSample.Ci <- microSample[partition==C.i & dz.mask,]
  microSample.Cj <- microSample[partition==C.j & dz.mask,]
  w.tests <- list()
  for(P in colnames(microSample.Ci)){
    w.tests[[P]] <- tryCatch(wilcox.test(microSample.Ci[,P],microSample.Cj[,P],conf.int = T),
                             error = function(e){return(list(estimate=0,p.value=NA))})
  }
  #w.tests <- lapply(colnames(microSample.Ci),function(P) wilcox.test(microSample.Ci[,P],microSample.Cj[,P],conf.int = T))
  diffs <- sapply(w.tests,function(w) w$estimate)
  # adjust for multiple comparisons
  pvals <- sapply(w.tests,function(w) w$p.value)
  ggplot() + geom_col(aes(y=diffs,x=colnames(microSample.Ci)))+ theme(axis.text.x = element_text(angle=90))
  
  p2 <- imagesc(microSample.Ci) + ggtitle(paste('Cluster',C.i,'AD-Syn')) + 
    theme(axis.text.x = element_text(angle=90),plot.title = element_text(hjust=0.5))
  p4 <- imagesc(microSample.Cj) + ggtitle(paste('Cluster',C.j,'AD-Syn')) + 
    theme(axis.text.x = element_text(angle=90),plot.title = element_text(hjust=0.5))
  
  p <- plot_grid(plotlist = list(p2,p4))
  ggsave(plot = p,filename = paste0(savedir,'Cluster',C.i,'vs',C.j,'ADSynPathology.pdf'),
         width = 18, height = 6, units = "cm")
  
  pathScores <- c('0','Rare','1+','2+','3+')
  p2 <- imagesc(region.by.item.matrix(microSample.Ci,fxn = 'median'),cmap='Blues',clim=c(1,5),caxis_labels = pathScores) + ggtitle(paste('Cluster',C.i,'AD-Syn')) + 
    theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),plot.title = element_text(hjust=0.5)) +
    nice_cbar(pos='right')
  p4 <- imagesc(region.by.item.matrix(microSample.Cj,fxn = 'median'),cmap='Blues',clim=c(1,5),caxis_labels = pathScores) + ggtitle(paste('Cluster',C.j,'AD-Syn')) + 
    theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),plot.title = element_text(hjust=0.5)) +
    nice_cbar(pos='right')
  
  diffs <- matrix(diffs,nrow=1,ncol=ncol(microSample),dimnames = list(NULL,colnames(microSample)))
  pvals <- matrix(pvals,nrow=1,ncol=ncol(microSample),dimnames = list(NULL,colnames(microSample)))
  
  # posthoc correct over only thio, tau, synuclein
  posthoc.mask <- Reduce('|',lapply(c('_Thio','_Syn','_Tau'),function(j) grepl(j,colnames(pvals))))
  pvals[posthoc.mask] <- p.adjust(pvals[posthoc.mask],method='fdr')
  pvals[!posthoc.mask] <- 1
  melted_p <- melt(t(region.by.item.matrix(pvals)))
  melted_p$value <- p.signif(melted_p$value)
  #melted_p$value <- p.signif(p.adjust(melted_p$value, method = 'fdr'))
  melted_p$value[melted_p$value!='ns'] <- '*'
  melted_p$value[melted_p$value=='ns'] <- ''
  
  clim <- c(-2,2)
  p3<- imagesc(region.by.item.matrix(diffs,fxn='median'),cmap='redblue',clim = clim,caxis_name = '') +ggtitle(paste('Cluster',C.i,'- Cluster',C.j)) + 
    theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),plot.title = element_text(hjust=0.5)) +
    geom_text(data=melted_p,aes(x=Var1,y=Var2,label=value),size=2.5,color='white') +
    nice_cbar(pos='right')
  
  p <- plot_grid(plotlist = list(p2,p4,p3),nrow=1)
  ggsave(plot = p,filename = paste0(savedir,'Cluster',C.i,'vs',C.j,'MedianADSynPathology.pdf'),
         width = 18, height = 4, units = "cm")
  
}