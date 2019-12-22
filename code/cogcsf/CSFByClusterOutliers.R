rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'csfcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

clusterColors <- getClusterColors(k)

###########
### CSF ###
###########

CSF.name <- 'Luminex'
CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% INDDIDs,]
CSF.vars <- c('LuminexTTau','LuminexPTau','LuminexAbeta42')

CSF.by.pt <- lapply(sort(unique(CSF$INDDID)), function(id) # order tests by date w/in subjs
  CSF[CSF$INDDID==id,c('CSFDate',CSF.vars)])
CSF.by.pt <- lapply(CSF.by.pt, function(X) X[order(X$CSFDate),-1])
#CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) colMeans(X)))
CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) X[1,])) # just use first sample
# for some patients, could compute a feature based on change in CSF proteins
#CSF.diff <- do.call('rbind',lapply(CSF.by.pt, function(X) X[nrow(X),] - X[1,]))
CSF.mean <- cbind(data.frame(INDDID=sort(unique(CSF$INDDID))),CSF.mean)
#CSF.mean <- CSF.mean[-(which(CSF.mean$LuminexTTau >1000)),] # ttau outliers

partitionSample <- as.character(partition[INDDIDs %in% unique(CSF.mean$INDDID)])
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))

CSF.mean <- CSF.mean[,-1]

# get sample size and effect size

clusterNames <- sort(unique(partitionSample))
results <- samp.size <- p.vals <- list()
for(CSF.protein in colnames(CSF.mean)){
  results[[CSF.protein]] <- samp.size[[CSF.protein]] <- p.vals[[CSF.protein]] <- 
    matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames)) 
  # store data in new variables
  CSF.mean.test <- CSF.mean
  partitionSample.test <- partitionSample
  for(k1 in clusterNames){
    for(k2 in clusterNames){      
      m <- wilcox.test(CSF.mean.test[partitionSample.test == k1,CSF.protein],CSF.mean.test[partitionSample.test == k2,CSF.protein],conf.int = TRUE)
      results[[CSF.protein]][k1,k2] <- m$estimate
      samp.size[[CSF.protein]][k1,k2] <- sum(partitionSample.test %in% c(k1,k2))
      p.vals[[CSF.protein]][k1,k2] <- m$p.value
    }
  }  
}
p.vals <- list.fdr.correct(p.vals)
min.beta <- min(results[['LuminexTTau']])
max.beta <- max(results[['LuminexTTau']])
p.full <- plot.allele.beta.matrix(results[['LuminexTTau']],p.vals[['LuminexTTau']],'Total Tau','Full Sample',min.beta,max.beta,clusterColors)
list[bmat,pmat] <- list(results[['LuminexTTau']],p.vals[['LuminexTTau']])
save(bmat,pmat,file = paste(savedir,'FigS9a_left_SourceData.RData',sep=''))


ymax <- 1000
clusterNames <- sort(unique(partitionSample))
results <- samp.size <- p.vals <- list()
for(CSF.protein in colnames(CSF.mean)){
  results[[CSF.protein]] <- samp.size[[CSF.protein]] <- p.vals[[CSF.protein]] <- 
    matrix(NA,ncol = k, nrow = k,dimnames = list(clusterNames,clusterNames)) 
  # store data in new variables
  CSF.mean.test <- CSF.mean
  partitionSample.test <- partitionSample
  # remove 2 outliers for total tau
  if(CSF.protein == 'LuminexTTau'){
    outlier.mask <- CSF.mean$LuminexTTau < ymax
    CSF.mean.test <- CSF.mean[outlier.mask,]
    partitionSample.test <- partitionSample[outlier.mask]
  }
  for(k1 in clusterNames){
    for(k2 in clusterNames){      
      m <- wilcox.test(CSF.mean.test[partitionSample.test == k1,CSF.protein],CSF.mean.test[partitionSample.test == k2,CSF.protein],conf.int = TRUE)
      results[[CSF.protein]][k1,k2] <- m$estimate
      samp.size[[CSF.protein]][k1,k2] <- sum(partitionSample.test %in% c(k1,k2))
      p.vals[[CSF.protein]][k1,k2] <- m$p.value
    }
  }  
}
p.vals <- list.fdr.correct(p.vals)
min.beta <- min(results[['LuminexTTau']])
max.beta <- max(results[['LuminexTTau']])
p.outlier <- plot.allele.beta.matrix(results[['LuminexTTau']],p.vals[['LuminexTTau']],'Total Tau','No Outliers',min.beta,max.beta,clusterColors)
list[bmat,pmat] <- list(results[['LuminexTTau']],p.vals[['LuminexTTau']])
save(bmat,pmat,file = paste(savedir,'FigS9a_right_SourceData.RData',sep=''))

p.all <- plot_grid(plotlist = list(p.full,p.outlier),nrow=1,align='hv')
ggsave(filename = paste(savedir,'TotalTauFullVsOutliersRemoved.pdf',sep=''),plot = p.all,
         height = 6,width=9,units='cm')