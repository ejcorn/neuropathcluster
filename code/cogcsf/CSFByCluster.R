rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'csfcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

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
print(paste('CSF + Path n =',length(partitionSample)))

CSF.mean <- CSF.mean[,-1]

df <- data.frame(y=as.vector(as.matrix(CSF.mean)),
                 x=rep(colnames(CSF.mean),each=nrow(CSF.mean)),
                 g=rep(partitionSample,ncol(CSF.mean)))
save(df,file = paste(params$sourcedata.dir,'Fig6a_SourceData.RData',sep=''))

# get sample size and effect size

clusterNames <- sort(unique(partitionSample))
results <- samp.size <- p.vals <- latex <- list()
for(CSF.protein in colnames(CSF.mean)){
  results[[CSF.protein]] <- samp.size[[CSF.protein]] <- p.vals[[CSF.protein]] <- latex[[CSF.protein]] <-
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
  p.vals[[CSF.protein]] <- p.vals[[CSF.protein]] * !diag(NA,k) # remove diagonals because they don't correspond to real tests
}
p.vals <- list.fdr.correct(p.vals)
for(CSF.protein in colnames(CSF.mean)){
  for(k1 in clusterNames){
    for(k2 in clusterNames){      
      latex[[CSF.protein]][k1,k2] <- paste0('$M_{',substr(k1,nchar(k1),nchar(k1)),'-',substr(k2,nchar(k2),nchar(k2)),'}=',
                                              round(results[[CSF.protein]][k1,k2],2),'$, $n=',samp.size[[CSF.protein]][k1,k2],'$, $p_','\\','mathrm{FDR}=',signif(p.vals[[CSF.protein]][k1,k2],2),'$')
    }
  }
}
get.mat.inds <- function(mat,inds) sapply(1:ncol(inds), function(j) mat[inds[1,j],inds[2,j]]) # iterate through unique combinations of clusters and get p-values from p-value matrix
cluster.combs <- combn(clusterNames,2)
p.table <- as.data.frame(sapply(p.vals, function(X) get.mat.inds(X,cluster.combs)))
rownames(p.table) <- sapply(1:ncol(cluster.combs), function(j) paste0(cluster.combs[1,j],' vs. ',cluster.combs[2,j]))
#signif(p.table,3)
colnames(p.table) <- c('Total Tau','Phosphorylated Tau','Amyloid-$\beta_{1-42}')
p.table.print <- signif(p.table,2)
p.table.print[p.table.print<0.001] <-'p < 0.001'
xtable(p.table.print,caption = '$p$-values for Figure \\ref{fig:figure6}',label = 'table:figure6_p2pvals')

# now make plots showing only the significant test

source('code/misc/facet_scale_override.R') # load custom functions to allow each facet to have independent y scales. From https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
which.sig.any <- which(rowSums(p.table<0.05)>1)
comps <- lapply(as.data.frame(cluster.combs[,which.sig.any]),function(x) as.character(x))

nearest.k <- function(x,k=25) ceiling(max(x)/k)*k
CSF.vars.sort <- sort(CSF.vars)
p <- ggplot(data=df,aes(x=g,y=y,fill=g)) + geom_boxplot(outlier.size=0.5,size=0.25) + theme_classic()+
  facet_wrap_custom(~x,scales='free',scale_overrides = list(
    scale_override(1,scale_y_continuous(breaks=seq(0,nearest.k(df$y[df$x==CSF.vars.sort[1]]),by=50))),
    scale_override(2,scale_y_continuous(breaks=seq(0,nearest.k(df$y[df$x==CSF.vars.sort[2]]),by=20))),
    scale_override(3,scale_y_continuous(breaks=seq(0,nearest.k(df$y[df$x==CSF.vars.sort[3]]),by=50)))
  ))+ 
  scale_fill_manual(values=clusterColors,name='') + ylab('CSF Protein (pg/ml)') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=2,position=5,comparisons = comps)
# left off: need to either manually delete non-significant comparisons or rewrite this script substantially
p
#p <- ggplot_build(p)
p <- mult.comp.ggpubr(p)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste(savedir,'/CSF',CSF.name,'BoxplotsbyClusterLouvain.pdf',sep=''),height = unit(2.5,'in'),width=unit(3.5,'in'),useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()