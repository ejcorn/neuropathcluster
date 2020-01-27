rm(list = setdiff(ls(), c("params","dz.exc","n.dx","exc.cl")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'genecluster/',sep='')
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

load(file=paste(savedir,'AlleleTablesCluster.RData',sep=''))

#######################
### Exclude disease ###
#######################

Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(Genes$INDDID),],dz.exc,n.dx)
# exclude dz.exc from corresponding cluster with highest representation
# i.e. Patients without Braak-CERAD AD OR patients not in cluster 2
Rem.Mask <- Rem.Mask | !partitionSample %in% paste('Cluster',exc.cl) 

for(g.i in names(Allele.Tables)){
  Allele.Tables[[g.i]] <- Allele.Tables[[g.i]][Rem.Mask,]
}
partitionSample <- partitionSample[Rem.Mask]

# display number of patients left in each cluster
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Genes: Exclude',dz.exc))
print(cluster.counts)

# exclude clusters with < 20 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample,n=9)
# apply mask to relevant variables
list[clusterColors,clusterNames] <- 
  lapply(list(clusterColors,clusterNames), function(X) X.n.exclude(X,cluster.counts,n=9))
k <- length(clusterNames)

for(g.i in names(Allele.Tables)){
  Allele.Tables[[g.i]] <- Allele.Tables[[g.i]][Rem.Cl,]
}
partitionSample <- flexdim.rowmask(partitionSample,Rem.Cl)

print(paste('Exclude',dz.exc,'n =',length(partitionSample)))
# run linear model to measure allele-cluster relationship
results <- betas <- pvals <- deg.freedom <- list()
for(g.i in names(Allele.Tables)){
  A <- as.data.frame(Allele.Tables[[g.i]])
  results[[g.i]] <- list()
  deg.freedom[[g.i]] <- matrix(NA,nrow = length(clusterNames),ncol =length(clusterNames),dimnames=list(clusterNames,clusterNames))
  for(k.1 in clusterNames){
    results[[g.i]][[k.1]] <- list()
    for(k.2 in clusterNames){
      # select
      # make data frame with two partitions at a time, pairwise test diffs
      df.A.p <- A[partitionSample %in% c(k.1,k.2),-1,drop=FALSE]
      # results[[g.i]][[k.1]][[k.2]] contains regression for Gene_i 
      # where coefficients represent log odds that cluster is k.1, not k.2, given # of allele_i
      partition <- as.numeric(partitionSample[partitionSample %in% c(k.1,k.2)] == k.1) 
      m <- summary(glm(partition ~ ., data=df.A.p,family='binomial'))$coef
      rownames(m)[1] <- colnames(A)[1] #name rows by allele, intercept is WT allele
      # if you couldn't estimate a coefficient because no alleles present, then add a row of NAs
      if(nrow(m) < ncol(A)){m <- rbind(m,matrix(NA,nrow=1,ncol=ncol(m),dimnames=list(colnames(A)[!colnames(A) %in% rownames(m)],colnames(m))))}
      results[[g.i]][[k.1]][[k.2]] <- m
      deg.freedom[[g.i]][k.1,k.2] <- length(partition) - nrow(m) # df is n_observations - n_parameters
    }
  }
  betas[[g.i]] <- extract.by.coef(results[[g.i]],clusterNames,'Estimate')  
  pvals[[g.i]] <- extract.by.coef(results[[g.i]],clusterNames,'Pr(>|z|)')
}

######################
### FDR correction ###
######################

# perform FDR adjustment over all alleles for all genes
# first, make lower tri and diagonal of pval matrices all NA b/c you only correct over 
# unique tests of interest

pvals <- lapply(pvals, function(G)
          lapply(G, function(A) A*ifelse(upper.tri(A),1,NA)))
pvals <- list.fdr.correct(pvals)

# replace NAs on lower tri with reflected upper tri values, replace diag with 1's
for(g.i in names(pvals)){
  for(a.i in names(pvals[[g.i]])){
    pvals[[g.i]][[a.i]][!upper.tri(pvals[[g.i]][[a.i]])] <- 0 # remove NAs    
    pvals[[g.i]][[a.i]] <- pvals[[g.i]][[a.i]] + t(pvals[[g.i]][[a.i]]) # reflect p-vals over diagonal
    pvals[[g.i]][[a.i]][as.logical(diag(TRUE,nrow(pvals[[g.i]][[a.i]])))] <- 1 # make diag 1's
  }
}

# plot betas
max.beta <- max(as.vector(unlist(betas)),na.rm = T)
min.beta <- min(as.vector(unlist(betas)),na.rm = T)
plots <- latex <- list()
for(g.i in names(betas)){
  plots[[g.i]] <- latex[[g.i]] <- list()
  for(a.i in names(betas[[g.i]])){
    plots[[g.i]][[a.i]] <- plot.allele.beta.matrix(betas[[g.i]][[a.i]],pvals[[g.i]][[a.i]],g.i,a.i,min.beta,max.beta,clusterColors)
    list[bmat,pmat] <- list(betas[[g.i]][[a.i]],pvals[[g.i]][[a.i]])
    latex[[g.i]][[a.i]] <- matrix(paste0('$\\beta=',signif(bmat,2),'$, $df=',deg.freedom[[g.i]],'$, $p_\\mathrm{FDR}=',signif(pmat,2),'$'),
                                  k,k,dimnames = list(clusterNames,clusterNames))
    if(dz.exc == 'Alzheimer\'s disease'){save(bmat,pmat,file = paste(savedir,'FigS5c-e',g.i,'-',a.i,'_SourceData.RData',sep=''))}
    if(dz.exc == 'PSP'){save(bmat,pmat,file = paste(savedir,'FigS6a-b',g.i,'-',a.i,'_SourceData.RData',sep=''))}
  }
  p.all <- plot_grid(plotlist = plots[[g.i]], align = 'hv',nrow=1)
  ggsave(filename = paste(savedir,g.i,'BetasExclude',dz.exc,'FromC',paste0(exc.cl,collapse=','),'NDx',n.dx,'.pdf',sep=''),plot = p.all,
         height = 2,width=2.5*length(plots[[g.i]]),units='in')
}