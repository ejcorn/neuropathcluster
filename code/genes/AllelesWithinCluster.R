rm(list = setdiff(ls(), c("params","dz.exc","n.dx","exc.cl")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'genecluster/',sep='')
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

load(file=paste(savedir,'AlleleTablesCluster.RData',sep=''))

#######################
### Exclude disease ###
#######################

Rem.Mask <- exclude.dz(patientSample[INDDIDs %in% unique(Genes$INDDID),],dz.exc,n.dx)
# exclude dz.exc from corresponding cluster with highest representation
# Find AD within cluster 2
Rem.Mask <- Rem.Mask & partitionSample == paste('Cluster',exc.cl)
# Give them a separate label
partitionSample[Rem.Mask] <- 'tmp'
Rem.Mask <- partitionSample == 'tmp' | partitionSample == paste('Cluster',exc.cl)

# Relabel them as a cluster for temporary implementation ease
partitionSample[partitionSample == 'tmp'] <- paste('Cluster',(1:k)[-exc.cl][1])

for(g.i in names(Allele.Tables)){
  Allele.Tables[[g.i]] <- Allele.Tables[[g.i]][Rem.Mask,]
}

# Isolate all cluster 2 people
partitionSample <- partitionSample[Rem.Mask]

# display number of patients left in each cluster
cluster.counts <- cluster.count(partitionSample,k)
print(paste('Genes: Exclude',dz.exc))
print(cluster.counts)

# exclude clusters with < 2 members after Rem.Mask application
Rem.Cl <- cluster.exclude.mask(cluster.counts,partitionSample,n=20)
# apply mask to relevant variables
list[clusterColors,clusterNames] <- 
  lapply(list(clusterColors,clusterNames), function(X) X.n.exclude(X,cluster.counts))

partitionSample <- flexdim.rowmask(partitionSample,Rem.Cl)

print(paste('Exclude',dz.exc,'n =',length(partitionSample)))
# run linear model to measure allele-cluster relationship
results <- betas <- pvals <- list()
for(g.i in names(Allele.Tables)){
  A <- as.data.frame(Allele.Tables[[g.i]])
  results[[g.i]] <- list()
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
      results[[g.i]][[k.1]][[k.2]] <- m      
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
plots <- list()
for(g.i in names(betas)){
  plots[[g.i]] <- list()
  for(a.i in names(betas[[g.i]])){
    plots[[g.i]][[a.i]] <- plot.allele.beta.matrix(betas[[g.i]][[a.i]],pvals[[g.i]][[a.i]],min.beta,max.beta)
  }
  p.all <- plot_grid(plotlist = plots[[g.i]], align = 'hv',nrow=1)
  ggsave(filename = paste(savedir,g.i,'BetasWithinCluster',dz.exc,'NDx',n.dx,'.pdf',sep=''),plot = p.all,
         height = 2,width=2.5*length(plots[[g.i]]),units='in')
}
