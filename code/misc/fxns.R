normalize.ab <- function(x,a,b){
  # scale any vector x between a and b
  x <- (b-a)*(x-min(x))/(max(x)-min(x)) + a
  return(x)
}

imagesc.prep <- function(A){
  melted_cormat <- melt(A)
  melted_cormat$Var1 <- melted_cormat$Var1[length(melted_cormat$Var1):1]
  return(melted_cormat)
}

fliplr <- function(x){
  x <- x[length(x):1]
  return(x)
}

fisher.r.to.z <- function(r){
  r <- 0.5*(log(1+r) - log(1-r))
  return(r)
}

row.Which.Max <- function(x){
  rwm <- unlist(apply(x,1,function(y) which.max(y)))
  return(rwm)
}

rowMax <- function(x){
  rwm <- unlist(apply(x,1,function(y) max(y)))
  return(rwm)
}

col.Which.Max <- function(x){
  cwm <- unlist(apply(x,2,function(y) which.max(y)))
  return(cwm)
}

colMax <- function(x){
  cwm <- unlist(apply(x,2,function(y) max(y)))
  return(cwm)
}

rowMax <- function(x){
  return(apply(x,1,function(y) max(y,na.rm=T)))
}

rowSD <- function(x){
  return(apply(x,1,sd))
}

colSD <- function(x){
  return(apply(x,2,sd))
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

other.dz <- function(patientSample,NPDx = 'NPDx1'){
  # group together diseases with very low representation
  # make short names for diseases for plotting
  other <- list('Argyrophilic grain disease','Down\'s syndrome','Schizophrenia',
              'Pathological Aging','Tauopathy unclassifiable','FTLD-Other',
              'Cerebrovascular Disease','CAA','PART','Tauopathy, unclassifiable')
  for(O in other){
    patientSample[,NPDx][patientSample[,NPDx] == O] <- 'Other'
  }
  
  # index each potential non-Other item in NPDx column with a 'short' name
  long.short.list <- list()
  long.short.list[['Alzheimer\'s disease']] <- 'AD'
  long.short.list[['Amyotrophic Lateral Sclerosis']] <- 'ALS'
  long.short.list[['CBD']] <- 'CBD'
  long.short.list[['FTLD-TDP']] <- 'FTLD'
  long.short.list[['LBD']] <- 'LBD'
  long.short.list[['Multiple System Atrophy']] <- 'MSA'
  long.short.list[['Other']] <- 'Oth'
  long.short.list[['Parkinson\'s disease']] <- 'Pa.D'
  long.short.list[['Pick\'s disease']] <- 'Pi.D'
  long.short.list[['PSP']] <- 'PSP'
  #long.short.list[['FTLD-Other']] <- 'FTLD-Other'
  long.short.list[['Unremarkable adult']] <- 'UA'

  # diseases added for NPDx2
  long.short.list[['Hippocampal Sclerosis']] <- 'HC Scl.'
  long.short.list[['Lewy body disease, Amygdala-predominant']] <- 'LBD'
  long.short.list[['Lewy body pathology unclassifiable']] <- 'LBD'
  long.short.list[['None']] <- 'None'

  # diseases present in NPDx
  dz.pres <- sort(unique(as.character(patientSample[,NPDx])))

  dz.short <- unname(unlist(long.short.list[dz.pres]))
  return(list(ptsampleother=patientSample,dz.short=dz.short))
}

other.clindx <- function(patientSample){
  # group together clinical diagnoses with very low representation
  # make short names for diseases for plotting
  patientSample$ClinicalDx1 <- as.character(patientSample$ClinicalDx1)
  other <- list('','Down\'s syndrome','Hydrocephalus','Motor neuron disease NOS',
              'Generally Impaired','Primary lateral sclerosis',
              'Cerebrovascular Disease',
              'Posterior Cortical Atrophy','Progressive muscular atrophy','Other (specify)')
  for(O in other){
    patientSample$ClinicalDx1[patientSample$ClinicalDx1 == O] <- 'Other'
  }
  dz.short <- c('AD','ALS','bvFTD','CBD','DNOS','FTLD','LBD','MCI','MSA','Norm.','Oth','Pa.D','PPA','PSP','Psy.')  
  return(list(ptsampleother=patientSample,dz.short=dz.short))
}

exclude.dz.basic <- function(patientSample,dz.exc,n.dx){
  # excludes disease based on presence of any diagnosis
  NPDx <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
  NPDx <- lapply(NPDx,function(x) get(x,patientSample))
  NPDx <- lapply(NPDx, function(X) X != dz.exc)
  NPDx.Mask <- Reduce('&',NPDx)
  return(Rem.Mask)
}

exclude.dz <- function(patientSample,dz.exc,n.dx){
  # excludes disease based on diagnosis that is > "Low Probability"
  NPDx <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
  NPDx <- lapply(NPDx,function(x) get(x,patientSample))
  NPDx <- lapply(NPDx, function(X) X == dz.exc)  

  if(dz.exc == 'Alzheimer\'s Disease'){
    # find patients with either Braak03 <2 or CERAD < 2 -- these folks don't have real AD    
    BraakMask <- patientSample$Braak03 < 2
    CERADMask <- patientSample$CERAD < 2
    missingBraak <- is.na(patientSample$Braak03)
    missingCERAD <- is.na(patientSample$CERAD)
    # only look at people with Braak and CERAD scores present
    BraakMask[missingBraak | missingCERAD] <- FALSE
    CERADMask[missingBraak | missingCERAD] <- FALSE
    Rem.Mask <- BraakMask | CERADMask
    
    NPDx.lik <- sapply(1:n.dx, function(i) paste('NPDx',i,'Likelihood',sep=''))
    NPDx.lik <- lapply(NPDx.lik,function(x) get(x,patientSample))
    NPDx.lik <- lapply(NPDx.lik, function(X) X != "Possible (Low Probability)")
    # NOT anyone with non-low probability dx of dz.exc
    NPDx.comb <- mapply(function(NPDx,NPDx.lik){list(!(NPDx & NPDx.lik))},NPDx=NPDx,NPDx.lik=NPDx.lik)

    Lik.Mask <- Reduce('&',NPDx.comb)
    # if missing CERAD or Braak, use subjective likelihood
    Rem.Mask[missingBraak | missingCERAD] <- 
      Lik.Mask[missingBraak | missingCERAD]
    return(Rem.Mask)

  }
  else{
    # incorporate subjective confidence for non-AD diagnoses
    NPDx.lik <- sapply(1:n.dx, function(i) paste('NPDx',i,'Likelihood',sep=''))
    NPDx.lik <- lapply(NPDx.lik,function(x) get(x,patientSample))
    NPDx.lik <- lapply(NPDx.lik, function(X) X != "Possible (Low Probability)")
    # NOT anyone with non-low probability dx of dz.exc
    NPDx.comb <- mapply(function(NPDx,NPDx.lik){list(!(NPDx & NPDx.lik))},NPDx=NPDx,NPDx.lik=NPDx.lik)

    Rem.Mask <- Reduce('&',NPDx.comb)

    return(Rem.Mask)
  }
}

cluster.count <- function(partition,k){
  # Counts number of each cluster in partition, 
  # which is a character vector labeling subjects as "Cluster 1", "Cluster 2", etc
  cluster.counts <- sapply(1:k, function(k.i) sum(partition == paste('Cluster',k.i)))
  names(cluster.counts) <- sapply(1:k, function(i) paste('Cluster',i))
  print(cluster.counts)
  return(cluster.counts)
}

count.ejc <- function(X){
  # count number of appearances of item i in vector X
  items <- sort(unique(X))
  counts <- sapply(items, function(i) sum(X==i))
  names(counts) <- items
  return(counts)
}

cluster.exclude.mask <- function(cluster.counts,partition,n=2){
  # keep clusters with > n subjects
  # cluster.counts has number of subjects in each cluster
  # n is threshold for exclusion 
  # Rem.Cl is a mask of subjects not belonging to small clusters
  Rem.Cl <- sapply(which(cluster.counts > n), function(i) paste('Cluster',i))
  Rem.Cl <- partition %in% Rem.Cl
  return(Rem.Cl)

}

X.n.exclude <- function(X,counts,n=1){
  # function designed to exclude colors for clusters with <= n subject in them
  # after other exclusions have been applied
  return(X[which(counts>n)])
}

flexdim.rowmask <- function(X,mask){
  # apply mask to X flexibly depending on dimensions
  # if X is a vector, mask directly
  # if X is a matrix (usually subjects-by-features), mask rows
  if(is.vector(X)){
    return(X[mask])
  }
  else if(!is.vector(X)){
    return(X[mask,])
  }
}

lm.boot <- function(data,indices){
  data <- data[indices,]
  m <- lm(y~.,data=data)
  coefficients(m)[-1]
}

remove.Disconnected.Subjects <- function(X,DisconnectedSubjects){
  if(!is.null(ncol(X))){
    if(!is_empty(DisconnectedSubjects)){
      X <- X[-DisconnectedSubjects,]
      return(X)
    }
  } else if(is.null(ncol(X))){
    if(!is_empty(DisconnectedSubjects)){
      X <- X[-DisconnectedSubjects]
      return(X)
    }
  }
  return(X)
}

compute.centroids <- function(X,partition){
  # compute centroids from n-by-p data matrix X using partition of length p
  # n is # of observations, p is # of features
  k <- max(partition)
  centroids <- do.call('rbind',lapply(1:k, function(k.i)
    colMeans(na.rm=T,X[partition==k.i,])))
  rownames(centroids) <- sapply(1:k, function(i) paste('Cluster',i))
  return(centroids)
}

order.cluster.by.feature <- function(X,centroids){
  # get indices of clusters based on centroid weight on certain features
  # centroids: k-by-p matrix of centroids where p is # of features, k is # of clusters
  # X: n-by-p data matrix where n is # of observations 
  # this function basically ensures clusters will always be in same order
  # despite subtle variation in centroids
  # reason for sequential marching through feature.names is that sometimes
  # AD/Thio/Abeta cluster has more Tau than the tauopathy cluster so go through it first
  # as written, will only work for k=4 but won't be applied in ProcessCluster if not unique

  feature.names <- c('Thio','Tau','TDP43','Syn') # march through in this order
  feature.order <- c(2,1,3,4) # arrange clusters with sequential max of above features in this order
  k <- length(feature.names)

  cluster.init.idx <- lapply(feature.names, function(i) grep(i,colnames(X)))
  cluster.init.amount <- sapply(cluster.init.idx, function(i) rowMeans(centroids[,i],na.rm=T))
  # obtain indices to reorder by max pathology in the groups in cluster.init.names
  cluster.init.reorder <- col.Which.Max(cluster.init.amount)  
  cluster.init.reorder <- rep(NA,k)
  while(sum(is.na(cluster.init.reorder)) > 0){ # look at all clusters one feature at a time, then ditch cluster once matched
    cluster.init.idx <- lapply(feature.names, function(i) grep(i,colnames(X)))
    cluster.init.amount <- sapply(cluster.init.idx, function(i) rowMeans(centroids[,i],na.rm=T))    
    cluster.match <- which.max(cluster.init.amount[,1]) # match first feature to a cluster
    cluster.init.reorder[cluster.match] <- feature.order[1]  # store that cluster i's predetermined index in position i

    # now get rid matched features and set matched clusters to 0
    feature.names <- feature.names[-1]
    feature.order <- feature.order[-1]
    centroids[cluster.match,] <- 0

  }
  return(cluster.init.reorder)
}

order.cluster.by.feature.old <- function(X,centroids,feature.names){
  # get indices of clusters based on centroid weight on certain features
  # feature.names: character vector of feature names
  # centroids: k-by-p matrix of centroids where p is # of features, k is # of clusters
  # X: n-by-p data matrix where n is # of observations 
  
  cluster.init.idx <- lapply(feature.names, function(i) grep(i,colnames(X)))
  cluster.init.amount <- sapply(cluster.init.idx, function(i) rowMeans(centroids[,i,drop=FALSE],na.rm=T))
  # obtain indices to reorder by max pathology in the groups in cluster.init.names
  cluster.init.reorder <- col.Which.Max(cluster.init.amount)  
  return(cluster.init.reorder)
}

reorder.partition <- function(partition,shuffIdx){
  # assumes partition contains elements 1:k
  # shuffIdx contains the new cluster order in terms of old index, i.e. 4 2 3 1 means cluster 1 = current cluster 4
  # thix fxn reindexes clusters, such that cluster shuffIdx[1] becomes cluster 1
  k <- max(partition)
  if(k == 4){
    newPartition <- rep(NA,length(partition))
    # get indices of member of each of the initial partition
    partition.idx.by.k <- lapply(1:k, function(k.i) which(partition == k.i))
    for(k.i in 1:k){
      # reindex clusters such that cluster shuffIdx[k.i] becomes cluster k.i
      newPartition[partition.idx.by.k[[shuffIdx[k.i]]]] <- k.i
    }
    return(newPartition)
  } else if(k != 4){
    return(partition)
  }
}

remove.Partition.Gap <- function(partition){
  newPartition <- rep(NA,length(partition))
  oldClusters <- sort(unique(partition))
  newClusters <- rank(oldClusters)
  for(i in 1:length(oldClusters)){
    newPartition[partition==oldClusters[i]] <- newClusters[i]
  }
  return(newPartition)
}

annex.small.clusters <- function(partition,X,centroids,thrsh=0.01){
  # deal with a single small cluster
  # annex members to other clusters based on spearman correlation between corresponding data in X and centroids of other clusters
  rownames(X) <- 1:nrow(X) # makes linear row names -- needed for when using subsamples of full data matrix
  k <- max(partition)
  cl.size <- sapply(1:k, function(k.i) sum(partition==k.i))
  small.cl <- which(cl.size < thrsh*length(partition))
  if(is_empty(small.cl)){
    return(partition) # don't do anything if no small clusters
  } else if(sum(cl.size[small.cl] == 0) > 0){ # if any small cluster is empty
    return(partition) # don't do anything if any small cluster is actualy just empty
  } else if(sum(cl.size[small.cl] >0) == length(small.cl)){ # if all small clsuters are not empty
    for(small.cl.i in small.cl){ # loop through small clusters, assign to nearest other cluster
      small.cl.obs <- which(partition == small.cl.i)
      # compute correlation of small cluster observation with centroids
      s.mat <- cor(t(X[small.cl.obs,]),t(centroids),method='spearman',use='pairwise.complete.obs')
      s.mat[,small.cl.i] <- NA # don't allow reassignment to small cluster
      new.clusters <- row.Which.Max(s.mat) # find next best cluster
      for(subj in small.cl.obs){ # reassign small cluster observations to next best cluster
        partition[subj] <- new.clusters[as.character(subj)]
      }
    }
    return(partition)
  }

}

disconnect.small.clusters <- function(partition){
  k <- max(partition)
  cl.size <- sapply(1:k, function(k.i) sum(partition==k.i))
  small.cl <- which(cl.size < thrsh*length(partition))
  if(cl.size[small.cl] == 0){
    return(list(partition=partition,DisconnectedSubjects=NULL)) # don't do anything if the small cluster is actualy just empty
  } else if(cl.size[small.cl] >0){
    DisconnectedSubjects <- which(partition==small.cl)
    partition <- partition[partition!=small.cl]
    return(list(partition=partition,DisconnectedSubjects=DisconnectedSubjects))
  }
}


boot.reg <- function(df,nboot){
  # bootstrap linear regression to test whether coefficients for all predictors
  # are ≠ 0
  # df must have a column named 'y' and all other columns are predictors
  m <- lm(y ~ .,data=df)
  df.boot <- boot(data=df,statistic=lm.boot,R=nboot)

  # compute non-parametric p-value 
  # Ho: beta ≠ 0
  # p-val interp: P(beta < 0) or P(beta > 0)
  # per http://statweb.stanford.edu/~tibs/sta305files/FoxOnBootingRegInR.pdf
  # make it 2-tailed (I do this below), and add 1 to numerator and denom (i don't do this)
  p.vals <- 2*colMeans(df.boot$t < 0)
  # for negative betas, reflect p-value around 0.5
  p.vals[p.vals > 0.5] <- (2*colMeans(df.boot$t > 0))[p.vals > 0.5]
  ul.ll <- t(sapply(1:ncol(df.boot$t), function(i) quantile(x = df.boot$t[,i],c(0.025,0.975))))
  return(list(p.vals=p.vals,ul.ll=ul.ll,m=m))
}

extract.by.coef <- function(results,grp.comparison,value,remove.diag=TRUE){
  # results: list of glm model objects perform pairwise comparisons 
  # of predicting elements of grp.comparison.  
  # grp.comparison: character vector of group names, correspond to results list element names
  # value: 'Pr(>|z|)' for p-value or 'Estimate' for beta weight
  # returns list of matrices containing p-val or beta weight for each coefficient (including intercept)
  val <- list()
  coef.names <- rownames(results[[1]][[1]]) # names of coefficients
  val <- lapply(coef.names, function(a.i) # loop through alleles
  t(sapply(grp.comparison, function(k1) # transpose so that you have a matrix of log odds ratios, where ij is log odds cluster=i, not j
      sapply(grp.comparison, function(k2) # set to 0 if p < FDR threshold
        results[[k1]][[k2]][a.i,value]))))
  names(val) <- coef.names
  # make diagonal = NA b/c it's a singular regression, when getting log odds for grp1 vs grp1 
  # only doing that computation for ease of implementation  

  if(remove.diag){
   for(v in names(val)){
    val[[v]][as.logical(diag(TRUE,nrow(val[[v]])))] <- NA
   }
  }
  
  return(val)
}
