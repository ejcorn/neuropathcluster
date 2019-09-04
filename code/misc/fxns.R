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

roundUpToMultiple <- function(n,k){
  # round up n to nearest integer multiple of k
  # only works for n > 0
  x = ceil(n/k) * k
  return(x)
}

roundDownToMultiple <- function(n,k){
  # round down n to nearest integer multiple of k
  # only works for n > 0
  x = floor(n/k) * k + 1
  return(x)
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

rowMedians <- function(x){
  return(apply(x,1, function(y) median(y,na.rm = T)))
}

colMedians <- function(x){
  return(apply(x,2, function(y) median(y,na.rm = T)))
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

other.dz <- function(patientSample,NPDx = 'NPDx1'){
  # group together diseases with very low representation
  # make short names for diseases for plotting
  other <- list('Argyrophilic grain disease','Down\'s syndrome','Schizophrenia',
              'Pathological Aging','Tauopathy unclassifiable','FTLD-Other',
              'Cerebrovascular Disease','CAA','PART','Tauopathy, unclassifiable','')
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

exclude.dz.old <- function(patientSample,dz.exc,n.dx){
  # excludes disease based on diagnosis that is > "Low Probability"
  NPDx <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
  NPDx <- lapply(NPDx,function(x) get(x,patientSample))
  NPDx <- lapply(NPDx, function(X) X == dz.exc)  

  if(dz.exc == 'Alzheimer\'s disease'){
    # find patients with either Braak03 <2 or CERAD < 2 -- these folks don't have real AD    
    BraakMask <- patientSample$Braak03 < 2 # patients with low Braak scores
    CERADMask <- patientSample$CERAD < 2 # patients with low CERAD scores
    missingBraak <- is.na(patientSample$Braak03)
    missingCERAD <- is.na(patientSample$CERAD)
    # only look at people with Braak and CERAD scores present
    BraakMask[missingBraak | missingCERAD] <- FALSE
    CERADMask[missingBraak | missingCERAD] <- FALSE
    Rem.Mask <- BraakMask | CERADMask # patients with either low Braak or low CERAD and hence not meeting criteria
    
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

exclude.dz <- function(patientSample,dz.exc,n.dx){
  # for dz.exc == Alzheimer's disease, returns a mask of patients with Braak-CERAD scores both < 1
  # i.e. a mask selecting patients not meeting Braak-CERAD criteria
  # if Braak-CERAD score cannot be determined, we use NPDx Likelihood (subjective pathologist dx as far as I can tell)
  # instead. This occurs for 58/1389 patients
  # for other diseases, this function simply returns a mask selecting patients with low probability or no diagnosis of dz.exc

  NPDx <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
  NPDx <- lapply(NPDx,function(x) get(x,patientSample))
  NPDx <- lapply(NPDx, function(X) X == dz.exc)

  if(dz.exc == 'Alzheimer\'s disease'){
    # find patients with either Braak03 <2 or CERAD < 2 -- these folks don't have real AD    
    BraakMask <- patientSample$Braak03 < 2 # patients with low Braak scores
    CERADMask <- patientSample$CERAD < 2 # patients with low CERAD scores

    # only look at people with Braak and CERAD scores present
    Rem.Mask <- rep(NA,nrow(patientSample)) # initialize mask to select patients without Braak-CERAD AD
    Rem.Mask[BraakMask] <- TRUE # if patients have low Braak scores they can't have Braak-CERAD AD
    Rem.Mask[CERADMask] <- TRUE # if patients have low CERAD scores they can't have Braak-CERAD AD
    Rem.Mask[!BraakMask & !CERADMask] <- FALSE # if patients have high Braak & CERAD scores they have Braak-CERAD AD
    
    # Now all other remaining patients are NAs because they have missing scores
    # This is either patients missing both Braak & CERAD, or with high Braak or CERAD missing the other score
    # such that you cannot rule out Braak-CERAD AD
    # in these 58 patients we will diagnose them with AD based on NPDx1-n.dx having > low probability
    missingBraakCERAD <- is.na(Rem.Mask)

    NPDx.lik <- sapply(1:n.dx, function(i) paste('NPDx',i,'Likelihood',sep=''))
    NPDx.lik <- lapply(NPDx.lik,function(x) get(x,patientSample))
    NPDx.lik <- lapply(NPDx.lik, function(X) X == "Definite (High Probability)" | X == "Probable (Intermediate Probability)")# NOT anyone with non-low probability dx of dz.exc
    # NOT anyone with non-low probability dx of dz.exc...TRUE is only low or no prob of dz.exc. false is > low prob
    NPDx.comb <- mapply(function(NPDx,NPDx.lik){list(!(NPDx & NPDx.lik))},NPDx=NPDx,NPDx.lik=NPDx.lik)

    Lik.Mask <- Reduce('&',NPDx.comb) # mask selecting patients with only low prob or no diagnosis of dz.exc
    Dz.Mask <- Reduce('|',NPDx)  # mask selecting patients with any diagnosis of dz.exc in NPDx1 to n.dx

    # if missing CERAD or Braak, use subjective likelihood. See ExclusionSanityCheck.R for confusion matrix in non-missing people.

    Rem.Mask[missingBraakCERAD] <- Lik.Mask[missingBraakCERAD]
    return(Rem.Mask)

  }
  else{
    # incorporate subjective confidence for non-AD diagnoses
    NPDx.lik <- sapply(1:n.dx, function(i) paste('NPDx',i,'Likelihood',sep=''))
    NPDx.lik <- lapply(NPDx.lik,function(x) get(x,patientSample))
    NPDx.lik <- lapply(NPDx.lik, function(X) X == "Definite (High Probability)" | X == "Probable (Intermediate Probability)")
    # NOT anyone with non-low probability dx of dz.exc
    #NPDx.comb <- mapply(function(NPDx,NPDx.lik){list(!(NPDx & NPDx.lik))},NPDx=NPDx,NPDx.lik=NPDx.lik)
    NPDx.comb <- lapply(NPDx, function(x) !x) # don't incorporate confidence actually... err towards excluding more

    Rem.Mask <- Reduce('&',NPDx.comb)

    return(Rem.Mask)
  }
}

insert.BraakCERAD <- function(patientSample,BraakCERAD){
  # this function takes the logical opposite of the mask generated by exclude.dz, specifically for AD, and 
  # replaces all Alzheimer's mentions in NPDx to correspond only to Braak-CERAD criteria
  # this involves:
  # - removing AD diagnoses when patients don't meet Braak-CERAD criteria
  # - shifting other diagnoses up, i.e. if you delete NPDx1, make NPDx2 the new NPDx1
  # - adding AD diagnoses to NPDx1 if patients meet Braak-CERAD criteria
  # - shifting other diagnoses back, i.e. if you add AD to NPDx1, make old NPDx1 the new NPDx2
  # - move existing AD diagnoses to NPDx1

  n.dx <- 4
  NPDx <- sapply(1:(n.dx+1), function(n) paste('NPDx',n,sep=''))
  dz.swap <- 'Alzheimer\'s disease'
  # allow for new likelihood level "Braak-CERAD"
  #levels(patientSample$NPDx1Likelihood) <- c(levels(patientSample$NPDx1Likelihood),'Braak-CERAD')
  # convert all NPDx columns to characters because I hate factors
  for(NPDx.i in NPDx){
    patientSample[,NPDx.i] <- as.character(patientSample[,NPDx.i])
    patientSample[,paste(NPDx.i,'Likelihood',sep='')] <- as.character(patientSample[,paste(NPDx.i,'Likelihood',sep='')])
  }

  for(i in 1:nrow(patientSample)){ # iterate through each patient
    # find which NPDx houses AD diagnosis
    pt.dx <- sapply(1:n.dx, function(n) as.character(get(NPDx[n],patientSample[i,])))
    pt.lik <- sapply(1:n.dx, function(n) as.character(get(paste(NPDx[n],'Likelihood',sep=''),patientSample[i,])))
    AD.idx <- which(pt.dx == dz.swap) # if empty nothing happens... means pt doesn't have AD
    if(length(AD.idx) > 1){print('AD.idx > 2');break} # two identical diagnoses
    if(!BraakCERAD[i] & !is_empty(AD.idx)){ 
      # remove AD diagnoses when patients don't meet Braak-CERAD criteria      
      # by shifting other diagnoses up
      print(paste('Deleting AD diagnosis from patient',i))
      for(dx in AD.idx:n.dx){
        patientSample[i,NPDx[dx]] <- pt.dx[dx+1]
        patientSample[i,paste(NPDx[dx],'Likelihood',sep='')] <- pt.lik[dx+1]
      }
      new.pt.dx <- sapply(1:n.dx, function(n) as.character(get(NPDx[n],patientSample[i,])))
      new.pt.lik <- sapply(1:n.dx, function(n) as.character(get(paste(NPDx[n],'Likelihood',sep=''),patientSample[i,])))
      print('Old')
      print(pt.dx)
      print(pt.lik)
      print('New')
      print(new.pt.dx)
      print(new.pt.lik)

    } else if(BraakCERAD[i] & is_empty(AD.idx)){
      # add AD diagnoses to NPDx1 when patients meet Braak-CERAD...shift other dx's back
      print(paste('Adding AD diagnosis to patient',i))
      for(dx in n.dx:1){
        # shift NPDx's back, leaving NPDx1 open
        patientSample[i,NPDx[dx+1]] <- pt.dx[dx]
        patientSample[i,paste(NPDx[dx+1],'Likelihood',sep='')] <- pt.lik[dx]
      }
      patientSample[i,'NPDx1'] <- dz.swap # add AD to NPDx1
      patientSample[i,'NPDx1Likelihood'] <- 'Braak-CERAD'

      new.pt.dx <- sapply(1:n.dx, function(n) as.character(get(NPDx[n],patientSample[i,])))
      new.pt.lik <- sapply(1:n.dx, function(n) as.character(get(paste(NPDx[n],'Likelihood',sep=''),patientSample[i,])))
      print('Old')
      print(pt.dx)
      print(pt.lik)
      print('New')
      print(new.pt.dx)
      print(new.pt.lik)

    } else if(BraakCERAD[i] & !is_empty(AD.idx)){
      # move existing AD diagnoses to NPDx1
      print(paste('Moving existing AD diagnosis within patient',i))
      # set NPDx1 to existing diagnoses
      patientSample[i,'NPDx1'] <- dz.swap
      patientSample[i,'NPDx1Likelihood'] <- 'Braak-CERAD'
      pt.dx.old <- pt.dx # save old dx list for testing
      pt.lik.old <- pt.lik
      pt.dx <- pt.dx[-AD.idx] # remove AD from list of existing diagnoses
      pt.lik <- pt.lik[-AD.idx] # remove AD from likelihood list
      for(dx in 1:(n.dx-1)){
        # set NPDx2 through NPDx(n.dx) equal to other diagnoses
        patientSample[i,NPDx[dx+1]] <- pt.dx[dx]
        patientSample[i,paste(NPDx[dx+1],'Likelihood',sep='')] <- pt.lik[dx]
      }
      new.pt.dx <- sapply(1:n.dx, function(n) as.character(get(NPDx[n],patientSample[i,])))
      new.pt.lik <- sapply(1:n.dx, function(n) as.character(get(paste(NPDx[n],'Likelihood',sep=''),patientSample[i,])))
      print('Old')
      print(pt.dx.old)
      print(pt.lik.old)
      print('New')
      print(new.pt.dx)
      print(new.pt.lik)
    } else if(!BraakCERAD[i] & is_empty(AD.idx)){
      # if no diagnosis of AD and doesn't meet Braak-CERAD criteria just print
      print(paste('Not altering diagnoses in patient',i))
    }

  }
  return(patientSample)
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
