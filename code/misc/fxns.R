c.names <- function(x,xnames){
  # INPUTS:
  # x: vector
  # xnames: names for each element of x
  #
  # OUTPUT:
  # x with xnames as element names
  #
  # like the c() function but allows names in one line
  
  names(x) <- xnames
  return(x)
}
unit.test <- function(test,pass,fail){
  # test: logical expression to evaluate
  # pass: error to print if true
  # fail: error to print if false
  if(test){
    print(pass)
  } else{print(fail)}

}

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
  if(is.null(dim(x))){x <- x[length(x):1]}
  else if(length(dim(x))>1){x <- x[,ncol(x):1]}
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

CramerV.GOF<-function(X2.output){
  
  # INPUTS:
  # X2.output: chisq.test output list from a goodness of fit (one sample) test
  # 
  # OUTPUT:
  # cramer V statistic: sqrt(X2/(df*n_observations))
  
  return(sqrt(X2.output$statistic/ (sum(X2.output$observed)*(length(X2.output$observed)-1)) ))
}

other.dz <- function(patientSample,NPDx = 'NPDx1'){
  # INPUTS:
  # patientSample: dataframe with NPDx columns
  # NPDx: character name of which NPDx column
  #
  # group together diseases with very low representation into two categories:
  # 'Tau-Other': miscellaneous, unclassified tauopathies
  # 'Other': everything else
  
  # make short names for diseases for plotting

  other.tau <- list('Argyrophilic grain disease','Tauopathy unclassifiable',
              'Tauopathy, unclassifiable','CTE','PART')
  for(O in other.tau){
    patientSample[,NPDx][patientSample[,NPDx] == O] <- 'Tau-Other'
  }
  other <- list('Down\'s syndrome','Schizophrenia',
              'Pathological Aging','FTLD-Other','LATE',
              'Cerebrovascular Disease','CAA','','Hippocampal Sclerosis','CTE')
  for(O in other){
    patientSample[,NPDx][patientSample[,NPDx] == O] <- 'Other'
  }
  
  # index each potential non-Other item in NPDx column with a 'short' name
  long.short.list <- list()
  long.short.list[['Alzheimer\'s disease']] <- 'AD'
  long.short.list[['ADNC - Low']] <- 'lAD'
  long.short.list[['ADNC - Intermediate']] <- 'iAD'
  long.short.list[['ADNC - High']] <- 'hAD'
  long.short.list[['Amyotrophic Lateral Sclerosis']] <- 'ALS'
  long.short.list[['CBD']] <- 'CBD'
  long.short.list[['CTE']] <- 'CTE'
  long.short.list[['FTLD-TDP']] <- 'FTLD'
  long.short.list[['LBD']] <- 'LBD'
  long.short.list[['LBD - Amygdala']] <- 'aLBD'
  long.short.list[['LBD - Brainstem']] <- 'bLBD'
  long.short.list[['LBD - Neocortical']] <- 'nLBD'
  long.short.list[['LBD - Limbic']] <- 'lLBD'
  long.short.list[['LBD - ?']] <- 'LBD'
  long.short.list[['Multiple System Atrophy']] <- 'MSA'
  long.short.list[['Other']] <- 'Oth'
  #long.short.list[['Parkinson\'s disease']] <- 'PD'
  long.short.list[['Pick\'s disease']] <- 'PiD'
  long.short.list[['PART']] <- 'PART'
  long.short.list[['PSP']] <- 'PSP'
  #long.short.list[['FTLD-Other']] <- 'FTLD-Other'
  long.short.list[['Unremarkable adult']] <- 'UA'

  # diseases added for NPDx2
  #long.short.list[['Hippocampal Sclerosis']] <- 'HC Scl.'
  long.short.list[['Lewy body disease, Amygdala-predominant']] <- 'LBD'
  long.short.list[['Lewy body pathology unclassifiable']] <- 'LBD'
  long.short.list[['None']] <- 'None'
  long.short.list[['Tau-Other']] <- 'T-O'

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
              'Cerebrovascular Disease','Dementia of undetermined etiology',
              'Posterior Cortical Atrophy','Progressive muscular atrophy','Other (specify)')
  for(O in other){
    patientSample$ClinicalDx1[patientSample$ClinicalDx1 == O] <- 'Other'
  }
  dz.short <- c('AD','ALS','bvFTD','CBD','DLB','FTLD','MCI','MSA','UA','Oth','PD','PDD','PPA','PSP','Psy.')  
  return(list(ptsampleother=patientSample,dz.short=dz.short))
}

grepl.npdx <- function(patientSample,term,n.dx){
  # INPUTS:
  # patientSample: dataframe of patient information containing NPDx columns
  # term: what to search for (character)
  # n.dx: how many diagnoses to search through
  #
  # OUTPUTS:
  # mask for patients with search term in any of npdx columns
  NPDx <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
  NPDx <- lapply(NPDx,function(x) get(x,patientSample))
  NPDx <- lapply(NPDx,function(X) grepl(term,X))
  return(Reduce('|',NPDx))
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
  # for dz.exc == Alzheimer's disease, returns a mask of patients with AD status of "low" or "none"
  # this is determined by Montine et al. 2012 ADNC grid in table 3 for ABC staging, can also be found in data/
  # i.e. a mask selecting patients who have little to know AD
  # if Braak-CERAD score cannot be determined, we use NPDx Likelihood (subjective pathologist dx as far as I can tell)
  # instead. 
  # for other diseases, this function simply returns a mask selecting patients with low probability or no diagnosis of dz.exc

  NPDx <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
  NPDx <- lapply(NPDx,function(x) get(x,patientSample))
  NPDx <- lapply(NPDx, function(X) X == dz.exc)

  if(dz.exc == 'Alzheimer\'s disease'){

    # only look at people with Braak, CERAD, and ABeta scores present
    Rem.Mask <- rep(NA,nrow(patientSample)) # initialize mask to select patients without Braak-CERAD AD
    Rem.Mask[patientSample$ADStatus == 'Low' | patientSample$ADStatus == 'None'] <- TRUE
    Rem.Mask[patientSample$ADStatus == 'Intermediate'| patientSample$ADStatus == 'High']  <- FALSE # if patients have high Braak & CERAD scores they have Braak-CERAD AD
    
    # Now all other remaining patients are NAs because they have missing scores
    # This is either patients missing both Braak & CERAD, or with high Braak or CERAD missing the other score
    # such that you cannot rule out Braak-CERAD AD
    # in these 58 patients we will diagnose them with AD based on NPDx1-n.dx having > low probability
    missingABC <- is.na(Rem.Mask)

    NPDx.lik <- sapply(1:n.dx, function(i) paste('NPDx',i,'Likelihood',sep=''))
    NPDx.lik <- lapply(NPDx.lik,function(x) get(x,patientSample))
    NPDx.lik <- lapply(NPDx.lik, function(X) X == "Definite (High Probability)" | X == "Probable (Intermediate Probability)")# NOT anyone with non-low probability dx of dz.exc
    # NOT anyone with non-low probability dx of dz.exc...TRUE is only low or no prob of dz.exc. false is > low prob
    NPDx.comb <- mapply(function(NPDx,NPDx.lik){list(!(NPDx & NPDx.lik))},NPDx=NPDx,NPDx.lik=NPDx.lik)

    Lik.Mask <- Reduce('&',NPDx.comb) # mask selecting patients with only low prob or no diagnosis of dz.exc
    Dz.Mask <- Reduce('|',NPDx)  # mask selecting patients with any diagnosis of dz.exc in NPDx1 to n.dx

    # if missing any of ABC criteria, use subjective likelihood. See ExclusionSanityCheck.R for confusion matrix in non-missing people.

    Rem.Mask[missingABC] <- Lik.Mask[missingABC]
    return(Rem.Mask)

  } else{
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

count.ejc <- function(X,items=c()){
  # count number of appearances of item i in vector X
  X <- as.character(X)
  X[is.na(X)] <- 'NA'
  if(length(items)==0){items <- sort(unique(X))}
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

#################################################
### Making labels for plots of feature values ###
#################################################

library(gsubfn)

get.feature.labels <- function(features.of.interest,all.features){
  # INPUTS:
  # features.of.interest: vector of character names of feature categories that are substrings of feature names
  # all.features: vector of all feature names (here, regional pathology scores, colnames(microSample))

  # OUTPUTS:
  # idx: integer vector to order features
  # labels: character vector with axis tick labels for pretty plotting
  features.of.interest <- rev(as.list(features.of.interest))
  idx <- sapply(1:length(features.of.interest), function(i) grep(features.of.interest[[i]], all.features))
  # remove missing feature categories
  features.of.interest <- features.of.interest[lapply(idx,length)>0]
  idx <- idx[lapply(idx,length)>0]
  labels <- list()
  for(i in 1:length(features.of.interest)){
    if(length(idx[[i]]) > 2){ # if multiple features in category, center label and pad with blanks
      labels[[i]] <- c(matrix("",floor(0.5*length(idx[[i]]))),features.of.interest[[i]], c(matrix("",ceiling(0.5*length(idx[[i]])-1))))
    } else if(length(idx[[i]])==1){ # if only 1 feature in category, place category label on that item
      labels[[i]] <- features.of.interest[[i]]
    } else if(length(idx[[i]])==2){ # if only 2 features in category, place category label on 2nd item with left padding
      labels[[i]] <- c("",features.of.interest[[i]])
    }
  }
  idx <- Reduce(c,idx)
  labels <- Reduce(c,labels)
  return(list(idx=idx,labels=labels))
}

micro.order.by <- function(micro,by='type'){
  # INPUTS:
  # micro: subject by feature pathology score matrix
  # by: 'type' or 'region'
  #
  # OUTPUTS:
  # micro.reordered: features ordered by region or by type
  list[pathItems.type,pathRegions.name] <- get.pathscore.names()
  if(by=='type'){
    labels <- pathItems.type
  } else if(by=='region'){
    labels <- pathRegions.name
  }
  indices <- unlist(lapply(labels, function(L) grep(L,colnames(micro))))
  return(micro[,indices])

}

get.pathscore.names <- function(vers = 'original'){
  
  pathItems.prettylab <- list("Neuron Loss","Gliosis","Angiopathy","Ubiquitin","Neuritic Plaques","TDP-43","Tau","Synuclein","Amyloid-beta")
  pathRegions.prettylab <- list("Cing","OC","SM","MF","Ang","CA1/Sub","EC","DG","Amyg","TS","CP","GP","SN","LC","Med","CB","Pons","Mesenc.")
  if(vers == 'original'){ # for use with data as it comes from INDD database in micro csv file
    pathItems.type <- list("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","ThioPlaques","TDP43","Tau","aSyn","Antibody")
    pathRegions.name <- list("Cing","OC","SMT","MF","Ang","CS","EC","DG","Amyg","TS","CP","GP","SN","LC","Med","CB","Pons","MB")
    return(list(pathItems.type=pathItems.type,pathRegions.name=pathRegions.name))
  }else if(vers == 'short'){ # for use with data in all scripts after GenerateSample
    pathItems.type <- list("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody")
    #pathRegions.name <- list("Ci","OC","SM","MF","An","CS","EC","DG","Am","TS","CP","GP","SN","LC","Me","CB","Po","MB")
    pathRegions.name <- list("Cing","OC","SMT","MF","Ang","CS","EC","DG","Amyg","TS","CP","GP","SN","LC","Med","CB","Pons","MB")
    return(list(pathItems.type=pathItems.type,pathRegions.name=pathRegions.name))
  }
  

}

region.by.item.matrix <- function(X,fxn='mean',vers = 'short'){
  # INPUTS:
  # X: matrix whose columns are each pathological feature (region and type)
  # fxn: median or mean
  # vers: version of column names to use, passed through to get.pathscore.names (see get.pathscore.names)
  #
  # OUTPUTS:
  # X.reshape: Column mean/medianss of X, reshaped into a region-by-type matrix
  # excluding any regions or types that are entirely NAs
  
  list[pathItems.type,pathRegions.name] <- get.pathscore.names(vers=vers)
  X.reshape <- matrix(NA,nrow=length(pathItems.type),ncol=length(pathRegions.name),
                                 dimnames = list(pathItems.type,pathRegions.name))
  for(region in pathRegions.name){
    for(item in pathItems.type){
      regionItemMask <- grepl(item,colnames(X)) & grepl(region,colnames(X))
      #print(paste0(region,'-',item,': ',sum(regionItemMask))) # make sure this process only selects one region
      if(fxn=='mean'){X.reshape[item,region] <- mean(X[,regionItemMask])}
      if(fxn=='nanmean'){X.reshape[item,region] <- mean(X[,regionItemMask],na.rm=TRUE)}
      if(fxn=='median'){X.reshape[item,region] <- median(X[,regionItemMask])}
    }
  }

  # remove any regions or types with all NAs
  X.reshape <- X.reshape[rowSums(is.na(X.reshape)) < ncol(X.reshape),colSums(is.na(X.reshape)) < nrow(X.reshape)]
  return(X.reshape)

}

####################################
### Cluster processing functions ###
####################################

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

compute.centroids <- function(X,partition,fxn='mean'){
  # compute centroids from n-by-p data matrix X using partition of length p
  # n is # of observations, p is # of features
  k.idx <- sort(unique(partition))
  if(fxn == 'mean'){
  centroids <- do.call('rbind',lapply(k.idx, function(k.i)
    colMeans(na.rm=T,X[partition==k.i,])))
  } else if(fxn == 'median'){
    centroids <- do.call('rbind',lapply(k.idx, function(k.i)
    colMedians(X[partition==k.i,])))
  }
  rownames(centroids) <- sapply(k.idx, function(i) paste('Cluster',i))
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

order.cluster.by.feature.old <- function(X,centroids,feature.names,max.min=rep('max',nrow(centroids)),fxn='mean'){
  # get indices of clusters based on centroid weight on certain features
  # feature.names: character vector of feature names or list of character vectors in the case of multiple features
  # max.min: whether cluster has max or min amount of path at the specified features
  # centroids: k-by-p matrix of centroids where p is # of features, k is # of clusters
  # X: n-by-p data matrix where n is # of observations 
  
  cluster.init.idx <- lapply(feature.names, function(i) unlist(sapply(i,function(j) grep(j,colnames(X)))))
  if(fxn=='mean'){
    cluster.init.amount <- sapply(cluster.init.idx, function(i) rowMeans(centroids[,i,drop=FALSE],na.rm=T))
  } else if(fxn == 'median'){
    cluster.init.amount <- sapply(cluster.init.idx, function(i) rowMedians(centroids[,i,drop=FALSE]))
  }
  # obtain indices to reorder by max pathology in the groups in cluster.init.names
  
  #cluster.init.reorder <- col.Which.Max(cluster.init.amount)  
  cluster.init.reorder <- rep(NA,nrow(centroids))
  for(k in 1:nrow(centroids)){
    if(max.min[k] == 'max'){
      cluster.init.reorder[k] <- which.max(cluster.init.amount[,k])
    } else if(max.min[k] == 'min'){
      cluster.init.reorder[k] <- which.min(cluster.init.amount[,k])
    }
  }
  return(cluster.init.reorder)
}

thresh.mat <- function(X,op='>',thresh){
  # INPUTS:
  # X: matrix
  # op: '<' or '>' to do < or > thresholding
  # abs_outside retains values whose absolute values are greater than thresh
  # thresh: scalar
  #
  # OUTPUTS:
  # X thresholded at thresh by operation
  # remaining values set to 0
  
  if(op=='<'){return(X*(X<thresh))} 
  else if(op=='>'){return(X*(X>thresh))}
  else if(op=='abs_outside'){return(X*(abs(X)>thresh))}
  else if(op=='abs_inside'){return(X*(abs(X)<thresh))}
  
}

reorder.partition <- function(partition,shuffIdx){
  # assumes partition contains elements 1:k
  # shuffIdx contains the new cluster order in terms of old index, i.e. 4 2 3 1 means cluster 1 = current cluster 4
  # thix fxn reindexes clusters, such that cluster shuffIdx[1] becomes cluster 1
  k <- max(partition)
  if(k == 4 | (length(unique(shuffIdx)) == 6)){ # only works for k=4 and k=6 with uniquely matched clusters
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
  # only check cluster indices that are in the partition (unique(partition))
  k.rng <- sort(unique(partition))
  cl.size <- sapply(k.rng, function(k.i) sum(partition==k.i))
  small.cl.idx <- which(cl.size < thrsh*length(partition)) # linear indices
  small.cl <- k.rng[small.cl.idx] # indices in partition
  if(length(small.cl)==0){
    return(partition) # don't do anything if no small clusters  
  } else if(sum(cl.size[small.cl.idx] >0) == length(small.cl)){ # if all small clsuters are not empty
    for(small.cl.i in small.cl){ # loop through small clusters, assign to nearest other cluster
      small.cl.obs <- which(partition == small.cl.i)
      # compute correlation of small cluster observation with centroids
      for(subj in small.cl.obs){
        s.mat <- cor(t(X[subj,,drop=FALSE]),t(centroids[-small.cl.idx,]),method='spearman',use='pairwise.complete.obs')
        new.cluster <- row.Which.Max(s.mat) # find next best cluster     
        if(length(new.cluster)==1){ # row.Which.Max returns integer(0) if correlations with all centroids can't be computed (might not have due to missing data + low variance in pathology)
          partition[subj] <- new.cluster[]
        } else if(length(new.cluster)==0){partition[subj] <- NA} # store NA if can't assign subject
      }          
    }
    return(partition)
  }

}

disconnect.small.clusters <- function(partition){
  k.rng <- sort(unique(partition))
  cl.size <- sapply(k.rng, function(k.i) sum(partition==k.i))
  small.cl.idx <- which(cl.size < thrsh*length(partition)) # linear indices
  small.cl <- k.rng[small.cl.idx] # indices in partition
  if(length(small.cl)==0){
    return(list(partition=partition,DisconnectedSubjects=NULL)) # don't do anything if the small cluster is actualy just empty
  } else if(sum(cl.size[small.cl.idx] >0) == length(small.cl)){
    DisconnectedSubjects <- which(partition %in% small.cl)
    partition <- partition[-DisconnectedSubjects]
    return(list(partition=partition,DisconnectedSubjects=DisconnectedSubjects))
  }
}

###############################
### Bootstrapping functions ###
###############################

lm.boot <- function(data,indices){
  data <- data[indices,]
  m <- lm(y~.,data=data)
  coefficients(m)[-1]
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
