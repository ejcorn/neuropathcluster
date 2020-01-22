rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Prep genetics ###
#####################

GOIs <- c('APOE','MAPTHaplotype')#,'C9orf72','LRRK2')
Genes <- read.csv(paste(params$opdir,'processed/Genetics_processed.csv',sep=''),stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% INDDIDs,]
Genes <- Genes[order(Genes$INDDID),]
missing.mask <- rowSums(Genes[,GOIs] != '') == length(GOIs)
Genes <- Genes[missing.mask,]  # remove missing
Genes.df <- Genes[,GOIs] # select genes of interest

# list alleles and each gene of interest
# put wild type allele first
Alleles <- list(APOE=c('E3','E2','E4'),MAPTHaplotype=c('H2','H1'))
Genotypes <- lapply(Genes.df, function(X) unique(X))
# for each gene, get an Allele-by-Genotype table counting # of each allele per genotype
n.alleles.per.genotype <- list()
Allele.Tables <- list()
for(g.i in names(Genotypes)){
  n.alleles.per.genotype[[g.i]] <- 
    sapply(Genotypes[[g.i]], function(gt.i)
      sapply(Alleles[[g.i]], function(A) str_count(pattern = A,string = gt.i)))
  Allele.Tables[[g.i]] <- t(sapply(Genes.df[,g.i], function(gt.i) n.alleles.per.genotype[[g.i]][,gt.i]))
  
  # remove wild type allele [-1] to avoid redundancy, i.e. if you're not H1 you are definitely H2 so 
  # including both as predictors is unnecessary
  Allele.Tables[[g.i]] <- as.data.frame(Allele.Tables[[g.i]][,Alleles[[g.i]][-1]])
  # replenish column names
  names(Allele.Tables[[g.i]]) <- Alleles[[g.i]][-1]
}


names(Allele.Tables) <- NULL # this is to prevent R from adding back in list item names to df columns

# add IDs back and concatenate across genes
Allele.Tables <- cbind(data.frame(INDDID=Genes$INDDID),do.call('cbind',Allele.Tables))

################
### Prep CSF ###
################

CSF.name <- 'Luminex'
CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% INDDIDs,]
CSF.vars <- c('LuminexTTau','LuminexPTau','LuminexAbeta42')

CSF.by.pt <- lapply(sort(unique(CSF$INDDID)), function(id) # order tests by date w/in subjs
  CSF[CSF$INDDID==id,c('CSFDate',CSF.vars)])
CSF.by.pt <- lapply(CSF.by.pt, function(X) X[order(X$CSFDate),-1])
#CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) colMeans(X)))
CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) X[1,])) # just use first available sample
# for some patients, could compute a feature based on change in CSF proteins
#CSF.diff <- do.call('rbind',lapply(CSF.by.pt, function(X) X[nrow(X),] - X[1,]))
CSF.mean <- cbind(data.frame(INDDID=sort(unique(CSF$INDDID))),CSF.mean)

##########################################
### Merge diagnosis and disease labels ###
##########################################

if(extralab == 'CSFGene'){
	# merge gene and CSF data
	finalINDDIDs <- merge(data.frame(INDDID=CSF.mean$INDDID),data.frame(INDDID=Allele.Tables$INDDID),by.x='INDDID')	 
	CSF.mean <- CSF.mean[CSF.mean$INDDID %in% as.numeric(unlist(finalINDDIDs)),]
	CSF.mean <- cbind(data.frame(INDDID=CSF.mean$INDDID),CSF.mean[,-1])
	Allele.Tables <- Allele.Tables[Allele.Tables$INDDID %in% as.numeric(unlist(finalINDDIDs)),]
	#
	if(!identical(Allele.Tables$INDDID,CSF.mean$INDDID)){
		print('IDs are misaligned!!!')
		break
	}
	df <- merge(CSF.mean,Allele.Tables,by.x='INDDID')
} else if(extralab == 'CSFOnly'){	
	df <- CSF.mean
} else if(extralab == 'GeneOnly'){
	df <- Allele.Tables
}

###########################
### Prep disease labels ###
###########################

patientSubsample <- patientSample[INDDIDs %in% df$INDDID,] # store with original names in all NPDx's
list[patientSample,dz.short]<- other.dz(patientSample)

dx1 <- data.frame(NPDx1 = patientSample$NPDx1,stringsAsFactors = F)
dx.orig <- dummy.data.frame(dx1,names = 'NPDx1')

# use "gold standard" AD dx criteria, based on Braak and CERAD
# this allows for some overlap but only with respect to AD. This approach is more conservative
# in that it assumes some prior knowledge of copathologic syndromes outside of clustering framework
True.AD.mask <- !exclude.dz(patientSample,dz.exc='Alzheimer\'s disease',n.dx=4)
dx.orig[,'NPDx1Alzheimer\'s disease'] <- as.numeric(True.AD.mask)

colnames(dx.orig) <- gsub('NPDx1','',dz.short)  # format column names to short disease labels
dx.orig <- dx.orig[INDDIDs %in% df$INDDID,] # exclude patients without data of interest
thresh <- 10 # threshold for sample size for each disease
dx <- dx.orig[,colSums(dx.orig)>=thresh] # get rid of poorly represented diseases

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
partitionSample <- partition[INDDIDs %in% df$INDDID]
clusters <- dummy.data.frame(as.data.frame(as.character(partitionSample)))
colnames(clusters) <- sapply(1:k, function(k.i) paste('Cluster',k.i))
# save this to ultimately exclude some diseases
print(paste(extralab,', n = ',nrow(df),', q	= ',ncol(df)-1,sep=''))
save(df,dx,clusters,patientSubsample,partitionSample,k,file = paste(savedir,'dzpredict_data',extralab,'.RData',sep=''))

########################################################
### Plot diagnostic composition of prediction sample ###
########################################################

list[patientSubsample,dz.short.subsample]<- other.dz(patientSubsample,NPDx='NPDx1')
list[p.k.dz,df.plt] <- plot.dz.by.cluster(patientSubsample$NPDx1,partitionSample,dz.short.subsample,getClusterColors(k),'Primary\nHistopathologic Diagnosis')

ggsave(filename = paste0(savedir,'PredictionSampleNPDx1Composition',extralab,'.pdf'),plot = p.k.dz,
       height = 5,width=7,units='cm')

