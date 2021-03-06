# Difference between this script and pd_prepdata.R:
# 1. for traditional diagnoses, this script prepares a dataframe that
# includes traditional diagnoses of any probability (possible-intermediate-definite)
# found in NPDx1-4, and (like pd_prepdata.R) uses Braak CERAD diagnosis for AD
# 2. also, if using CSFOnly, this script will add CSF data from about 200 individuals designated as
# putativelly "normal" but without sufficient pathology data to confirm

rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/processfxns.R')
source('code/misc/trainfxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSampleABC.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Prep genetics ###
#####################

GOIs <- c('APOE','MAPTHaplotype')#,'C9orf72','LRRK2')
Alleles <- list(APOE=c('E3','E2','E4'),MAPTHaplotype=c('H2','H1'))
Genes <- read.csv(paste(params$opdir,'processed/Genetics_processed.csv',sep=''),stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% INDDIDs,]
Genes <- Genes[order(Genes$INDDID),]

list[Allele.Tables,Genes] <- process.Gene(Genes,GOIs,Alleles,remove.wt=TRUE)

names(Allele.Tables) <- NULL # this is to prevent R from adding back in list item names to df columns

# add IDs back and concatenate across genes
Allele.Tables <- cbind(data.frame(INDDID=Genes$INDDID),do.call('cbind',Allele.Tables))

################
### Prep CSF ###
################

data.date <- '122219'
CSF.name <- 'Luminex'
CSF.vars <- get.CSF.vars(CSF.name) # get column names corresponding to a-beta, phosphotau and total tau for given type of CSF assay

CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% INDDIDs,]
CSF.sample <- extract.CSF.sample(CSF,CSF.name,n.sample='first')

#################
### Prep MMSE ###
#################

# get an MMSE that occuirred within 1 year of CSF testing
mmse.window <- 365 # set threshold for window 
mmse <- read.csv(paste0(params$homedir,'data/INDD_MMSE',data.date,'.csv'),stringsAsFactors=F)
mmse$TestDate <- as.Date(mmse$TestDate,format='%m/%d/%Y')
# add in normals to CSF data here so you can get MMSE data for every possible subject
# later, the normals will be removed by the merge function if not used
if(grepl('CSF',extralab)){ # if using CSF data, with or without gene data, base MMSE off of earliest CSF sample
	list[CSF.normal,throwaway] <- addnormal.CSF(CSF.sample,CSF.name,data.date)
	CSF.normal <- rbind(CSF.sample,CSF.normal)
	df.mmse <- merge(CSF.normal,mmse,by='INDDID')
	df.mmse$DiffTime <- as.numeric(difftime(df.mmse$TestDate,df.mmse$CSFDate,units='days'))
	# loop through subjects with CSF testing and get all of their MMSEs within 1 year of CSF sample
	df.mmse <- lapply(CSF.normal$INDDID, function(ID) df.mmse[df.mmse$INDDID == ID & abs(df.mmse$DiffTime)<mmse.window,])
	# loop through subjects and extract the single MMSE closest to CSF
	df.mmse <- do.call('rbind',lapply(df.mmse, function(X) X[which.min(abs(X$DiffTime)),]))
} else if(grepl('GeneOnly',extralab)){ # if using gene data only, just take earliest MMSE
	list[Alleles.normal,throwaway] <- addnormal.Genes(Genes,GOIs,Alleles,data.date)
	Alleles.normal <- rbind(Allele.Tables,Alleles.normal)
	df.mmse <- merge(Alleles.normal,mmse,by='INDDID')
	# loop through subjects with genetics and mmse and extract mmse data
	df.mmse <- lapply(unique(df.mmse$INDDID),function(ID) df.mmse[df.mmse$INDDID == ID,]) 
	# loop through subject data and extract single earliest MMSE 
	df.mmse <- do.call('rbind',lapply(df.mmse, function(X) X[which.min(X$TestDate),,drop=FALSE]))
}
df.mmse <- df.mmse[,c('INDDID','MMSETotal')] # only save INDDID and MMSE score. can also save the difference in dates

##########################################
### Merge diagnosis and disease labels ###
##########################################

nl.INDDIDs <- c() # set nl.INDDIDs to empty vector unless extralab contains
if(grepl('CSFGene',extralab)){
	# merge gene and CSF data
	finalINDDIDs <- merge(data.frame(INDDID=CSF.sample$INDDID),data.frame(INDDID=Allele.Tables$INDDID),by.x='INDDID')	 
	CSF.sample <- CSF.sample[CSF.sample$INDDID %in% as.numeric(unlist(finalINDDIDs)),]
	Allele.Tables <- Allele.Tables[Allele.Tables$INDDID %in% as.numeric(unlist(finalINDDIDs)),]
	#
	if(!identical(Allele.Tables$INDDID,CSF.sample$INDDID)){
		print('IDs are misaligned!!!')
		break
	}
	df <- merge(CSF.sample,Allele.Tables,by.x='INDDID')
	if(grepl('AddNormal',extralab)){
		list[CSF.normal,nl.INDDIDs.CSF] <- addnormal.CSF(df,CSF.name,data.date)
		list[Genes.normal,nl.INDDIDs.Gene] <- addnormal.Genes(df,GOIs,Alleles,data.date)
		df.nl <- merge(CSF.normal,Genes.normal,by='INDDID')
		nl.INDDIDs <- df.nl$INDDID
		df <- rbind(df,df.nl)
	}
	df <- df[,c('INDDID',CSF.vars,colnames(Allele.Tables)[-1])]  # colnames(Allele.Tables)[-1] gets the names of non wild type alleles
} else if(grepl('CSFOnly',extralab)){	
	df <- CSF.sample
	# add normals if CSF only
	if(grepl('AddNormal',extralab)){
		list[CSF.normal,nl.INDDIDs] <- addnormal.CSF(df,CSF.name,data.date)
		df <- rbind(df,CSF.normal)

	}
	df <- df[,c('INDDID',CSF.vars)]

} else if(grepl('GeneOnly',extralab)){
	df <- Allele.Tables
	if(grepl('AddNormal',extralab)){
		list[Genes.normal,nl.INDDIDs] <- addnormal.Genes(df,GOIs,Alleles,data.date)
		df <- rbind(df,Genes.normal)
	}
}

if(grepl('MMSE',extralab)){
	df <- merge(df,df.mmse,by='INDDID')
	nl.INDDIDs <- nl.INDDIDs[nl.INDDIDs %in% df$INDDID]
}

###########################
### Prep disease labels ###
###########################
rownames(patientSample) <- INDDIDs # give dx data rownames by INDDID of patient to ensure reliable indexing in correct order
patientSubsample <- patientSample[as.character(df$INDDID[!df$INDDID %in% nl.INDDIDs]),] # store with original names in all NPDx's
list[patientSample,dz.short]<- other.dz(patientSample)

# start with patientSample NPDx1 which contains all the major diagnoses and is aligned with short name functions
dx1 <- patientSample[,'NPDx1',drop=F]
dx.orig <- dummy.data.frame(dx1,names = 'NPDx1')

# replace each disease's binary column with an indicator of whether
# patients have any intermediate-high probability diagnosis of that disease
# as opposed to the default, which is just whether diagnosis is in NPDx1
dzs <- unique(as.character(patientSample$NPDx1))
n.dx <- 4
for(dz in dzs){
	# generate disease labels based on ABC stage for AD and NPDx1-4 for all other diseases
	# this will allow subjects to have multiple diagnoses, as is the case in traditional disease model
  True.dz.mask <- !exclude.dz(patientSample,dz.exc=dz,n.dx=n.dx)
	dx.orig[,paste0('NPDx1',dz)] <- as.numeric(True.dz.mask)
}
colnames(dx.orig) <- gsub('NPDx1','',dz.short)  # format column names to short disease labels
# don't attempt to predict "other" diseases as one category
dx.orig <- dx.orig[,-which(dz.short %in% c('T-O','Oth'))]
dz.short <- dz.short[-which(dz.short %in% c('T-O','Oth'))] 

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
names(partition) <- INDDIDs # name each cluster assignment by INDDID

if(grepl('RandomClusters',extralab)){ # as a negative control, randomly shuffle cluster assignments if indicated by extralab
	clusters <- dummy.data.frame(data.frame(x=sample(as.character(partition),replace = F),row.names=INDDIDs))
} else{
	clusters <- dummy.data.frame(data.frame(x=as.character(partition),row.names=INDDIDs))
}
colnames(clusters) <- sapply(1:k, function(k.i) paste('Cluster',k.i))

if(grepl('AddNormal',extralab)){
	# add observations to the dx matrix that are negative for all diagnoses
	tmp.dx <- as.data.frame(matrix(0,nrow=length(nl.INDDIDs),ncol=ncol(dx.orig),dimnames = list(nl.INDDIDs,colnames(dx.orig))))
	tmp.cl <- as.data.frame(matrix(0,nrow=length(nl.INDDIDs),ncol=ncol(clusters),dimnames = list(nl.INDDIDs,colnames(clusters))))
	# now make a new column in dx corresponding to normal with 1's for added observation
	dx.orig$Normal <- 0
	tmp.dx$Normal <- 1	
	dx.orig <- rbind(dx.orig,tmp.dx)
	clusters <- rbind(clusters,tmp.cl)
	
}

clusters <- clusters[as.character(df$INDDID),]
dx.orig <- dx.orig[as.character(df$INDDID),] # exclude patients without data of interest AND ensure ordering of subjects matches with CSF, etc data
thresh <- 10 # threshold for sample size of positives for each disease
dx <- dx.orig[,colSums(dx.orig)>=thresh] # get rid of poorly represented diseases
clusters <- clusters[,colSums(clusters)>=thresh] # remove any clusters with < 10 observations
partitionSample <- partition[as.character(df$INDDID[!df$INDDID %in% nl.INDDIDs])]

unit.test(identical(rownames(dx),as.character(df$INDDID)),'INDDIDs matched between features and diagnoses','ERROR: INDDIDs not matched')
unit.test(identical(rownames(clusters),as.character(df$INDDID)),'INDDIDs matched between features and cluster assignments','ERROR: INDDIDs not matched')
# exclude cluster 2
#df <- df[clusters[,'Cluster 2'] ==0,]
#clusters <- clusters[clusters[,'Cluster 2'] ==0,]
#clusters <- clusters[,c('Cluster 1','Cluster 3')]

# save this to ultimately exclude some diseases
print(paste(extralab,', n = ',nrow(df),', q	= ',ncol(df)-1,sep=''))
save(df,dx,clusters,patientSubsample,partitionSample,k,nl.INDDIDs,file = paste(savedir,'dzpredict_data',extralab,'.RData',sep=''))

########################################################
### Plot diagnostic composition of prediction sample ###
########################################################

list[patientSubsample,dz.short.subsample]<- other.dz(patientSubsample,NPDx='NPDx1')
if(grepl('RandomClusters',extralab)){partitionSample <- sample(partitionSample,replace=F)} # as a negative control, randomly shuffle cluster assignments if indicated by extralab
list[p.k.dz,df.plt] <- plot.dz.by.cluster(patientSubsample$NPDx1,partitionSample,dz.short.subsample,getClusterColors(k),'Primary\nHistopathologic Diagnosis')

ggsave(filename = paste0(savedir,'PredictionSampleNPDx1Composition',extralab,'.pdf'),plot = p.k.dz,
       height = 5,width=7,units='cm')

