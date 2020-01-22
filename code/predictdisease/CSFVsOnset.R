# Plot CSF sample timing relative to disease onset

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
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors=F)[,-(1)] # Get rid of index column and INDDIDs
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
list[CSF.normal,throwaway] <- addnormal(CSF.sample,CSF.name,data.date)
df.mmse <- merge(CSF.normal,mmse,by='INDDID')
df.mmse$DiffTime <- as.numeric(difftime(df.mmse$TestDate,df.mmse$CSFDate,units='days'))
# loop through subjects with CSF testing and get all of their MMSEs within 1 year of CSF sample
df.mmse <- lapply(CSF.normal$INDDID, function(ID) df.mmse[df.mmse$INDDID == ID & abs(df.mmse$DiffTime)<mmse.window,])
# loop through subjects and extract the single MMSE closest to CSF
df.mmse <- do.call('rbind',lapply(df.mmse, function(X) X[which.min(abs(X$DiffTime)),]))

df.mmse <- df.mmse[,c('INDDID','MMSETotal')] # only save INDDID and MMSE score. can also save the difference in dates

########################################################################
### Plot distribution of disease onset dates relative to CSF testing ###
########################################################################

onset <- read.csv(paste0('data/INDD_GlobalOnset',data.date,'.csv'),stringsAsFactors = F)
CSF.onset <- merge(CSF.sample,onset,by='INDDID')
CSF.onset$CSFYear <- year(CSF.onset$CSFDate)
CSF.onset$CSFMinusOnset <- CSF.onset$CSFYear - CSF.onset$GlobalYearOnset

names(partition) <- INDDIDs # name each cluster assignment by INDDID
CSF.onset$Cluster <- paste('Cluster',partition[as.character(CSF.onset$INDDID)])
ggplot(CSF.onset) + geom_boxplot(aes(x=as.character(Cluster),y=CSFMinusOnset),alpha=0.5,color='blue') + theme_classic() +
  theme(axis.text.x=element_text(angle=90)) + xlab('')

CSF.onset <- merge(CSF.onset,patientSample,by='INDDID')

ggplot(CSF.onset) + geom_boxplot(aes(x=as.character(NPDx1),y=CSFMinusOnset),alpha=0.5,color='blue') + theme_classic() +
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) + xlab('')

