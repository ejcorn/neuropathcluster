rm(list = setdiff(ls(), "params"))
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

#################
### Gene data ###
#################

GOIs <- c('APOE','MAPTHaplotype')#,'C9orf72','LRRK2')
Genes <- read.csv(paste(params$opdir,'processed/Genetics_processed.csv',sep=''),stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% INDDIDs,]
missing.mask <- rowSums(Genes[,GOIs] != '') == length(GOIs)
Genes <- Genes[missing.mask,]  # remove missing
Genes.df <- Genes[,GOIs] # select genes of interest

############################
### Process into alleles ###
############################

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
  rownames(Allele.Tables[[g.i]]) <- NULL
}

partitionSample <- partition[INDDIDs %in% unique(Genes$INDDID)]
partitionSample <- sapply(partitionSample, function(i) paste('Cluster',i))
clusterNames <- sort(unique(partitionSample))
save(Allele.Tables,partitionSample,clusterNames,Genes,file=paste(savedir,'AlleleTablesCluster.RData',sep=''))
