# Plot CSF sample timing relative to disease onset

rm(list = setdiff(ls(), c("params")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/patientcharacteristics/',sep='')
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
list[CSF.normal,throwaway] <- addnormal.CSF(CSF.sample,CSF.name,data.date)
CSF.normal <- rbind(CSF.sample,CSF.normal)
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

clusterNames <- paste('Cluster',1:k)
clusterColors <- getClusterColors(k)
comps <- combn(clusterNames,m = 2,simplify = FALSE)
p <- ggplot(data=CSF.onset,aes(x=as.character(Cluster),y=CSFMinusOnset,fill=as.character(Cluster))) + geom_boxplot(outlier.size=0.5) + theme_classic()+
  scale_y_continuous(breaks=seq(-10,20,by = 10)) + ggtitle('CSF Sample - Onset')+
  scale_fill_manual(values=clusterColors,name='') + ylab('Years') + xlab('') +  
  theme(legend.position = 'none',legend.key.size = unit(0.1,'in')) + #c(0.1,0.75)
  theme(text= element_text(size=8),plot.title=element_text(size=8,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,color = clusterColors)) + 
  stat_compare_means(size=1.75,comparisons = comps,position=5,method = "wilcox.test") 
p <- mult.comp.ggpubr(p,method='none')
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file = paste0(savedir,'CSFSampleDateMinusDiseaseOnsetByCluster.pdf'),height = unit(6/2.54,'in'),width=unit(4.5/2.54,'in'),useDingbats = F)
plot(ggplot_gtable(p)) #+ scale_y_continuous(limits=c(0,ymax))
dev.off()
save(CSF.onset, file = paste0(params$sourcedata.dir,'FigS5d_SourceData_CSFvsDiseaseOnset.RData'))
