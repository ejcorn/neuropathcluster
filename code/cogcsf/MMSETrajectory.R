rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1)] # Get rid of index column only

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
data.date <- '122219'
mmse <- read.csv(paste0(params$homedir,'data/INDD_MMSE',data.date,'.csv'),stringsAsFactors=F)
demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics',data.date,'.csv',sep=''),stringsAsFactors = F)

mmse$TestDate <- as.Date(mmse$TestDate,format='%m/%d/%Y')
df.all <- merge(patientSample,mmse,by='INDDID')
df.all <- merge(df.all,demo,by='INDDID')
df.all$DeathMinusTest <- df.all$AgeatDeath - df.all$AgeatTest # calculate date of mmse's relative to autopsy
# add cluster assignement to dataframe
df.all$Cluster <- NA
for(INDDID in unique(df.all$INDDID)){
  df.all$Cluster[df.all$INDDID == INDDID] <- partition[INDDIDs == INDDID]
}
df.all$Cluster <- paste('Cluster',df.all$Cluster)
# remove MMSE if score is outside of possible range
df.all <- df.all[df.all$MMSETotal <=30,]
# remove if ages are unreasonable
df.all <- df.all[df.all$AgeatTest < 115,]

p <- ggplot(df.all[df.all$Cluster=='Cluster 2',]) + geom_line(aes(x=DeathMinusTest,y=MMSETotal,alpha=INDDID)) + theme_classic()+
  scale_x_reverse() + guides(alpha=FALSE)+#scale_color_manual(values=getClusterColors(k)) + 
  xlab('Years Prior to Death') + ylab('MMSE Total') + theme(text=element_text(size=8))
p

