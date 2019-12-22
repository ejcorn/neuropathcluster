library(tidyverse)
library(viridis)
library(reshape2)
library(caret)
library(R.matlab)
library(RColorBrewer)
library(lm.beta)
library(brainwaver)
library(dummies)
library(Hmisc)
library(cowplot)
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)

############################
### Cluster by diagnoses ###
############################

clusterColors <- getClusterColors(k)
list[patientSample,dz.short]<- other.clindx(patientSample)
clindx <- sapply(sort(unique(patientSample$ClinicalDx1)),function(X) sum(patientSample$ClinicalDx1==X))
names(clindx) <- sort(unique(patientSample$ClinicalDx1))
clindx


dx <- data.frame(ClinicalDx1 = patientSample$ClinicalDx1,stringsAsFactors = F)
dx <- dummy.data.frame(dx,names = 'ClinicalDx1')
colnames(dx) <- gsub('ClinicalDx1','',colnames(dx))  # format column names
k.dz <- lapply(1:k, function(k.i) colMeans(dx[partition==k.i,]))
k.dz <- do.call('cbind',k.dz)

df.plt <- data.frame(pr = as.vector(k.dz),
                     Dz = rep(rownames(k.dz),k),
                     Dz.short = rep(dz.short,k),
                     cl = as.vector(sapply(1:k, function(k.i) rep(paste('Cluster',k.i),ncol(dx)))))
th <- 0.05 # only label diseases with > 0.5 % of a cluster
df.plt$Dz.short[df.plt$pr < th] <- ''
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(unique(df.plt$Dz)))

p.k.dz <- ggplot(data=df.plt) + geom_col(aes(x=cl,y=pr,fill=Dz),position=position_stack(1)) +
  geom_text(aes(x=cl,y=pr,label=Dz.short,group = Dz),position=position_stack(vjust=.5),size=1.5) +
  scale_y_continuous(expand=c(0,0)) + xlab('') + ylab('Proportion') +
  scale_fill_manual(values = dz.colors,name='Clinical Diagnosis') +
  theme_classic() + theme(text = element_text(size=6),legend.key.size = unit(0.01,'in')) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,color=clusterColors))
p.k.dz
ggsave(filename = paste(savedir,'ClustersByClinicalDiagnosisLouvain.pdf',sep=''),plot = p.k.dz,
       height = 2,width=3,units='in')