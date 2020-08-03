rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'micekmeans/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(savedir,'microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]
summary(microSample.imp)
microSample.imp <- complete(microSample.imp,1)
colnames(microSample.imp) <- colnames(microSample)

W <- polychoric(microSample.imp) # compute feature-by-feature polychoric correlation matrix

save(W,file=paste0(savedir,'PolychoricFeatureCorrelations.RData'))