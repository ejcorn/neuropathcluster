rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
source('code/misc/fxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]
#microSample <- as.data.frame(scale(microSample,center = T))
if(sum((colSums(is.na(microSample)) == nrow(microSample))) > 0){
  break
}
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
microSample <- complete(microSample.imp,1)

if(params$dist.met == 'spearman'){
	writeMat(paste(params$opdir,'processed/pathDataForClustering.mat',sep=''),X=as.matrix(microSample-.1)) # subtract 0.1 to make R save as double
} else if(grepl('polychor',params$dist.met)){
	dir.create(paste0(params$opdir,'optimcluster'),recursive=TRUE)
  library(psych)
  DisconnectedSubjects <- which(rowSD(microSample)==0)
	subject.cormat <- polychoric(t(microSample[-DisconnectedSubjects,]))
#   X <- as.data.frame(t(microSample[-DisconnectedSubjects,]))
#   for(j in 1:ncol(X)){X[,j] <- as.character(X[,j])}
#   subject.cormat <- polycor::hetcor(data=X,ML=TRUE)
# 	subject.cormat$rho <- subject.cormat$correlations
	save(subject.cormat,DisconnectedSubjects,file = paste0(params$opdir,'processed/SubjectPolychoricMatrix.RData'))
	writeMat(paste0(params$opdir,'optimcluster/subjectCorrMat.mat'),W=subject.cormat$rho,DisconnectedSubjects=as.matrix(DisconnectedSubjects)) # subtract 0.1 to make R save as double
}
