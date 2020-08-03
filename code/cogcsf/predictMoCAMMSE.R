# relationship between path scores and cognition
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/predictcog/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/trainfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1)] # Get rid of index column only
load(paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
microSample.comp <- complete(microSample.imp)
colnames(microSample.comp) <- colnames(microSample)

INDDIDs.all <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample <- remove.Disconnected.Subjects(patientSample,DisconnectedSubjects)
microSample.comp <- remove.Disconnected.Subjects(microSample.comp,DisconnectedSubjects)

if(sum(duplicated(INDDIDs.all))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs.all,DisconnectedSubjects)

# load MMSE data
data.date <- '122219'
mmse <- read.csv(paste0(params$homedir,'data/INDD_MMSE',data.date,'.csv'),stringsAsFactors=F)
subj.scores <- lapply(INDDIDs, function(id) mmse[mmse$INDDID == id,]) # get each subject's scores
mmse <- sapply(subj.scores, function(df) ifelse(length(df)>0,yes=df$MMSETotal[which.max(df$AgeatTest)],no=NaN)) # either NaN if no scores or take most recent
mmse[mmse>30] <- NA # some very large mmse scores in here...get rid of them

# load MoCA data
cog <- read.csv(paste(params$opdir,'processed/MoCA_processed.csv',sep=''),stringsAsFactors = F)[,-1] # remove index column
subj.scores <- lapply(INDDIDs, function(id) cog[cog$INDDID == id,]) # get each subject's scores
#cog <- sapply(subj.scores, function(df) ifelse(length(df)>0,yes=df[which.max(df$AgeatTest),grep('Total',colnames(cog))],no=NaN)) # either NaN if no scores or take most recent
cog <- lapply(subj.scores, function(df) df[which.max(df$AgeatTest),grep('Total',colnames(cog))])
empty.df <- cog[[which(sapply(cog,nrow)==0)[1]]]
empty.df[1,] <- NA
cog[sapply(cog,nrow)==0] <- list(empty.df)
cog <- cbind(INDDID=INDDIDs,do.call('rbind',cog))

###################################
### Predict cognition from path ###
###################################
# predict different measures of cognition using different summaries of pathology

# define dependent variables to iterate through
dependent.vars <- list(VisuospatialTotal=cog[,'VisuospatialTotal',drop=FALSE],OrientationTotal=cog[,'OrientationTotal',drop=FALSE],MMSE=data.frame(MMSE=mmse))

# define independent variables to iterate through
cluster.dummy <- dummy.data.frame(data=data.frame(x=paste('Cluster',partition)))[,-1] # [,-1] because you have an intercept for those with all 0's

patientSample$ADStatus[patientSample$ADStatus==''] <- 'None'
ADNC.dummy <- dummy.data.frame(data=data.frame(x=paste(patientSample$ADStatus)))[,-1]

patientSample$DLBType[patientSample$DLBType %in% c('N/A','')] <-'None'
LBD.dummy <- dummy.data.frame(data=data.frame(x=paste(patientSample$DLBType)))[,-1]

microSample.region <- region.item.average(microSample,fxn='nanmean',region.item = 'region') # regional pathology averages
microSample.item <- region.item.average(microSample,fxn='nanmean',region.item = 'item') # pathology type averages

load(file = paste0(params$opdir,'micekmeans/FactorAnalysisNamed.RData')) # load EFA scores
rownames(scores) <- as.character(INDDIDs.all)
scores.efa <- scores[as.character(INDDIDs),] # select 895 subjects from 901

load(file = paste0(params$opdir,'micekmeans/PolychoricPCA.RData')) # load EFA scores
rownames(scores) <- as.character(INDDIDs.all)
scores.pca <- scores[as.character(INDDIDs),] # select 895 subjects from 901

independent.vars <-list(AllPathology=microSample.comp,Clusters=cluster.dummy,RegionalPathology=microSample.region,PathologyType=microSample.item,ADNC=ADNC.dummy,LBD=LBD.dummy,EFA=scores.efa[,1:k],PCA=scores.pca[,1:k])

y.names <- names(dependent.vars)
x.names <- names(independent.vars)

# set up caret so i can use rf, lasso, glm, and then ideally separately use univariate polychoric correlation
k.folds <- 5; nreps <- 50
ctrl <- trainControl(method = "repeatedcv", number = k.folds, repeats=nreps,returnResamp='all',classProbs=TRUE,savePrediction='all')
for(p.meth in c('cor')){#c('rf','glmnet','cor')){
  # define caret prediction method - enet or leapSeq gives decent results for a linear model, regular lm sucks b/c so many noisy predictors. RF is the best
  # fit model using all of pathology data
  m.all <- list()
  #load(file = paste0(savedir,'PredictCogFromPath_Model_',p.meth,'.RData'))
  
  for(x.name in x.names){
    m.all[[x.name]] <- list()
    for(y.name in y.names){
      print(paste0(x.name,'-',y.name))
      df <- as.data.frame(cbind(independent.vars[[x.name]],y=dependent.vars[[y.name]][,y.name])) # join imputed pathology with given dependent cognitive variable
      #df <- df[!is.na(df$y),]
      df.fit <- df[rowSums(is.na(df))==0,] # remove missing subjects
      if(p.meth != 'cor'){ # for predictive modeling
        m.all[[x.name]][[y.name]] <- train(y~., data=df.fit,method=p.meth,returnResamp=TRUE,trControl=ctrl) # train model and store .. for rf can add importance=T to get %IncMSE but r=0.91 with IncNodePurity + much comp time
      } else if(p.meth == 'cor'){ # for mass univariate  correlation. do partial corr
        r.tests <- corr.test(df.fit[,'y'],df.fit[,setdiff(colnames(df.fit),'y')],method='spearman',adjust = 'fdr') # correlate y with each feature
        m.all[[x.name]][[y.name]]$r <- r.tests$r
        m.all[[x.name]][[y.name]]$p <- r.tests$p # adjusted
      }
    }
  }
  #x.names <- names(independent.vars)
  save(m.all,ctrl,p.meth,x.names,y.names,file = paste0(savedir,'PredictCogFromPath_Model_',p.meth,'.RData'))
}
