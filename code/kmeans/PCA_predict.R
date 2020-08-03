# relationship between path scores and cognition
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/predict_pcs/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/trainfxns.R')
source('code/misc/processfxns.R')

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
mmse <- cbind(INDDID=INDDIDs,mmse)
# load MoCA data
cog <- read.csv(paste(params$opdir,'processed/MoCA_processed.csv',sep=''),stringsAsFactors = F)[,-1] # remove index column
subj.scores <- lapply(INDDIDs, function(id) cog[cog$INDDID == id,]) # get each subject's scores
#cog <- sapply(subj.scores, function(df) ifelse(length(df)>0,yes=df[which.max(df$AgeatTest),grep('Total',colnames(cog))],no=NaN)) # either NaN if no scores or take most recent
cog <- lapply(subj.scores, function(df) df[which.max(df$AgeatTest),grep('Total',colnames(cog))])
empty.df <- cog[[which(sapply(cog,nrow)==0)[1]]]
empty.df[1,] <- NA
cog[sapply(cog,nrow)==0] <- list(empty.df)
cog <- cbind(INDDID=INDDIDs,do.call('rbind',cog))

# impute scores as average if total present but subscores missing
moca.key <- c(VisuospatialTotal=5,NamingNamingTotal=3,AttentionDigitTotal=2,LanguageRepeatTotal=3,AbstractionTotal=2,RecallTotal=5,OrientationTotal=6)
na.subscores <- !is.na(cog$MoCATotal) & is.na(rowSums(cog[,grep('Total',colnames(cog))]))
for(j in which(na.subscores)){ 
  cog[j,which(is.na(cog[j,]))] <- moca.key[names(cog[j,is.na(cog[j,]),drop=F])]*cog[j,'MoCATotal']/30 # proportional to total score
}

sum(!is.na(cog$MoCATotal))

# load gene data
GOIs <- c('APOE','MAPTHaplotype')#,'C9orf72','LRRK2')
Alleles <- list(APOE=c('E3','E2','E4'),MAPTHaplotype=c('H2','H1'))
Genes <- read.csv(paste(params$opdir,'processed/Genetics_processed.csv',sep=''),stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% INDDIDs,]
Genes <- Genes[order(Genes$INDDID),]

list[Allele.Tables,Genes] <- process.Gene(Genes,GOIs,Alleles,remove.wt=TRUE)

names(Allele.Tables) <- NULL # this is to prevent R from adding back in list item names to df columns

# add IDs back and concatenate across genes
Allele.Tables <- cbind(data.frame(INDDID=Genes$INDDID),do.call('cbind',Allele.Tables))

# load CSF data
CSF.name <- 'Luminex'
CSF.vars <- get.CSF.vars(CSF.name) # get column names corresponding to a-beta, phosphotau and total tau for given type of CSF assay
CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% INDDIDs,]
CSF.sample <- extract.CSF.sample(CSF,CSF.name,n.sample='last')

# load pathology pca
load(file = paste0(params$opdir,'micekmeans/PCA/PolychoricPCA.RData')) # load EFA scores
rownames(scores) <- as.character(INDDIDs.all)
scores.pca <- scores[as.character(INDDIDs),] # select 895 subjects from 901
scores.pca <- cbind(INDDID=INDDIDs,scores.pca)
save(loadings,scores,explained, file=paste0(params$sourcedata.dir,'FigS11a-f_SourceData_PCALoadingsScores.RData'))
# define two dataframes
df.nomoca <- merge(x = mmse,y=Allele.Tables,by.x='INDDID')
df.nomoca <- merge(x=df.nomoca,y=CSF.sample[,c('INDDID',CSF.vars)],by.x='INDDID')
df.nomoca <- merge(df.nomoca,scores.pca,by.x='INDDID')

df.moca <- merge(x = mmse,y=Allele.Tables,by.x='INDDID')
df.moca <- merge(x=df.moca,y=CSF.sample[,c('INDDID',CSF.vars)],by.x='INDDID')
df.moca <- merge(x=df.moca,y=cog,by.x = 'INDDID')
df.moca <- merge(df.moca,scores.pca,by.x='INDDID')

independent.vars <- list(NoMoCA=df.nomoca,MoCA=df.moca)

# predict a bunch of stuff
k.folds <- 5; nreps <- 10
p.meth <- 'rf'
ctrl <- trainControl(method = "repeatedcv", number = k.folds, repeats=nreps,returnResamp='all',classProbs=TRUE,savePrediction='all')
for(x.name in names(independent.vars)){
  df <- independent.vars[[x.name]]
  m.disc <- m.cont <- list()
  for(pc.pred in paste0('PC',1:3)){
    pc.remove <- setdiff(paste0('PC',1:ncol(scores.pca)),pc.pred) # remove PCs you are not predicting
    df.fit <- df[!is.na(rowMeans(df)),setdiff(colnames(df),c(pc.remove,'INDDID'))] # remove NAs and INDDID
    colnames(df.fit)[which(colnames(df.fit)==pc.pred)] <- 'y'
    caret_model <- train(y~., data=df.fit,method=p.meth,returnResamp=TRUE,trControl=ctrl) # train model and store .. for rf can add importance=T to get %IncMSE but r=0.91 with IncNodePurity + much comp time
    m.cont[[pc.pred]] <- caret_model
    
    df.fit$y <- factor(as.numeric(df.fit$y > sd(scores.pca[,pc.pred])),levels = c(0,1),labels = make.names(c(0,1))) # greater than 1 SD
    caret_model_clf <- train(y~., data=df.fit,method=p.meth,returnResamp=TRUE,trControl=ctrl,metric = 'ROC') # train model and store .. for rf can add importance=T to get %IncMSE but r=0.91 with IncNodePurity + much comp time
    caret_model_clf$ROC <- roc(response=as.numeric(caret_model_clf$finalModel$y),predictor=as.numeric(caret_model_clf$finalModel$predicted))
    m.disc[[pc.pred]] <- caret_model_clf 
  }
  
  # plot
  ln.col <- "#3F5151"
  plt.col <- "#9B110E"
  
  p.list <- lapply(m.disc, function(X) ggroc(X$ROC,legacy.axes=TRUE,color=ln.col) + 
                     xlab('FPR') + ylab('TPR')+ geom_abline(slope=1,intercept=0,size=0.5,linetype='dashed',alpha= 0.5) +
                     scale_y_continuous(breaks=c(0,0.5,1)) + scale_x_continuous(breaks=c(0,0.5,1))+
                     theme_classic() + annotate(size=2.5,geom='text',x=-Inf,y=Inf,vjust=1,hjust=-1,label=paste0('AUC = ',signif(X$ROC$auc,2)))+
                     theme(text=element_text(size=6)))
  p <- plot_grid(plotlist=p.list,align='hv',nrow=1)
  ggsave(filename = paste(savedir,'ROC_ClassifyPCs',x.name,'.pdf',sep=''),plot = p,
         height = 6,width=18,units='cm')
  save(m.disc,file=paste0(params$sourcedata.dir,'FigS11g-i_SourceData_PCAROC.RData'))
  p.list <- lapply(m.cont, function(X) ggplot(data=data.frame(x=X$finalModel$predicted,y=X$finalModel$y)) + 
                     geom_point(aes(x=x,y=y),stroke=0,alpha=0.8,color=plt.col,size=1) + geom_smooth(aes(x=x,y=y),method='lm',fill=ln.col,col=ln.col)+
                     ggtitle(paste0(expression(R^2),' = ',signif(mean(X$results$Rsquared),2))) + theme_classic()+ xlab('Predicted') + ylab('Actual')+
                     theme(text=element_text(size=6),plot.title = element_text(size=6,hjust=0.5)))
  #p.xyc(x = X$finalModel$predicted,y=X$finalModel$y,ca='grey'))
  p <- plot_grid(plotlist=p.list,align='hv',nrow=1)
  ggsave(filename = paste(savedir,'Regression_PredictPCs',x.name,'.pdf',sep=''),plot = p,
         height = 3,width=9,units='cm',useDingbats=F)
  save(m.cont,file=paste0(params$sourcedata.dir,'FigS11g-i_SourceData_PCAROC.RData'))
}