# analyze factor scores
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'micekmeans/efa/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

#######################
### factor analysis ###
#######################

microSample.imp <- complete(microSample.imp,1)
load(file = paste0(params$opdir,'micekmeans/FactorAnalysis.RData'))

featureTypes <- list(Syn='Syn',Tau='Tau',TDP43='TDP43',Thio='Thio',
                     NeocorticalAtrophy=c('MFNeuronLoss','MFGliosis'),
                     LimbicAtrophy=c('CSNeuronLoss','CSGliosis'),
                     BasalGangliaAtrophy=c('CPNeuronLoss','CPGliosis'),
                     MedullaAtrophy=c('MedNeuronLoss','MedGliosis'),
                     CerebellumAtrophy=c('CBNeuronLoss','CBGliosis'),
                     MidbrainAtrophy=c('MBNeuronLoss','MBGliosis'),
                     Angiopathy='Angiopathy') # find factors that have highest mean weights on each type of feature
#featureTypes.name <- lapply(featureTypes,function(X) paste(X,collapse = ''))
featureTypes.name <- names(featureTypes)
featureTypes.idx <- lapply(featureTypes,function(types) unlist(lapply(types, function(type) grep(type,rownames(loadings)))))
featureMeans.loadings <- sapply(featureTypes.idx, function(X) colMeans(loadings[X,]))
featureMeans.max <- col.Which.Max(featureMeans.loadings) # find factor with maximum weight on each feature
names(featureMeans.max) <- featureTypes.name
# rename factors that load heavily on single pathological features if unique matches made
if(length(unique(featureMeans.max)) == length(featureMeans.max)){
  colnames(scores)[featureMeans.max] <- names(featureMeans.max)
  colnames(loadings)[featureMeans.max] <- names(featureMeans.max)
}

p <- imagesc(cor(scores),cmap='redblue',clim=c(-1,1)) + coord_equal() +
  theme(axis.text.x  = element_text(angle=90,hjust=1,vjust=0.5),text=element_text(size=6)) +
  theme(legend.key.size = unit(0.1,'cm'))
ggsave(filename = paste0(savedir,'FactorCorrelation.pdf'),plot = p,
       height = 6,width=6,units='cm')

###########################################################################
### Relate factor scores to CSF variables in multiple linear regression ###
###########################################################################

featureTypes.name <- colnames(scores)
# load CSF data and merge with pathology factor scores
CSF.name <- 'Luminex'
CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% INDDIDs,]
CSF.sample <- extract.CSF.sample(CSF,CSF.name,n.sample='first')
CSF.vars <- get.CSF.vars(CSF.name)
df.CSF <- merge(cbind(data.frame(INDDID=INDDIDs),scores[,featureTypes.name]),CSF.sample,by='INDDID')

# predict CSF test value from pathology
df.reg <- lapply(CSF.vars, function(v) cbind(data.frame(y=df.CSF[,v]),df.CSF[,featureTypes.name]))
m.list <- lapply(df.reg, function(df) summary(lm(y~.,data=df)))
names(m.list) <- CSF.vars

# predict pathology from CSF test value

# first regress all factors out of each other to measure independence
df.resid <- lapply(featureTypes.name, function(n) residuals(lm(df.CSF[,n]~.,data=df.CSF[,setdiff(featureTypes.name,n)])))
df.reg <- lapply(df.resid, function(y) cbind(data.frame(y=y),df.CSF[,CSF.vars]))
m.list <- lapply(df.reg, function(df) summary(lm(y~.,data=df)))
names(m.list) <- featureTypes.name

r <- sapply(m.list, function(X) X$r.squared)
p <- ggplot() + geom_col(aes(x=names(r),y=r)) + theme_pubr()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),text=element_text(size=8)) +
  ylab('R^2') + xlab('')
ggsave(filename = paste0(savedir,'FactorsCSFRegress.pdf'),plot = p,
       height = 6,width=6,units='cm')

# repeat but don't regress
df.reg <- lapply(featureTypes.name, function(n) cbind(data.frame(y=df.CSF[,n]),df.CSF[,CSF.vars]))
m.list <- lapply(df.reg, function(df) summary(lm(y~.,data=df)))
names(m.list) <- featureTypes.name

r <- sapply(m.list, function(X) X$r.squared)
p <- ggplot() + geom_col(aes(x=names(r),y=r)) + theme_pubr()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),text=element_text(size=8)) +
  ylab('R^2') + xlab('')
ggsave(filename = paste0(savedir,'FactorsCSF.pdf'),plot = p,
       height = 6,width=6,units='cm')

