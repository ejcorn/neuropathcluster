rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

#################
### Load data ###
#################

CSF.name <- 'Luminex'
CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% INDDIDs,]
CSF.vars <- c('LuminexTTau','LuminexPTau','LuminexAbeta42')

CSF.by.pt <- lapply(sort(unique(CSF$INDDID)), function(id) # order tests by date w/in subjs
  CSF[CSF$INDDID==id,c('CSFDate',CSF.vars)])
CSF.by.pt <- lapply(CSF.by.pt, function(X) X[order(X$CSFDate),-1])
#CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) colMeans(X)))
CSF.mean <- do.call('rbind',lapply(CSF.by.pt, function(X) X[nrow(X),])) # just use most recent sample
# for some patients, could compute a feature based on change in CSF proteins
#CSF.diff <- do.call('rbind',lapply(CSF.by.pt, function(X) X[nrow(X),] - X[1,]))

CSF.mean$LuminexPTauAbetaRatio <- CSF.mean$LuminexPTau / CSF.mean$LuminexAbeta42
CSF.mean$LuminexTTauAbetaRatio <- CSF.mean$LuminexTTau / CSF.mean$LuminexAbeta42
CSF.mean <- as.data.frame(scale(CSF.mean,center=T))

load(file = paste(savedir,'dzpredict_dataCSFOnly.RData',sep=''))

CSF.ROC <- list()

for(dz in names(dx)){
	CSF.ROC[[dz]] <- list()
	for(CSF.feat in names(CSF.mean)){
		CSF.ROC[[dz]][[CSF.feat]] <- list()
		CSF.ROC[[dz]][[CSF.feat]]$ROC <- roc(dx[,dz],CSF.mean[,CSF.feat])
		# confusion matrix for thresholding CSF feature at mean value
		m <- confusionMatrix(data=as.factor(as.numeric(CSF.mean[,CSF.feat]>0)),reference=as.factor(dx[,dz]))
		CSF.ROC[[dz]][[CSF.feat]]$Sensitivity<- m$byClass['Sensitivity']
		CSF.ROC[[dz]][[CSF.feat]]$Specificity<- m$byClass['Specificity']
	}
}

p.dz.roc <- list()
for(dz in names(dx)){
	p.dz.roc[[dz]] <- ggroc(lapply(CSF.ROC[[dz]],function(X) X$ROC),legacy.axes=T) + ggtitle(dz) +
		theme(plot.title= element_text(size=8,hjust=0.5),text=element_text(size=8),
			legend.position='none') 
	print(paste(names(CSF.ROC[[dz]]),'AUC =',sapply(CSF.ROC[[dz]], function(X) signif(auc(X$ROC),2))))
}

p.all <- plot_grid(plotlist = p.dz.roc,align = 'hv',nrow=1,axis = 'b')
ggsave(filename = paste(savedir,'CSFSingleFeatureDisease.pdf',sep=''),plot = p.all,
       height = 5,width=18,units='cm')


for(cl in names(clusters)){
	CSF.ROC[[cl]] <- list()
	for(CSF.feat in names(CSF.mean)){
		CSF.ROC[[cl]][[CSF.feat]] <- list()
		CSF.ROC[[cl]][[CSF.feat]]$ROC <- roc(clusters[,cl],CSF.mean[,CSF.feat])
		# confusion matrix for thresholding CSF feature at mean value
		m <- confusionMatrix(data=as.factor(as.numeric(CSF.mean[,CSF.feat]>0)),reference=as.factor(clusters[,cl]))
		CSF.ROC[[cl]][[CSF.feat]]$Sensitivity<- m$byClass['Sensitivity']
		CSF.ROC[[cl]][[CSF.feat]]$Specificity<- m$byClass['Specificity']
	}
}

p.dz.roc <- list()
for(cl in names(clusters)){
	p.dz.roc[[cl]] <- ggroc(lapply(CSF.ROC[[cl]],function(X) X$ROC),legacy.axes=T) + ggtitle(cl) +
		theme(plot.title= element_text(size=8,hjust=0.5),text=element_text(size=8),
			legend.position='none') 
	print(cl)
	print(paste(names(CSF.ROC[[cl]]),'AUC =',sapply(CSF.ROC[[cl]], function(X) signif(auc(X$ROC),2))))
}

p.all <- plot_grid(plotlist = p.dz.roc,align = 'hv',nrow=1,axis = 'b')
ggsave(filename = paste(savedir,'CSFSingleFeatureCluster.pdf',sep=''),plot = p.all,
       height = 5,width=18,units='cm')



idx <- lapply(CSF.ROC, function(CSF.ROC.cl) 
		lapply(CSF.ROC.cl,function(X) which.max(X$ROC$sensitivities+X$ROC$specificities)))
lapply(names(CSF.ROC), function(cl)
	lapply(names(CSF.ROC[[cl]]), function(CSF.feat) CSF.ROC[[cl]][[CSF.feat]]$sensitivities[idx[[cl]][[CSF.mean]]]))

