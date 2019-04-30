rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')
source('code/misc/plottingfxns.R')
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

##################
### Make plots ###
##################

# traditional disease label
list[patientSample,dz.short]<- other.dz(patientSample)
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.short))

load(file = paste(savedir,'predictcluster_Lassoperf',extralab,'.RData',sep=''))
clusterColors <- getClusterColors(length(cluster.res))
# plot feature weights

weight <- lapply(cluster.res, function(C) 
			do.call('cbind',lapply(C, function(C.i) as.matrix(C.i$FeatureImportance))))
grp <- lapply(names(weight), function(W)
				rep(W,nrow(weight[[W]])-1))
weight <- lapply(weight, function(a) #[-1] removes intercept
	cbind(rowMeans(a,na.rm=T)[-1],
		rowMeans(a,na.rm=T)[-1] + 1.96*rowSD(a)[-1]/sqrt(length(a)),
		rowMeans(a)[-1] - 1.96*rowSD(a)[-1]/sqrt(length(a))))

# compute mean and 95% CI

df.plt <- do.call('rbind',weight)
feat <- rownames(df.plt)
rownames(df.plt) <- 1:nrow(df.plt)
df.plt <- as.data.frame(df.plt)
names(df.plt) <- c('feature.mean','feature.ul','feature.ll')
df.plt$feat <- feat
df.plt$grp <- unlist(grp)
ttl <- ''
df.plt$lab <- paste(signif(df.plt$feature.ll,2),'-',signif(df.plt$feature.ul,2),sep='')

p.dx.dz <- ggplot(data=df.plt) + geom_col(aes(x=grp,y=feature.mean,fill=feat),position=position_dodge()) +
	geom_errorbar(aes(x=grp,ymin=feature.ll,ymax=feature.ul,group=feat),position=position_dodge()) +
	#geom_text(aes(x=dz,y=0.2,label=paste(signif(feature.mean,2),'\n',lab)),size=2,hjust=0.5) +
	scale_y_continuous(expand=c(0,0)) +
	scale_fill_brewer(palette = 'Dark2') +
	xlab('') + ylab('Feature Weight') + ggtitle(ttl) +
	theme_classic() + theme(text = element_text(size=8), plot.title = element_text(hjust=0.5)) +
	theme(legend.position = 'none',axis.text.x = element_text(hjust=0.5,angle=90,color=clusterColors))
p.dx.dz