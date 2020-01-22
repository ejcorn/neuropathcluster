#extralab <- 'CSFOnlyMMSE'
rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/CSFSpace/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/processfxns.R')
source('code/misc/trainfxns.R')
source('code/misc/plottingfxns.R')

#################
### Load data ###
#################
CSF.name <- 'Luminex'
load(file = paste(params$resultsdir,'predictdisease/dzpredict_data',extralab,'.RData',sep=''))
# add to partitionSample if normal CSF samples have been added
#partitionSample <- c(partitionSample,rep(k+1,nrow(df)-length(partitionSample)))
nl.subjs <- rep(k+1,length(nl.INDDIDs))
names(nl.subjs) <- nl.INDDIDs
partitionSample <- c(partitionSample,nl.subjs)
partitionSample <- partitionSample[as.character(df$INDDID)]

library(scatterplot3d) 
clusterColors <- getClusterColors(k)
clusterColors <- c(clusterColors,'grey')

pdf(paste0(savedir,'CSFTestByCluster3D.pdf'),width=3,height = 3)
scatterplot3d(x=df$LuminexAbeta42,y=df$LuminexTTau,z=df$LuminexPTau, pch = 20, color=clusterColors[partitionSample],cex.symbols = 0.2,cex.axis = 0.25,cex.lab = 0.5)
legend("topright", legend = levels(factor(unique(as.character(partitionSample)))),
     col = clusterColors, pch = 20,cex=0.5,pt.cex = 1)
dev.off()

# disease clusters
p <- ggplot() + geom_point(aes(x=df$LuminexPTau,y=df$LuminexAbeta42,color=as.character(partitionSample)),alpha=0.5,stroke=0) + 
	xlab('CSF Phosphorylated Tau') + ylab('CSF Amyloid-Beta') +
	scale_color_manual(values=clusterColors,limits=as.character(1:k),name='') +
	theme_classic() #+ theme(legend.position='none')
p
ggsave(filename = paste(savedir,'CSFTestByCluster.pdf',sep=''),plot = p,
       height = 12,width=12,units='cm')

# plot each disease cluster vs. all other points to illustrate one vs. all task
for(k.i in 1:k){
	partitionSample.tmp <- partitionSample
	partitionSample.tmp[partitionSample !=k.i] <- 0
	p <- ggplot() + geom_point(aes(x=df$LuminexPTau,y=df$LuminexAbeta42,color=as.character(partitionSample.tmp)),alpha=0.5,stroke=0) + 
		xlab('CSF Phosphorylated Tau') + ylab('CSF Amyloid-Beta') +
		scale_color_manual(values=c(clusterColors[k.i],'grey'),limits=c(as.character(k.i),'0'),label=c(as.character(k.i),'Other'),name='') +
		theme_classic() #+ theme(legend.position='none')
	p
	ggsave(filename = paste(savedir,'CSFTestCluster',k.i,'VsAll.pdf',sep=''),plot = p,
	       height = 12,width=12,units='cm')
}
# disease clusters

dx$Other <- as.numeric(rowSums(dx) == 0) # if patient has diagnosis not highly represented, call other
dx.Dual <- rowSums(dx) == 2 # if patient has two major diagnoses, label separately
dx[dx.Dual,] <- 0 # first set dual diagnoses to nothing
dx$Dual <- as.numeric(dx.Dual) # then set to dual dx logical vector

dz.names <- names(dx) # get all disease names
dz.labels <- sapply(1:nrow(dx), function(i) dz.names[which(dx[i,] == 1)]) # get disease assignments for each patient

pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.names)) # assign a color to each disease
names(dz.colors) <- dz.names

pdf(paste0(savedir,'CSFTestByDisease3D.pdf'),width=3,height = 3)
scatterplot3d(x=df$LuminexAbeta42,y=df$LuminexTTau,z=df$LuminexPTau, pch = 20, color=dz.colors[dz.labels],cex.symbols = 0.2,cex.axis = 0.25,cex.lab = 0.5)
legend("topright", legend = dz.names,
     col = dz.colors, pch = 20,cex=0.5,pt.cex = 1)
dev.off()

p <- ggplot() + geom_point(aes(x=df$LuminexPTau,y=df$LuminexAbeta42,color=dz.labels),alpha=0.5,stroke=0) + 
	xlab('CSF Phosphorylated Tau') + ylab('CSF Amyloid-Beta') +
	scale_color_manual(values=dz.colors,limits=unique(dz.labels),name='') +
	theme_classic()
p
ggsave(filename = paste(savedir,'CSFTestByDisease.pdf',sep=''),plot = p,
       height = 12,width=12,units='cm')


load(file = paste(params$resultsdir,'predictdisease/dzpredict_data',extralab,'.RData',sep=''))
for(dz in colnames(dx)){
	p <- ggplot() + geom_point(aes(x=df$LuminexPTau,y=df$LuminexAbeta42,color=as.character(dx[,dz])),alpha=0.5,stroke=0) + 
		xlab('CSF Phosphorylated Tau') + ylab('CSF Amyloid-Beta') +
		scale_color_manual(values=unname(c('grey','blue')),limits=c('0','1'),labels=c('Other',dz),name='') +
		theme_classic()
	p
	ggsave(filename = paste(savedir,'CSFTest',dz,'VsAll.pdf',sep=''),plot = p,
	       height = 12,width=12,units='cm')
}


## repeat the one vs. all plots with first 2 PCs of CSF data (97% of variance)
load(file = paste(params$resultsdir,'predictdisease/dzpredict_data',extralab,'.RData',sep=''))

# add to partitionSample if normal CSF samples have been added
partitionSample <- c(partitionSample,rep(k+1,nrow(df)-length(partitionSample)))

library(scatterplot3d) 
clusterColors <- getClusterColors(k)
clusterColors <- c(clusterColors,'grey')

CSF.vars <- get.CSF.vars(CSF.name)
if(grepl('MMSE',extralab)){CSF.vars <- c(CSF.vars,'MMSETotal')}
pc <- princomp(df[,CSF.vars],scores=TRUE)

CSF.pcs <- as.data.frame(pc$scores[,1:2])
for(dz in colnames(dx)){
	p <- ggplot() + geom_point(aes(x=CSF.pcs$Comp.1,y=CSF.pcs$Comp.2,color=as.character(dx[,dz])),alpha=0.5,stroke=0) + 
		xlab('CSF PC1') + ylab('CSF PC2') +
		scale_color_manual(values=unname(c('grey','blue')),limits=c('0','1'),labels=c('Other',dz),name='') +
		theme_classic()
	p
	ggsave(filename = paste(savedir,'CSFPC1-2Test',dz,'VsAll.pdf',sep=''),plot = p,
	       height = 12,width=12,units='cm')
}

for(k.i in 1:k){
	p <- ggplot() + geom_point(aes(x=CSF.pcs$Comp.1,y=CSF.pcs$Comp.2,color=as.character(clusters[,k.i])),alpha=0.5,stroke=0) + 
		xlab('CSF PC1') + ylab('CSF PC2') +
		scale_color_manual(values=c(clusterColors[k.i],'grey'),limits=c('1','0'),label=c(as.character(k.i),'Other'),name='') +
		theme_classic() #+ theme(legend.position='none')
	p
	ggsave(filename = paste(savedir,'CSFPC1-2Cluster',k.i,'VsAll.pdf',sep=''),plot = p,
	       height = 12,width=12,units='cm')
}