extralab <- 'CSFOnly'
rm(list = setdiff(ls(), c("params","extralab")))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'predictdisease/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/trainfxns.R')
source('code/misc/plottingfxns.R')

#################
### Load data ###
#################

load(file = paste(savedir,'dzpredict_data',extralab,'.RData',sep=''))

library(scatterplot3d) 
clusterColors <- getClusterColors(k)

pdf(paste(savedir,'CSFTestByCluster3D.pdf'),width=3,height = 3)
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

pdf(paste(savedir,'CSFTestByDisease3D.pdf'),width=3,height = 3)
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