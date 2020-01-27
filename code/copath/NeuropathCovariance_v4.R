rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'copath/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]
#microSample <- scale(microSample,center=T)
if(sum((colSums(is.na(microSample)) == nrow(microSample))) > 0){
  break
}

###########################################
##### Neuropath correlation structure #####
###########################################

pathItems.type <- list("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody")
pathItems.prettylab <- list("Neuron Loss","Gliosis","Angiopathy","Ubiquitin","Neuritic Plaques","TDP-43","Tau","Synuclein","Amyloid-beta")
pathItems.index <- sapply(1:length(pathItems.type), function(i) grep(pathItems.type[[i]], colnames(microSample)))
pathItems.labels <- sapply(1:length(pathItems.type), function(i) c(matrix("",floor(0.5*length(pathItems.index[[i]]))),
                                                                   pathItems.prettylab[[i]], c(matrix("",ceiling(0.5*length(pathItems.index[[i]])-1)))))
pathItems.index <- Reduce(c,pathItems.index)
pathItems.labels <- Reduce(c,pathItems.labels)

microSample.typeorder <- microSample[,pathItems.index]


A <- rcorr(as.matrix(microSample.typeorder),type='spearman')
p <- A$P
# remove NAs -- happens when not enough obs present across subjects for feature pair
# or when only remaining obs for a feature pair are constant
# set NAs to 0 and ignore p-value

p.calc.fdr <- na.exclude(p[upper.tri(p)])
p.thrsh <- compute.FDR(p.calc.fdr,q=0.05) # fdr adjust for all unique comparisons
p
p <- p < p.thrsh
A <- A$r*p
A[diag(1,nrow = nrow(A))] <- 1 # make diag 1
#A[is.na(A)] <- 0

# identify pairwise available sample sizes
SS <- sapply(microSample.typeorder, function(i) 
        sapply(microSample.typeorder, function(j) sum(!is.na(i) & !is.na(j))))

poorlySampled <- c('Ubiquitin','OC','DG')
poorlySampled.index <- unlist(sapply(poorlySampled, function(i) grep(i,colnames(SS))))
summary(SS[-poorlySampled.index,-poorlySampled.index][!is.na(A)[-poorlySampled.index,-poorlySampled.index]])
summary(SS[poorlySampled.index,poorlySampled.index][!is.na(A)[poorlySampled.index,poorlySampled.index]])


pathItems.rect.idx <- unlist(lapply(1:length(pathItems.type), function(i) 
  length(grep(pathItems.type[[i]],colnames(microSample.typeorder)))))
pathItems.rect.idx <- fliplr(pathItems.rect.idx)

melted_cormat <- melt(A)
melted_cormat$Var1 <- fliplr(melted_cormat$Var1)
cs.rect <- cumsum(pathItems.rect.idx)

df <- data.frame(xmx = 1+cs.rect, xmi = c(1,1+cs.rect[-length(cs.rect)]))
pal <- colorRampPalette(c('#8B0000','#c23b22','#ffffff','#779ecb','#00008b'))
pal <- pal(100)
p1 <- ggplot() + 
  geom_tile(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
  # boxes on diagonal only
  geom_rect(data=df,aes(xmin = xmi, xmax = xmx, ymin = nrow(A)-xmi+1, ymax = nrow(A)-xmx +1),fill=NA,color='black',size=0.4) +
  # box off every feature combo
  #geom_rect(data=df,aes(xmin = xmi, xmax = xmi, ymin = 0, ymax = nrow(A)),fill=NA,color='grey80',size=0.2) +
  #geom_rect(data=df,aes(xmin = 0, xmax = ncol(A), ymin = nrow(A)-xmi+1, ymax = nrow(A)-xmx+1),fill=NA,color='grey80',size=0.2) +
  scale_y_discrete(labels = pathItems.labels,expand=c(0,0)) + 
  scale_x_discrete(labels = pathItems.labels[length(pathItems.labels):1],expand=c(0,0)) + coord_fixed() +
  scale_fill_gradientn(name = "r",limits = c(-1,1),breaks=c(-1,0,1),colours = pal,na.value = 'grey80') + 
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=0.5,size = 6,color='black'),
        axis.text.y = element_text(size = 6,color='black'), text = element_text(size = 6,color='black'),
        axis.ticks = element_blank(),legend.key.height = unit(0.15,'in'),legend.key.width = unit(0.05,'in'),
        axis.line = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1,linetype='solid'))
p1

# by region

pathRegions.name <- list("Ci","OC","SM","MF","An","CS","EC","DG","Am","TS","CP","GP","SN","LC","Me","CB","Po","MB")
pathRegions.prettylab <- list("Cing","OC","SM","MF","Ang","CA1/Sub","EC","DG","Amyg","TS","CP","GP","SN","LC","Med","CB","Pons","Mesenc.")
pathRegions.index <- sapply(1:length(pathRegions.name), function(i) grep(pathRegions.name[[i]], substr(colnames(microSample),start = 1,stop=2)))
pathRegions.labels <- sapply(1:length(pathRegions.name), function(i) c(matrix("",floor(0.5*length(pathRegions.index[[i]]))),
                                                                       pathRegions.prettylab[[i]], c(matrix("",ceiling(0.5*length(pathRegions.index[[i]])-1)))))
pathRegions.index <- Reduce(c,pathRegions.index)
pathRegions.labels <- Reduce(c,pathRegions.labels)

#pathItems.region <- list("Neocortical","MF","SN","Subcortical","Amyg","Brainstem","Ang")
microSample.regionorder <- microSample[,pathRegions.index]

A <- rcorr(as.matrix(microSample.regionorder),type='spearman')
p <- A$P
# remove NAs -- happens when not enough obs present across subjects for feature pair
# or when only remaining obs for a feature pair are constant
# set NAs to 0 and ignore p-value

p.calc.fdr <- na.exclude(p[upper.tri(p)])
p.thrsh <- compute.FDR(p.calc.fdr,q=0.05) # fdr adjust for all unique comparisons
p
p <- p < p.thrsh
A <- A$r*p
A[diag(1,nrow = nrow(A))] <- 1 # make diag 1
#A[is.na(A)] <- 0

melted_cormat <- melt(A)
melted_cormat$Var1 <- melted_cormat$Var1[length(melted_cormat$Var1):1]

# rectangle around each region
pathRegion.rect.idx <- unlist(lapply(1:length(pathRegions.name), function(i) length(grep(pathRegions.name[[i]],substr(colnames(microSample.regionorder),start=1,stop=2)))))
pathRegion.rect.idx <- fliplr(pathRegion.rect.idx)
cs.rect <- cumsum(pathRegion.rect.idx)
df.rect <- data.frame(xmx = 1+cs.rect, xmi = c(1,1+cs.rect[-length(cs.rect)]))
# demarcate brainstem, subcortex, cortex
border.regions <- sapply(c('SN','TS','CS'), function(r) which(fliplr(pathRegions.name)==r))
df.border <- data.frame(x = df.rect$xmx[border.regions])

pal <- colorRampPalette(c('#8B0000','#c23b22','#ffffff','#779ecb','#00008b'))
pal <- pal(100)

p2 <- ggplot() + 
  geom_tile(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
  # boxes on diagonal only
  geom_rect(data=df.rect,aes(xmin = xmi, xmax = xmx, ymin = nrow(A)-xmi+1, ymax = nrow(A)-xmx +1),fill=NA,color='black',size=0.4) +
  geom_vline(data=NULL,xintercept=df.border$x,linetype=1,size=0.4) +
  geom_hline(data=NULL,yintercept=nrow(A) - df.border$x+1,linetype=1,size=0.4) +
  # box off every feature combo
  #geom_rect(data=df,aes(xmin = xmi, xmax = xmi, ymin = 0, ymax = nrow(A)),fill=NA,color='grey80',size=0.2) +
  #geom_rect(data=df,aes(xmin = 0, xmax = ncol(A), ymin = nrow(A)-xmi+1, ymax = nrow(A)-xmx+1),fill=NA,color='grey80',size=0.2) +
  scale_y_discrete(labels = pathRegions.labels,expand=c(0,0)) + 
  scale_x_discrete(labels = pathRegions.labels[length(pathRegions.labels):1],expand=c(0,0)) + coord_fixed() +
  scale_fill_gradientn(name = "r",limits = c(-1,1),breaks=c(-1,0,1),colours = pal,na.value='grey80') + 
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=0.5,size = 6,color='black'),
        axis.text.y = element_text(size = 6,color='black'), text = element_text(size = 6,color='black'),
        axis.ticks = element_blank(),legend.key.height = unit(0.15,'in'),legend.key.width = unit(0.05,'in'),
        axis.line = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1,linetype='solid'))
p2
p1p2 <- plot_grid(plotlist = list(p1,p2),align="hv",nrow=1)
ggsave(p1p2,filename = file.path(savedir,'neuropathcorrTypeRegion.pdf'),units = 'in',height = 3.5,width = 7)

##############

# List highest values by region, excluding angiopathy, gliosis, and neuron loss
pathItems.type <- list("Tau","Thio","Ubiquitin","Antibody","TDP43","Syn")
pathItems.index <- sapply(1:length(pathItems.type), function(i) grep(pathItems.type[[i]], colnames(microSample)))
pathItems.index <- Reduce(c,pathItems.index)
microSample.typeorder <- microSample[,pathItems.index]

pathRegions.name <- list("Ci","OC","SM","MF","An","CS","EC","DG","Am","TS","CP","GP","SN","LC","Me","CB","Po","MB")
pathRegions.prettylab <- list("Cing","OC","SM","MF","Ang","CA1/Sub","EC","DG","Amyg","TS","CP","GP","SN","LC","Med","CB","Pons","Mesenc.")
pathRegions.index <- lapply(1:length(pathRegions.name), function(i) grep(pathRegions.name[[i]], substr(colnames(microSample.typeorder),start = 1,stop=2)))

# for each region, examine co-occurrence between difference forms of neuropath
copath.by.region <- lapply(pathRegions.index, function(idx) 
  cor(microSample.typeorder[,idx],method = 'spearman',use= 'pairwise.complete.obs') * !diag(length(idx)))  #remove diag
p.by.region <- lapply(pathRegions.index, function(idx) 
  (rcorr(as.matrix(microSample.typeorder[,idx]),type='spearman')$P * !diag(length(idx))) < p.thrsh)  #remove diag
# get index of path type pair for each region with highest abs. correlation
path.pairs <- list()
pair.cors <- list()
top.n <- 3
for(R in 1:length(pathRegions.index)){
  # get top n correlation values
  cp.mat <- copath.by.region[[R]] * p.by.region[[R]]
  top.copath <- tail(sort(unique(as.vector(abs(cp.mat)))),n=top.n)
  # find index
  top.idx <- t(sapply(top.copath, function(CP) which(abs(cp.mat) == CP,arr.ind=T)[1,]))
  # get name pair
  path.pairs[[R]] <- sapply(1:nrow(top.idx), function(i) 
    paste(rownames(copath.by.region[[R]])[top.idx[i,1]],'-',colnames(copath.by.region[[R]])[top.idx[i,2]]))
  # get r value
  pair.cors[[R]] <- sapply(1:nrow(top.idx), function(i) 
    copath.by.region[[R]][top.idx[i,1],top.idx[i,2]])
}

idx.x <- unlist(lapply(pathRegions.prettylab,function(R) rep(R,top.n)))
idx.y <- unlist(lapply(pathRegions.prettylab, function(R) 3:1))
path.pairs <- unlist(path.pairs)
pair.cors <- unlist(pair.cors)

pal <- colorRampPalette(c('red','white','blue'))(100)
p <- ggplot() + geom_text(aes(x=idx.y,y=idx.x,label=path.pairs,color=pair.cors),size=2.5) +
  scale_x_continuous(limits=c(0.5,3.5)) + #scale_color_gradientn(colors = pal,name='r') +
  scale_y_discrete(limits=idx.x) +
  scale_color_viridis(option = 'plasma',name='r') + theme_classic() +
  theme(axis.line.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.title.y = element_blank(),legend.key.size = unit(0.1,units = 'in'),
        legend.text = element_text(size=6),text = element_text(size=8))
p
ggsave(plot = p,filename = file.path(savedir,'CopathologyByRegionWords.pdf'),height = 3,width = 8,units = 'in')

idx.max.copath <- lapply(1:length(pathRegions.index), function(R)
  which(abs(copath.by.region[[R]]) == max(abs(copath.by.region[[R]]),na.rm = T),arr.ind=T)[1,])
# get name of path involved in pair
path.pairs <- lapply(1:length(pathRegions.index), function(R)
  paste(rownames(copath.by.region[[R]])[idx.max.copath[[R]][1]],'-',
        colnames(copath.by.region[[R]])[idx.max.copath[[R]][2]]))
# get value of correlation
pair.cors <- lapply(1:length(pathRegions.index), function(R)
  copath.by.region[[R]][idx.max.copath[[R]][1],idx.max.copath[[R]][2]])
