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

p.calc.fdr <- na.exclude(p[upper.tri(p)])
p.thrsh <- compute.FDR(p.calc.fdr,q=0.05) # fdr adjust for all unique comparisons
p <- p < p.thrsh
A <- A$r*p
A[diag(1,nrow = nrow(A))] <- 1 # make diag 1
A[is.na(A)] <- 0

Ctx.Limbic <- list("Ci","OC","SM","MF","An","CS","EC","DG","Am")
Ctx.Limbic.index <- unlist(sapply(Ctx.Limbic, function(X) grep(X, substr(colnames(microSample),start = 1,stop=2))))
Subctx.brainstem <- list("TS","CP","GP","SN","LC","Me","CB","Po","MB")
Subctx.brainstem.index <- unlist(sapply(Subctx.brainstem, function(X) grep(X, substr(colnames(microSample),start = 1,stop=2))))

A.ctx.limbic <- A[Ctx.Limbic.index,Ctx.Limbic.index] 
#A.ctx.limbic <- fisher.r.to.z(A.ctx.limbic)
A.subctx.brainstem <- A[Subctx.brainstem.index,Subctx.brainstem.index]
#A.subctx.brainstem <- fisher.r.to.z(A.subctx.brainstem)
A.ctx.limbic <- A.ctx.limbic[upper.tri(A.ctx.limbic)]
A.subctx.brainstem <- A.subctx.brainstem[upper.tri(A.subctx.brainstem)]

wilcox.test(A.ctx.limbic,
            A.subctx.brainstem,
            conf.int=T)

wilcox.test(A.ctx.limbic[A.ctx.limbic!=0],
            A.subctx.brainstem[A.subctx.brainstem!=0],
            conf.int=T)

