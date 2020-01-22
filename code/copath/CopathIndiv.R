rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'copath/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,1]

if(sum((colSums(is.na(microSample)) == nrow(microSample))) > 0){
  break
}

#########################################
### Compute copath within individuals ###
#########################################

pathItems.type <- list("Tau","Thio","Antibody","TDP43","Syn")
nTypes <- length(pathItems.type)
pathItems.index <- lapply(pathItems.type, function(X) grep(X, colnames(microSample)))
type.by.patient <- lapply(1:nrow(microSample), function(Pt)
  do.call('cbind',lapply(pathItems.index, function(idx) as.numeric(microSample[Pt,idx]))))
# Shuffle rrs w/o replacement
#type.by.patient <- lapply(type.by.patient, function(Pt)
#  do.call('cbind',lapply(1:ncol(Pt), function(i) sample(Pt[,i],replace = F))))

copath.by.patient <- 
  lapply(type.by.patient, function(X) rcorr(X,type='spearman'))

copath.by.patient <- do.call('rbind', lapply(copath.by.patient, function(X) as.vector(t(X$r))))

copath.names <- sapply(pathItems.type, function(A) # name each element for combo
  sapply(pathItems.type, function(B) paste(A,'-',B,sep='')))
colnames(copath.by.patient) <- as.vector(t(copath.names))
#get rid of diagonal (1) and keep only utri (ltri is duplicate because symmetric cor-mat)
utri <- which(as.vector(t(upper.tri(matrix(1,nTypes,nTypes)))))
copath.by.patient <- as.data.frame(copath.by.patient[,utri])

save(copath.by.patient,file=paste(savedir,'CopathByPatient.RData',sep=''))

###################################
### Test if indiv. copath r ≠ 0 ###
###################################

p.tests <- lapply(copath.by.patient, function(X) wilcox.test(X,conf.int=T))
#p.tests <- lapply(copath.by.patient, function(X) t.test(X))
p.pv <- sapply(p.tests, function(X) X$p.value)
p.esz <- sapply(p.tests, function(X) X$estimate)
sig <- p.signif(p.adjust(p.pv,method = 'fdr'))

ncp <- rep(colnames(copath.by.patient),each=nrow(copath.by.patient))
p <- ggplot() + geom_boxplot(aes(x = ncp,y= unlist(copath.by.patient)),fill='#5F4B8B',alpha=0.6,outlier.stroke = 0) +
  xlab('') + ylab('Inidividual Copathology (r)') + theme_classic() + 
  theme(text=element_text(size=8), axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) + 
  annotate(geom='text',x=names(sig),y=1.02,label=sig)
p
ggsave(plot = p,filename = paste(savedir,'IndividualCopathologyByType.pdf',sep=''),height = 2,width = 4,units = 'in')

# pairwise tests of differences in copathology

p.tests <- lapply(copath.by.patient, function(X) 
  lapply(copath.by.patient, function(Y) wilcox.test(X,Y,conf.int=T)))
p.pv <- sapply(p.tests, function(X) 
  sapply(X, function(Y) Y$p.value))
p.th <- compute.FDR(p.pv[upper.tri(p.pv)],q = 0.05)
p.esz <- sapply(p.tests, function(X) 
  sapply(X, function(Y) Y$estimate))
p.esz <- p.esz * (p.pv < p.th)
rownames(p.esz) <- colnames(copath.by.patient)
colnames(p.esz) <- colnames(copath.by.patient)

pdf(file = paste(savedir,"CopathologyDifferencesEffectSize.pdf",sep=''),width = unit(2.5,'inches'), height = unit(2.5,'inches'))
corrplot(p.esz,method = 'color',
         tl.col = 'black',
         tl.cex=0.5,
         cl.cex=0.5,
         cl.lim =(c(min(p.esz),max(p.esz))))
dev.off()

