rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'micekmeans/PCA/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(params$opdir,'micekmeans/microSampleImputedmiceRF.RData'))
summary(microSample.imp)
microSample.imp <- complete(microSample.imp,1)
colnames(microSample.imp) <- colnames(microSample)

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]


###########
### PCA ###
###########
# load polychoric correlation matrix
load(file=paste0(params$opdir,'micekmeans/PolychoricFeatureCorrelations.RData'))
# perform pca on that correlation matrix using psych package
nf <- 15
pc <- principal(r=W$rho,rotate='none',nf)
explained <- pc$Vaccounted['Proportion Explained',]
# project the pathology onto those components
scores <- factor.scores(microSample.imp,pc)
loadings <- scores$weights
scores <- scores$scores
scores <- as.data.frame(scores)
save(pc,loadings,scores,explained,file = paste0(savedir,'PolychoricPCA.RData'))

##################################
### analyze component loadings ###
##################################
n.pcs <- 3
p.list <- lapply(1:n.pcs, function(PC) imagesc(region.by.item.matrix(t(loadings[,paste0('PC',PC)])),cmap = 'redblue_asymmetric')+
                   theme(legend.key.height=unit(0.3,'cm'),legend.key.width = unit(0.1,'cm'),legend.position = 'right')+
                   ggtitle(paste0('PC',PC,': ',expression(R^2),' = ',signif(explained[PC],2))) + 
                   theme(plot.title = element_text(hjust=0.5,size=8),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),plot.margin = unit(c(0,0,0,0),'cm')))
p <- plot_grid(plotlist = p.list,align = 'hv',nrow=1)
ggsave(filename = paste0(savedir,'PCALoadings.pdf'),plot = p,
       width = 20,height=3,units='cm')

# plot scores by diagnosis
list[patientSample,dz.short]<- other.dz(patientSample)
dx.order <- sort(unique(patientSample$NPDx1))
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.short))

pvals <- lapply(paste0('PC',1:n.pcs), function(PC) sapply(dx.order, function(dx)  # test if scores in each group differ from 0
  wilcox.test(scores[patientSample$NPDx1==dx,PC])$p.value))
pvals <- list.fdr.correct(pvals)
labs <- lapply(pvals,p.signif)
names(labs) <- paste0('PC',1:n.pcs)

p.list <- lapply(paste0('PC',1:n.pcs), function(PC) ggplot() + geom_boxplot(aes(y=scores[,PC],x=patientSample$NPDx1,fill=patientSample$NPDx1),outlier.size = 0.25,size=0.25) + 
                   theme(axis.text.x=element_text(angle=90)) + scale_x_discrete(limits=dx.order,labels=dz.short)+
                   annotate(geom='text',x=dx.order,y=Inf,vjust=1,label=labs[[PC]],size=2)+
                   theme_classic()+ scale_fill_manual(values=dz.colors)+xlab('') +ylab(PC)+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
                   scale_y_continuous(limits=max(abs(scores[,PC]))*c(-1,1))+theme(legend.position = 'none',text=element_text(size=8)))
p <- plot_grid(plotlist = p.list,align = 'hv',nrow=1)
ggsave(filename = paste0(savedir,'PCAScoresByDiagnosis.pdf'),plot = p,
       width = 18,height=3,units='cm')
