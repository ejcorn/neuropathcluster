rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'micekmeans/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs
load(file = paste0(savedir,'microSampleImputedmiceRF.RData'))

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

#######################
### factor analysis ###
#######################

summary(microSample.imp)
microSample.imp <- complete(microSample.imp,1)
colnames(microSample.imp) <- colnames(microSample)

library(psych)
library(GPArotation)
nf <- 15 # number of factors to extract
fa.path <- fa(r=micro.order.by(as.matrix(microSample.imp)),cor = 'poly',nfactors = nf)
loadings <- fa.path$loadings # get loadings
scores <- fa.path$scores
explained <- fa.path$Vaccounted['Proportion Var',]

order.by.varexp <- order(explained,decreasing = T) # reorder factor loadings by variance explained
loadings <- loadings[,order.by.varexp]
scores <- scores[,order.by.varexp]
explained <- explained[order.by.varexp]

colnames(loadings) <- paste0('F',1:nf)
colnames(scores) <- paste0('F',1:nf)
names(explained) <- paste0('F',1:nf)

save(fa.path,loadings,scores,explained,file = paste0(savedir,'FactorAnalysis.RData'))

###########
### PCA ###
###########
W <- polychoric(microSample.imp)
pc <- principal(r=W$rho,rotate='none',nf)
explained <- pc$Vaccounted['Proportion Explained',]
scores <- factor.scores(microSample.imp,pc)
loadings <- scores$weights
scores <- scores$scores
colnames(scores) <- paste0('F',1:nf)
scores <- as.data.frame(scores)
save(pc,loadings,scores,explained,file = paste0(savedir,'PolychoricPCA.RData'))
ggplot() + geom_boxplot(aes(y=scores[,'F3'],x=patientSample$NPDx1)) + theme(axis.text.x=element_text(angle=90))
### below is not used ###

# imagesc(cor(scores),clim = c(-1,1))
# 
# CSF.name <- 'Luminex'
# CSF <- read.csv(paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''),stringsAsFactors = F)
# CSF <- CSF[CSF$INDDID %in% INDDIDs,]
# CSF.sample <- extract.CSF.sample(CSF,CSF.name,n.sample='first')
# 
# df.CSF <- merge(cbind(data.frame(INDDID=INDDIDs),fa.path$scores),CSF.sample,by='INDDID')
# imagesc(cor(df.CSF[,-1]))
# 
# ###########
# ### PCA ###
# ###########
# 
# pca.mdl <- princomp(micro.order.by(microSample.imp),scores = T)
# loadings <- pca.mdl$loadings[,1:nf] # get loadings
# scores <- pca.mdl$scores[,1:nf]
# colnames(loadings) <- paste0('PC',1:nf)
# colnames(scores) <- paste0('PC',1:nf)
# 
# p <- imagesc(t(loadings)) + theme(axis.text.x = element_text(angle=90,size=4,hjust=1,vjust=0.5))
# ggsave(filename = paste0(savedir,'PCALoadings.pdf'),plot = p,
#        height = 12,width=18,units='cm')
# 
# df.plot <- lapply(1:nf, function(fac) data.frame(dx=patientSample$NPDx1,Score=scores[,fac]))
# p.list <- lapply(df.plot, function(df) ggplot(df) + geom_boxplot(aes(x=dx,y=Score)) + theme_classic()+
#                    theme(text=element_text(size=8),axis.text.x = element_blank()) + xlab(''))
# p.list[[nf]] <- p.list[[nf]] + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
# p <- plot_grid(plotlist = p.list,ncol=1)
# ggsave(filename = paste0(savedir,'PCAScoresByDiagnosis.pdf'),plot = p,
#        height = 18,width=18,units='cm')
# 
# #########################
# ### PCA -- polychoric ###
# #########################
# 
# polycor <- polychoric(micro.order.by(microSample.imp))
# ed <- eigen(polycor$rho)
# loadings <- ed$vectors[,1:nf]
# colnames(loadings) <- paste0('PC',1:nf)
# rownames(loadings) <- colnames(micro.order.by(microSample.imp))
# p <- imagesc(t(loadings),clim=c(-0.25,0.25)) + theme(axis.text.x = element_text(angle=90,size=4,hjust=1,vjust=0.5))
# ggsave(filename = paste0(savedir,'PolychoricRMatrixEigenvectorLoadings.pdf'),plot = p,
#        height = 12,width=18,units='cm')
# 
# 
# ##################
# ### Clustering ###
# ##################
# library(cluster)
# subject.cormat <- polychoric(t(microSample.imp))
# imagesc(subject.cormat$rho)
# #clust <- hclust(dist(microSample.imp/max(microSample.imp)))
# X <- microSample.imp/max(microSample.imp)
# k.rng <- 2:15
# clust <- lapply(k.rng, function(k) kmeans(x = X,centers = k,nstart = 100))
# sil <- lapply(clust, function(K) silhouette(x=K$cluster,dist=dist(X)))
# sil.mean <- sapply(sil,function(S) mean(S[,3]))
# p <- ggplot() + geom_line(aes(x=k.rng,y=sil.mean)) + #theme_classic() +
#   scale_x_continuous(breaks = k.rng) + xlab('k') + ylab('Mean Silhouette')
# p
# clust <- clust[which(k.rng <= 6)]
# p.list <- lapply(clust,function(X) imagesc(t(fliplr(micro.order.by(X$centers)))) + 
#                    theme(axis.text.y = element_text(size=4,hjust=1,vjust=0.5)))
# p <- plot_grid(plotlist = p.list)
# ggsave(filename = paste0(savedir,'KMeansCentroids.pdf'),plot = p,
#        height = 18,width=18,units='cm')
# for(k in c(4,6)){
#   partition <- clust[[which(k.rng==k)]]
#   DisconnectedSubjects <- c()
#   list[pathItems.index,pathItems.labels] <- get.feature.labels(pathItems.type,colnames(microSample))
#   # order microSample columns to appear in Fig. 4d displaying cluster centroids
#   microSample <- microSample[,pathItems.index]
#   
#   save(partition,centroids,k,DisconnectedSubjects,pathItems.labels,
#        file = paste(savedir,'subjLouvainPartitionReordered.RData',sep=''))
# }
# 
# ####################################################
# ### Latent dirichlet allocation (topic modeling) ###
# ####################################################
# 
# library(topicmodels)
# 
# 
