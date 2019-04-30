library(tidyverse)
library(viridis)
library(reshape2)
library(glmnet)
library(caret)
library(R.matlab)
library(lm.beta)
library(brainwaver)
library(car)
library(boot)
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'copath/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,1]
#microSample <- as.data.frame(scale(microSample,center = T))
if(sum((colSums(is.na(microSample)) == nrow(microSample))) > 0){
  break
}

##########################################################################################
### Can you explain vascular path, neuron loss, gliosis with copath vs. individ. path? ###
##########################################################################################

pathItems.type <- list("Tau","Thio","Antibody","TDP43","Syn")
nTypes <- length(pathItems.type)
pathItems.index <- lapply(1:nTypes, function(i) grep(pathItems.type[[i]], colnames(microSample)))
type.by.patient <- lapply(1:nrow(microSample), function(Pt)
  do.call('cbind',lapply(pathItems.index, function(idx) as.numeric(microSample[Pt,idx]))))

copath.by.patient <- 
  lapply(type.by.patient, function(X) cor(X,method='spearman',use = 'pairwise.complete.obs'))
copath.by.patient <- do.call('rbind', lapply(copath.by.patient, function(X) as.vector(t(X))))

copath.names <- sapply(pathItems.type, function(A) # name each element for combo
  sapply(pathItems.type, function(B) paste(A,'-',B,sep='')))
colnames(copath.by.patient) <- as.vector(t(copath.names))
#get rid of diagonal (1) and keep only utri (ltri is duplicate because symmetric cor-mat)
utri <- which(as.vector(t(upper.tri(matrix(1,nTypes,nTypes)))))
copath.by.patient <- as.data.frame(copath.by.patient[,utri])
mean.copath <- colMeans(copath.by.patient) # mean copath 
sd.copath <- apply(X = copath.by.patient,MARGIN = 2,FUN='sd') # variance in copath

################
### Vascular ###
################

vasc.path <- data.frame(y = rowMeans(na.rm = T,microSample[,grep('Angiopathy',colnames(microSample))]))
ind.path <- as.data.frame(lapply(pathItems.type, function(ptype) data.frame(rowMeans(na.rm=T,microSample[,grep(ptype,colnames(microSample))]))))
colnames(ind.path) <- pathItems.type
df <- cbind(vasc.path,ind.path)

m <- lm(y ~ .,data=df)
df.boot <- boot(data=df,statistic=lm.boot,R=1000)

p.vals <- colMeans(df.boot$t < 0)
p.vals[p.vals > 0.5] <- colMeans(df.boot$t > 0)[p.vals > 0.5]
ul.ll <- t(sapply(1:ncol(df.boot$t), function(i) quantile(x = df.boot$t[,i],c(0.025,0.975))))

df <- data.frame(n=unlist(pathItems.type),b=coefficients(m)[-1],ul=ul.ll[,'97.5%'],ll=ul.ll[,'2.5%'])
sig <- ifelse(p.adjust(p.vals,method = 'fdr') < 0.025, yes= '*',no='')

p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#E41A1C',color='black',alpha=0.3) +
  geom_errorbar(aes(x=n,ymin=ll,ymax=ul),width=0.5)+
  ylab(expression(beta)) + xlab('') + ggtitle('Angiopathy') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.1*max(df$ul),label=sig,color='black') +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5))
p
ggsave(plot = p,filename = paste(savedir,'AngiopathyVsPath.pdf',sep=''),height = 1.5,width = 2.25,units = 'in')

df.full <- cbind(data.frame(y=vasc.path$y),ind.path,copath.by.patient)

m2 <- lm(y ~ .,data=df.full)
df.boot <- boot(data=df.full,statistic=lm.boot,R=1000)

p.vals <- colMeans(df.boot$t < 0)
p.vals[p.vals > 0.5] <- colMeans(df.boot$t > 0)[p.vals > 0.5]
ul.ll <- t(sapply(1:ncol(df.boot$t), function(i) quantile(x = df.boot$t[,i],c(0.025,0.975))))

df <- data.frame(n=names(coefficients(m2)[-1]),b=coefficients(m2)[-1],ul=ul.ll[,'97.5%'],ll=ul.ll[,'2.5%'])
sig <- ifelse(p.adjust(p.vals,method = 'fdr') < 0.025, yes= '*',no='')

p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#E41A1C',color='black',alpha=0.3) +
  geom_errorbar(aes(x=n,ymin=ll,ymax=ul),width=0.5)+
  ylab(expression(beta)) + xlab('') + ggtitle('Angiopathy') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.2*max(df$ul),label=sig) +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=0.5))
p
ggsave(plot = p,filename = paste(savedir,'AngiopathyVsCoPath.pdf',sep=''),height = 2,width = 2.25,units = 'in')

#################################################
### Test interactions between pathology types ###
#################################################

# y <- residuals(lm(vasc.path$Angiopathy ~ Thio + Syn + Antibody + Tau,data=ind.path))
# summary(lm.beta(lm(y~ Ubiquitin*TDP43,data=ind.path)))
# summary(lm.beta(lm(vasc.path$Angiopathy~ .,data=df.full)))

###################
### Neuron Loss ###
###################

neuron.death <- data.frame(y = rowMeans(na.rm=T,microSample[,grep('NeuronLoss',colnames(microSample))]))
df <- cbind(neuron.death,ind.path)
m <- lm(y ~ .,data=df)
df.boot <- boot(data=df,statistic=lm.boot,R=1000)

p.vals <- colMeans(df.boot$t < 0)
p.vals[p.vals > 0.5] <- colMeans(df.boot$t > 0)[p.vals > 0.5]
ul.ll <- t(sapply(1:ncol(df.boot$t), function(i) quantile(x = df.boot$t[,i],c(0.025,0.975))))

df <- data.frame(n=unlist(pathItems.type),b=coefficients(m)[-1],ul=ul.ll[,'97.5%'],ll=ul.ll[,'2.5%'])
sig <- ifelse(p.adjust(p.vals,method = 'fdr') < 0.025, yes= '*',no='')

p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#377EB8',color='black',alpha=0.3) +
  geom_errorbar(aes(x=n,ymin=ll,ymax=ul),width=0.5)+
  ylab(expression(beta)) + xlab('') + ggtitle('Neuron Loss') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.1*max(df$ul),label=sig) +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5))
p
ggsave(plot = p,filename = paste(savedir,'NeuronLossVsPath.pdf',sep=''),height = 1.5,width = 2.25,units = 'in')

df.full <- cbind(data.frame(y=neuron.death$y),ind.path,copath.by.patient)

m2 <- lm(y ~ .,data=df.full)
df.boot <- boot(data=df.full,statistic=lm.boot,R=1000)

p.vals <- colMeans(df.boot$t < 0)
p.vals[p.vals > 0.5] <- colMeans(df.boot$t > 0)[p.vals > 0.5]
ul.ll <- t(sapply(1:ncol(df.boot$t), function(i) quantile(x = df.boot$t[,i],c(0.025,0.975))))

df <- data.frame(n=names(coefficients(m2)[-1]),b=coefficients(m2)[-1],ul=ul.ll[,'97.5%'],ll=ul.ll[,'2.5%'])
sig <- ifelse(p.adjust(p.vals,method = 'fdr') < 0.025, yes= '*',no='')

p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#E41A1C',color='black',alpha=0.3) +
  geom_errorbar(aes(x=n,ymin=ll,ymax=ul),width=0.5)+
  ylab(expression(beta)) + xlab('') + ggtitle('Neuron Loss') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.2*max(df$ul),label=sig) +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=0.5))
p
ggsave(plot = p,filename = paste(savedir,'NeuronLossVsCoPath.pdf',sep=''),height = 2,width = 2.25,units = 'in')

###############
### Gliosis ###
###############

gliosis <- data.frame(y = rowMeans(na.rm=T,microSample[,grep('Gliosis',colnames(microSample))]))
df <- cbind(gliosis,ind.path)
m <- lm(y ~ .,data=df)
df.boot <- boot(data=df,statistic=lm.boot,R=1000)

p.vals <- colMeans(df.boot$t < 0)
p.vals[p.vals > 0.5] <- colMeans(df.boot$t > 0)[p.vals > 0.5]
ul.ll <- t(sapply(1:ncol(df.boot$t), function(i) quantile(x = df.boot$t[,i],c(0.025,0.975))))

df <- data.frame(n=unlist(pathItems.type),b=coefficients(m)[-1],ul=ul.ll[,'97.5%'],ll=ul.ll[,'2.5%'])
sig <- ifelse(p.adjust(p.vals,method = 'fdr') < 0.025, yes= '*',no='')

p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#4DAF4A',color='black',alpha=0.3) +
  geom_errorbar(aes(x=n,ymin=ll,ymax=ul),width=0.5)+
  ylab(expression(beta)) + xlab('') + ggtitle('Gliosis') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.2*max(df$ul),label=sig) +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5))
p
ggsave(plot = p,filename = paste(savedir,'GliosisVsPath.pdf',sep=''),height = 1.5,width = 2.25,units = 'in')

df.full <- cbind(data.frame(y=gliosis$y),ind.path,copath.by.patient)

m2 <- lm(y ~ .,data=df.full)
df.boot <- boot(data=df.full,statistic=lm.boot,R=1000)

p.vals <- colMeans(df.boot$t < 0)
p.vals[p.vals > 0.5] <- colMeans(df.boot$t > 0)[p.vals > 0.5]
ul.ll <- t(sapply(1:ncol(df.boot$t), function(i) quantile(x = df.boot$t[,i],c(0.025,0.975))))

df <- data.frame(n=names(coefficients(m2)[-1]),b=coefficients(m2)[-1],ul=ul.ll[,'97.5%'],ll=ul.ll[,'2.5%'])
sig <- ifelse(p.adjust(p.vals,method = 'fdr') < 0.025, yes= '*',no='')

p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#E41A1C',color='black',alpha=0.3) +
  geom_errorbar(aes(x=n,ymin=ll,ymax=ul),width=0.5)+
  ylab(expression(beta)) + xlab('') + ggtitle('Angiopathy') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.2*max(df$ul),label=sig) +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=0.5))
p
ggsave(plot = p,filename = paste(savedir,'GliosisVsCoPath.pdf',sep=''),height = 2,width = 2.25,  units = 'in')
