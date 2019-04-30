#do PCA, correlation matrix on neuropathology
library(tidyverse)
library(viridis)
library(reshape2)
library(glmnet)
library(caret)
library(scatterplot3d)
library(igraph)
library(blockmodels)
library(R.matlab)
library(lm.beta)
rm(list = ls())
homedir <- '~/Dropbox/Neurodegeneration/PathCogClinDx/'
setwd(homedir)

microSample <- read.csv('data/microSample.csv')[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv('data/patientSample.csv')[,-(1:2)] # Get rid of index column and INDDIDs
INDDIDs <- read.csv('data/microSample.csv')[,2]

if(sum((colSums(is.na(microSample)) == nrow(microSample))) > 0){
  break
}

fisher.r.to.z <- function(r){
  r <- 0.5*(log(1+r) - log(1-r))
  return(r)
}

partial.resid <- function(mdl,b){
  # for lm object mdl, compute partial residuals WRT named coefficient (string in b)
  s <- summary(mdl)$coef
  try(p.r <- residuals(mdl) + s[b,'Estimate']*mdl$model[,b],silent = T)
  try(p.r <- residuals(mdl) + s[paste('`',b,'`',sep=''),'Estimate']*mdl$model[,b],silent=T)
  return(p.r)
}

###########################################################
### Can you explain MOCA with copath vs. individ. path? ###
###########################################################

cog <- read.csv('data/MOCA.csv',stringsAsFactors = F)
cog <- cog[cog$INDDID %in% INDDIDs,]
meanMOCA <- sapply(unique(cog$INDDID), function(ID) mean(cog$MoCATotal[cog$INDDID == ID],na.rm= T))
microSample <- microSample[INDDIDs %in% unique(cog$INDDID),]
patientSample <- patientSample[INDDIDs %in% unique(cog$INDDID),]
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
copath.by.patient <- fisher.r.to.z(copath.by.patient)
copath.by.patient[copath.by.patient == Inf] <- NA # remove any cors. exactly = to 1
#get rid of diagonal (1) and keep only utri (ltri is duplicate because symmetric cor-mat)
utri <- which(as.vector(t(upper.tri(matrix(1,nTypes,nTypes)))))
copath.by.patient <- as.data.frame(copath.by.patient[,utri])

d <- which(duplicated(copath.by.patient$`Syn-TDP43`))
copath.by.patient$`Syn-TDP43`[d]
fisher.r.to.z(cor(type.by.patient[[d[1]]],method = 'spearman',use='pairwise.complete.obs'))

########################################
### Examine MOCAs by path and copath ###
########################################

# Plot vascular pathology by single path #
ind.path <- as.data.frame(lapply(pathItems.type, function(ptype) data.frame(rowMeans(na.rm=T,microSample[,grep(ptype,colnames(microSample))]))))
colnames(ind.path) <- pathItems.type
m <- lm.beta(lm(meanMOCA ~ .,data=ind.path))
ip.tab <- summary(m)$coef[-1,]

df <- data.frame(n=rownames(ip.tab),b=ip.tab[,'Standardized'])
sig <- ifelse(p.adjust(ip.tab[,'Pr(>|t|)'],method = 'fdr') < 0.05, yes= '*',no='')
p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#5F4B8B',color='black',alpha=0.3) +
  ylab(expression(beta)) + xlab('') + ggtitle('MOCA') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.1*max(df$b),label=sig) +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5))
p
ggsave(plot = p,filename = paste('results/MOCAVsPath.pdf',sep=''),height = 1.5,width = 2.25,units = 'in')

df <- data.frame(x=m$model[,'TDP43'],y=partial.resid(mdl = m,b='TDP43'))
ggplot(data=df,aes(x=x,y=y)) + geom_point()

# Regress out single path, plot vasc path by copath #
df.full <- cbind(ind.path,copath.by.patient)
m2 <- lm.beta(lm(meanMOCA ~ .,data=df.full))
ip.tab <- summary(m2)$coef[-1,]
sig <- ifelse(p.adjust(ip.tab[,'Pr(>|t|)'],method = 'fdr') < 0.05, yes= '*',no='')
sig <- sig[!rownames(ip.tab) %in% colnames(ind.path)]
ip.tab <- ip.tab[!rownames(ip.tab) %in% colnames(ind.path),]
df <- data.frame(n=rownames(ip.tab),b=ip.tab[,'Standardized'])

sig <- ifelse(p.adjust(ip.tab[,'Pr(>|t|)'],method = 'fdr') < 0.05, yes= '*',no='')
p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill='#5F4B8B',color='black',alpha=0.3) +
  ylab(expression(beta)) + xlab('') + ggtitle('MOCA') + theme_classic() +
  annotate(geom='text',x=df$n,y= 1.1*max(df$b),label=sig) +
  theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),
        axis.text.x = element_text(angle=45,vjust=0.5))
p
ggsave(plot = p,filename = paste('results/MOCAVsCoPath.pdf',sep=''),height = 2,width = 2.25,units = 'in')

df <- data.frame(x=m2$model[,'Syn-TDP43'],y=partial.resid(mdl = m2,b='Syn-TDP43'))
ggplot(data=df,aes(x=x,y=y)) + geom_point()

cor(ind.path,copath.by.patient,use='pairwise.complete.obs')
