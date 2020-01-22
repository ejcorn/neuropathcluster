rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=TRUE)

################################
### Load and Preprocess Data ###
################################

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics122219.csv',sep=''))

demo <- demo[demo$INDDID %in% INDDIDs,]	# only look at patients in our sample

age <- demo$AgeatDeath # extract age at death

age <- age[-which(age == 0 | is.na(age))] # remove entries with age = 0 or not available

print('Age At Death')
print(summary(age))
print(paste('SD of age at death:',sd(age)))

sex <- demo$Sex
print(paste('Female:',100*mean(sex=='Female'),'%'))
print(paste('Male:',100*mean(sex=='Male'),'%'))

NPDx1 <- patientSample$NPDx1
NPDx1.unique <- unique(NPDx1)
mean.age <- sapply(NPDx1.unique, function(dx) mean(demo$AgeatDeath[NPDx1==dx],na.rm=T)) # get mean age of death for each dx
df <- data.frame(y=mean.age,x=NPDx1.unique)
p <- ggplot(df) + geom_col(aes(x=x,y=y,fill=x)) + theme_classic() + xlab('') + ylab('Age at Death') +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),legend.position='none') + scale_y_continuous(expand= c(0,0))
p
ggsave(p,filename = paste(savedir,'AgeAtDeathByDxMean.pdf',sep=''),units= 'in',height = 4,width=4)

df <- data.frame(y=demo$AgeatDeath,x=NPDx1)
p <- ggplot(df) + geom_jitter(aes(x=x,y=y,color=x),alpha=0.5,stroke=0) + theme_classic() + xlab('') + ylab('Age at Death') +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),legend.position = 'none') + scale_y_continuous(expand= c(0,0))
p

ggsave(p,filename = paste(savedir,'AgeAtDeathByDxJitter.pdf',sep=''),units= 'in',height = 6,width=6)

# load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
# partition <- sapply(partition, function(i) paste('Cluster',i))
# 
# df <- data.frame(y=demo$AgeatDeath,x=partition)
# p <- ggplot(df) + geom_boxplot(aes(x=x,y=y,color=x),alpha=0.5) + theme_classic() + xlab('') + ylab('Age at Death') +
#   theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),legend.position = 'none') + scale_y_continuous(expand= c(0,0))
# p
