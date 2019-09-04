# look at non-AD diagnoses in Cluster 2 and differences in age at death within cluster 2

rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'analyzecluster/',sep='')
dir.create(savedir,recursive=TRUE)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

################################
### Load and Preprocess Data ###
################################

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

demo <- read.csv(paste(homedir,'data/INDD_GlobalDemographics6819.csv',sep=''))

demo <- demo[demo$INDDID %in% INDDIDs,]	# only look at patients in our sample

age <- demo$AgeatDeath # extract age at death

###############################################
### Examine Braak-CERAD vs. 1' histopath dx ###
###############################################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
partitionSample <- sapply(partition, function(i) paste('Cluster',i))

list[patientSample,dz.short]<- other.dz(patientSample)
clusterColors <- getClusterColors(k)

dz.exc <- 'Alzheimer\'s disease'
exc.cl <- 2
n.dx <- 4
Rem.Mask <- exclude.dz(patientSample,dz.exc,n.dx) # mask all patients meeting Braak-CERAD AD

# Verify that all remaining subjects don't have CERAD & Braak > 1
print(paste(sum(patientSample$Braak03[Rem.Mask] > 1 & patientSample$CERAD[Rem.Mask] > 1,na.rm=T),'patients remaining with Braak & CERAD > 1'))
print(sum(patientSample$Braak03[Rem.Mask] > 1,na.rm=T))

# excludes disease based on diagnosis that is > "Low Probability"
NPDx <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
NPDx <- lapply(NPDx,function(x) get(x,patientSample))
NPDx <- lapply(NPDx, function(X) X == dz.exc)  

BraakMask <- patientSample$Braak03 < 2 # patients with low Braak scores
CERADMask <- patientSample$CERAD < 2 # patients with low CERAD scores
missingBraak <- is.na(patientSample$Braak03)
missingCERAD <- is.na(patientSample$CERAD)

# only look at people with Braak and CERAD scores present
Rem.Mask <- rep(NA,nrow(patientSample)) # initialize mask to select patients without Braak-CERAD AD
Rem.Mask[BraakMask] <- TRUE # if patients have low Braak scores they can't have Braak-CERAD AD
Rem.Mask[CERADMask] <- TRUE # if patients have low CERAD scores they can't have Braak-CERAD AD
Rem.Mask[!BraakMask & !CERADMask] <- FALSE # if patients have high Braak & CERAD scores they have Braak-CERAD AD
# confirm:

print(paste(sum(patientSample$Braak03[Rem.Mask]>1 & patientSample$CERAD[Rem.Mask]>1,na.rm=T),
	'patients with Braak-CERAD AD in mask supposed to highlight patients without Braak-CERAD AD'))
print(paste(sum(patientSample$Braak03[!Rem.Mask]<2 & patientSample$CERAD[!Rem.Mask]<2,na.rm=T),
	'patients without Braak-CERAD AD in mask supposed to highlight patients with Braak-CERAD AD'))
print(paste(sum(is.na(Rem.Mask)),'patients who cannot be evaluated with Braak-CERAD AD criteria'))

# Now all other remaining patients are NAs because they have missing scores
# This is either patients missing both Braak & CERAD, or with high Braak or CERAD missing the other score
# such that you cannot rule out Braak-CERAD AD
# in these 58 patients we will diagnose them with AD based on NPDx1-n.dx having > low probability
missingBraakCERAD <- is.na(Rem.Mask)

NPDx.lik <- sapply(1:n.dx, function(i) paste('NPDx',i,'Likelihood',sep=''))
NPDx.lik <- lapply(NPDx.lik,function(x) get(x,patientSample))
NPDx.lik <- lapply(NPDx.lik, function(X) X == "Definite (High Probability)" | X == "Probable (Intermediate Probability)")

# NOT anyone with non-low probability dx of dz.exc...TRUE is only low or no prob of dz.exc. false is > low prob
NPDx.comb <- mapply(function(NPDx,NPDx.lik){list(!(NPDx & NPDx.lik))},NPDx=NPDx,NPDx.lik=NPDx.lik)

Lik.Mask <- Reduce('&',NPDx.comb) # mask selecting patients with only low prob or no diagnosis of dz.exc
Dz.Mask <- Reduce('|',NPDx)  # mask selecting patients with any diagnosis of dz.exc in NPDx1 to n.dx

# Does anybody who is missing either Braak or CERAD have an AD dx with greater than low probability?
print(paste(sum(!Lik.Mask[missingBraakCERAD]),' out of ',length(Lik.Mask[missingBraakCERAD]),
	' patients (',100*mean(!Lik.Mask[missingBraakCERAD]),'%) missing Braak or CERAD w/ > low prob of ',dz.exc,sep=''))
# yes... for 22 of 58 patients we will be saying they have AD

# How well do NPDxLikelihood line up with Braak-CERAD? Does anybody not meeting Braak-CERAD criteria have high prob of AD?

print(paste(sum(!Lik.Mask[Rem.Mask & !missingBraakCERAD]),' out of ',sum(Rem.Mask & !missingBraakCERAD),
	' non-Braak-CERAD AD patients (',100*mean(!Lik.Mask[Rem.Mask & !missingBraakCERAD]),
	'%) with w/ > low prob of ',dz.exc,sep=''))

# Does everybody meeting Braak-CERAD criteria have a high prob of AD

print(paste(sum(!Lik.Mask[!Rem.Mask & !missingBraakCERAD]),' out of ',sum(!Rem.Mask &!missingBraakCERAD),
	' Braak-CERAD AD patients (',100*mean(!Lik.Mask[!Rem.Mask & !missingBraakCERAD]),
	'%) with w/ > low prob of ',dz.exc,sep=''))

BraakCerad <- factor(ifelse(!Rem.Mask[!missingBraakCERAD],yes='AD',no='non-AD'))
HighLikelihood <- factor(ifelse(!Lik.Mask[!missingBraakCERAD],yes='AD',no='non-AD'))
print('Confusion matrix for NPDx Diagnoses predicting Braak-CERAD scores')
print(confusionMatrix(data=HighLikelihood,reference=BraakCerad))
# NPDx diagnoses are 87% sensitive and 96% specific. This means you are more likely to miss Braak-CERAD AD, but not overestimate it
# So when we “impute” the Braak-CERAD scores for these 58 patients, we may be missing some of them with AD.

# if missing CERAD or Braak, use subjective likelihood
Rem.Mask[missingBraakCERAD] <- Lik.Mask[missingBraakCERAD]

####################
### Age at death ###
####################

Cluster.Mask <- partitionSample == paste('Cluster',exc.cl) # look only at Cluster of interest

NonAD.Cluster.Mask <- Rem.Mask & Cluster.Mask

# load the patientSample that has Braak-CERAD AD applied
patientSample.BraakCERAD <- read.csv(paste(params$opdir,'processed/patientSampleBraakCERAD.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
diagnoses <- patientSample.BraakCERAD$NPDx1[NonAD.Cluster.Mask]

# see how many patients have low-level AD change (CERAD > 1
noise.x <- rnorm(sum(NonAD.Cluster.Mask),sd=0.1)
noise.y <- rnorm(sum(NonAD.Cluster.Mask),sd=0.1)
plot(patientSample$Braak03[NonAD.Cluster.Mask] + noise.x,patientSample$CERAD[NonAD.Cluster.Mask] + noise.y)

# get the raw diagnoses of this group before applying the Braak-CERAD criteria
raw.diagnoses <- patientSample$NPDx1[NonAD.Cluster.Mask]
raw.likelihood <- patientSample$NPDx1Likelihood[NonAD.Cluster.Mask]
raw.likelihood[raw.diagnoses == dz.exc]
# ^--- *some* the patients in cluster 2 without are like the "false positives" of the 
# the majority are other diagnoses and largest constituent is PSP

# so these should be the patients who are FALSE for the likelihood mask (indicating high estimated probability of AD)
# and TRUE for Rem.Mask (indicating negative for Braak-CERAD)
# question is: are these just early AD patients? or on the border? I bet they have 

CNDR.vs.BC.FalsePositives <- !Lik.Mask & Rem.Mask
print(paste((100*sum(CNDR.vs.BC.FalsePositives[NonAD.Cluster.Mask]) / sum(CNDR.vs.BC.FalsePositives)),
	'% of CNDR \'false positives\' are in non-AD cluster 2',sep=''))
print(paste(100*mean(CNDR.vs.BC.FalsePositives[NonAD.Cluster.Mask]),'% of non-AD Cluster 2 are CNDR \'false positives\''))

p <- ggplot() + geom_bar(aes(x=diagnoses)) + theme_classic() +
  scale_y_continuous(expand= c(0,0))+ ylab('Count') + xlab('Primary Histological diagnosis') + 
  theme(text = element_text(size =8), axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5)) +
  theme(legend.position = 'none')
p

age.nonad <- age[Rem.Mask & Cluster.Mask] # age of non-AD cluster 2
age.ad <- age[!Rem.Mask & Cluster.Mask] # age of AD cluster 2
age.ad <- age.ad[age.ad != 0 & !is.na(age.ad)] # get rid of missing data
age.nonad <- age.nonad[age.nonad != 0 & !is.na(age.nonad)] # get rid of missing data


