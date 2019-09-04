rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=TRUE)
source('code/misc/fxns.R')

################################
### Load and Preprocess Data ###
################################

micro <- read.csv(file = "data/INDD_Micro12119.csv",stringsAsFactors = FALSE)
patients <- read.csv(file = "data/INDD_Patients12119.csv",stringsAsFactors = FALSE)
micro <- micro[order(micro$INDDID),] # make sure data is ordered by INDDID
patients <- patients[order(patients$INDDID),]

## Remove qualitative or unusable columns
micro <- micro[,-grep("Other",colnames(micro))]
micro <- micro[,-grep("Update",colnames(micro))]
micro <- micro[,-grep("Specify",colnames(micro))]
micro <- micro[,-grep("SynucleinStaining",colnames(micro))]
micro <- micro[,-grep("Neocortical",colnames(micro))]
micro <- micro[,-grep("Brainstem",colnames(micro))]
micro <- micro[,-grep("Subcortical",colnames(micro))]
micro <- micro[,-grep('SC',colnames(micro))] # remove spinal cord because poorly represented in normals
## Remove subjects with old pathology scoring system
# For these subjects, scores exceed 5+
# It's possible the below method will fail to eliminate a subject with 
# old system scores < 5 for all regions assessed
# will update but likely won't affect


# micro > 5 will leave in "Not Done", "Not Avail", "Rare", "Presumed 0"
# don't exclude now as those will be dealt with later
microNumeric <- as.data.frame(lapply(micro, function(micro) as.numeric(micro))) # this will only capture the numeric old scores
#oldScoringSystemMask <- sapply(1:ncol(micro), function(C) micro[,C] > 5 & !grepl("Not Done",micro[,C]) & 
#                                 !grepl("Rare",micro[,C]) & !grepl("Not Avail",micro[,C]) & !grepl("Presumed",micro[,C]))
oldScoringSystemMask <- microNumeric != 0  # only 0 is part of the new system 
oldScoringSystemMask[,1] <- FALSE # Don't count INDDIDs, which are all > 5
oldScoringSystemMask <- !(rowSums(oldScoringSystemMask,na.rm = T) > 0) # get subjects utilizing old score system

micro <- micro[oldScoringSystemMask,]

#*** PERHAPS ADD SUBST. NIG. LOC. CER. from patients data *****#

# Replace 0,1+ etc with 1:10 to make scoring sequence quantitative

pathscores.qual <- c("0","Rare","1+","2+","3+")
# only want things that are 0, Rare, 1+, 2+, 3+
pathscores.qual <- pathscores.qual[length(pathscores.qual):1]
pathscores.quant <- seq(1:length(pathscores.qual))[length(pathscores.qual):1]

# make new matrix microSample full of NAs, then fill in iteratively
microSample <- as.data.frame(matrix(NA,nrow=nrow(micro),ncol=ncol(micro)))
colnames(microSample) <- colnames(micro)
microSample$INDDID <- micro$INDDID # copy INDDIDs
for(C in colnames(micro)[-1]){ # iterate through columns, excluding INDDIDs
  # call presumed 0 a 1 (lowest value in new ordinal scale)
  microSample[which(grepl("Presumed",micro[,C])),C] <- 1
  # make not available or not done NAs
  micro[which(grepl("Not",micro[,C])),C] <- NA
  # replace qualitative score with new quantitative score
  for(P in 1:length(pathscores.qual)){
    microSample[which(micro[,C]==pathscores.qual[P]),C] <- pathscores.quant[P]
  } 
  microSample[,C] <- as.numeric(microSample[,C])
}

DataRepresentation <- sapply(seq(0,1,length.out = 100), function(i) sum((rowSums(is.na(microSample)) < i*ncol(microSample))))
# Y axis shows percent of subjects with less than X% of features containing missing data
p1 <- ggplot() + geom_point(aes(x = seq(0,100,length.out = 100), y = 100*DataRepresentation / nrow(microSample))) +
  xlab('max. % of features missing') + ylab('% of subjects') + theme_classic() + ggtitle('Data Representation') +
  theme(plot.title = element_text(size = 8,hjust=0.5),text=element_text(size=8))
p1
ggsave(p1,filename = paste(savedir,'DataRepresentationSubjects.pdf',sep=''),units= 'cm',height = 4.5,width=4.5)
save(DataRepresentation, file=paste(savedir,'FigS1a_SourceData.RData',sep=''))

DataRepresentation <- sapply(seq(0,1,length.out = 100), function(i) 100*mean((colSums(is.na(microSample)) < i*nrow(microSample))))
# Y axis shows percent of features with less than X% of subjects containing missing data
p1 <- ggplot() + geom_point(aes(x = seq(0,100,length.out = 100), y = DataRepresentation)) +
  xlab('max. % of subjects missing') + ylab('% of features') + theme_classic() + ggtitle('Data Representation') +
  theme(plot.title = element_text(size =8,hjust=0.5),text=element_text(size=8))
p1
ggsave(p1,filename = paste(savedir,'DataRepresentationFeatures.pdf',sep=''),units= 'cm',height = 4.5,width=4.5)
save(DataRepresentation, file=paste(savedir,'FigS1b_SourceData.RData',sep=''))
# Remove subjects with too many NAs -- discuss what this really ought to be given the graph above

missing.thrsh.r <- params$missing.thrsh.r
missing.Mask <- rowMeans(is.na(microSample[,-1])) < missing.thrsh.r
microSample <- microSample[missing.Mask,]

# Remove features with too many NAs -- discuss what this ought to be given the graph above
missing.thrsh.c <- params$missing.thrsh.c
missing.Mask <- colMeans(is.na(microSample)) < missing.thrsh.c
microSample <- microSample[,missing.Mask]

#####################################################
### Generate sample focusing on diagnostic subset ###
#####################################################

#remove patients without micro data
patients <- patients[patients$INDDID %in% microSample$INDDID,]

#remove patients without patient data
microSample <- microSample[microSample$INDDID %in% patients$INDDID,]

# remove certain diseases

#X<- patients[grepl('schizophrenia',patients$NPDx1),]
#X2 <- patients[grepl('Unremarkable adult',patients$NPDx1),]
#X2 <- microSample[grepl('Unremarkable adult',patients$NPDx1),]
rem.dz <- list('Fetal brain') #,'Unremarkable adult')
patientMask <- lapply(rem.dz, function(X) grepl(X,patients$NPDx1))
patientMask <- Reduce("+",patientMask) == 0 & !patients$NPDx1 == ''
# find schizophrenics without additional neuropathological diagnoses
#schizMask <- grepl('schizophrenia',patients$NPDx1) & patients$NPDx2 == ''
#patientMask <- patientMask & !schizMask

microSample <- microSample[patientMask,]
patients <- patients[patientMask,]

write.csv(x = microSample,file = paste(savedir,"microSample.csv",sep=''))

# Narrow down clinical diagnostic labels

n.dx <- 2
ClinDx.all <- sapply(1:n.dx, function(i) paste('ClinicalDx',i,sep=''))
for(ClinDx.i in ClinDx.all){
  patients[,ClinDx.i][which(grepl("Alzheimer",patients[,ClinDx.i]))] <- "Alzheimer's disease"
  patients[,ClinDx.i][which(grepl("Corticobasal",patients[,ClinDx.i]))] <- "CBD"
  patients[,ClinDx.i][which(grepl("supranuclear",patients[,ClinDx.i]))] <- "PSP"
  patients[,ClinDx.i][which(grepl("Dementia with Lewy",patients[,ClinDx.i]))] <- "LBD"
  patients[,ClinDx.i][which(grepl("PPA",patients[,ClinDx.i]))] <- "PPA"
  patients[,ClinDx.i][which(grepl("bvFTD",patients[,ClinDx.i]))] <- "bvFTD"
  patients[,ClinDx.i][which(grepl("Parkinson",patients[,ClinDx.i]))] <- "Parkinson's disease"
  patients[,ClinDx.i][which(grepl("Amyotrophic",patients[,ClinDx.i]))] <- "Amyotrophic Lateral Sclerosis"
  patients[,ClinDx.i][which(grepl("Mild",patients[,ClinDx.i]))] <- "MCI"
  patients[,ClinDx.i][which(grepl("Impaired, not MCI",patients[,ClinDx.i]))] <- "Generally Impaired"
  patients[,ClinDx.i][which(grepl("Multiple",patients[,ClinDx.i]))] <- "Multiple System Atrophy"
  patients[,ClinDx.i][which(grepl("ascular",patients[,ClinDx.i]))] <- "Cerebrovascular Disease"
  patients[,ClinDx.i][which(grepl("chizoph",patients[,ClinDx.i]) | grepl("Depression",patients[,ClinDx.i]) | grepl("psychiatric",patients[,ClinDx.i]))] <- "Psychiatric Illness"
  #patients[,ClinDx.i][which(grepl("Unremarkable adult",patients[,ClinDx.i]))] <- "Normal"
}

#Narrow down disease labels for primary - quaternary NPDx

n.dx <- 4
NPDx.all <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
for(NPdx.i in NPDx.all){
  patients[,NPdx.i][which(grepl("Alzheimer",patients[,NPdx.i]))] <- "Alzheimer's disease"
  patients[,NPdx.i][which(grepl("Corticobasal",patients[,NPdx.i]))] <- "CBD"
  patients[,NPdx.i][which(grepl("supranuclear",patients[,NPdx.i]))] <- "PSP"
  patients[,NPdx.i][which(grepl("dementia with Lewy",patients[,NPdx.i]))] <- "LBD"
  patients[,NPdx.i][which(grepl("PPA",patients[,NPdx.i]))] <- "PPA"
  patients[,NPdx.i][which(grepl("amyloid angiopathy",patients[,NPdx.i]))] <- "CAA"
  patients[,NPdx.i][which(grepl("FTLD-TDP",patients[,NPdx.i]))] <- "FTLD-TDP"
  patients[,NPdx.i][which(grepl("Frontotemporal",patients[,NPdx.i]))] <- "FTLD-Other"
  patients[,NPdx.i][which(grepl("Parkinson's disease dementia",patients[,NPdx.i]))] <- "LBD" #"Parkinson's Disease Dementia"
  patients[,NPdx.i][which(grepl("Parkinson's disease",patients[,NPdx.i]))] <- "Parkinson's disease"
  patients[,NPdx.i][which(grepl("Lewy",patients[,NPdx.i]))] <- "LBD"
  patients[,NPdx.i][which(grepl("Amyotrophic",patients[,NPdx.i]))] <- "Amyotrophic Lateral Sclerosis"
  patients[,NPdx.i][which(grepl("Mild",patients[,NPdx.i]))] <- "MCI"
  patients[,NPdx.i][which(grepl("Impaired, not MCI",patients[,NPdx.i]))] <- "Generally Impaired"
  patients[,NPdx.i][which(grepl("Multiple",patients[,NPdx.i]))] <- "Multiple System Atrophy"
  patients[,NPdx.i][which(grepl("ascular",patients[,NPdx.i]))] <- "Cerebrovascular Disease"
  patients[,NPdx.i][which(grepl("Prion",patients[,NPdx.i]))] <- "Prion Disease"
  patients[,NPdx.i][which(grepl("Primary age-related",patients[,NPdx.i]))] <- "PART"
  patients[,NPdx.i][which(grepl("Pathological ",patients[,NPdx.i]))] <- "Pathological Aging"
  patients[,NPdx.i][which(grepl("chizoph",patients[,NPdx.i]))] <- "Schizophrenia"
  patients[,NPdx.i][which(grepl("Unremarkable adult",patients[,NPdx.i]))] <- "Unremarkable adult"
  #patients[,NPdx.i][which(grepl("chizoph",patients[,NPdx.i]) | grepl("Depression",patients[,NPdx.i]) | grepl("psychiatric",patients[,NPdx.i]))] <- "Psychiatric Illness"
  
}

# there is one patient who has Alzheimer's listed as NPDx1 and NPDx2
DoubleDx <- 120215
patients$NPDx2[which(patients$INDDID == DoubleDx)] <- ''
patients$NPDx2Likelihood[which(patients$INDDID == DoubleDx)] <- ''

# process Braak scores. Braak03 = floor(Braak06/2) for those without Braak03
patients$Braak03 <- as.numeric(patients$Braak03)
patients$Braak06 <- as.numeric(patients$Braak06)
patients$Braak03[is.na(patients$Braak03)] <- floor(patients$Braak06[is.na(patients$Braak03)]/2)

write.csv(x = patients,file = paste(savedir,"patientSample.csv",sep=''))

# Replace Alzheimer's disease with Braak-CERAD scores

dz.exc <- 'Alzheimer\'s disease'
n.dx <- 4
Rem.Mask <- exclude.dz(patients,dz.exc,n.dx) # get mask for all patients meeting Braak-CERAD AD

BraakCERAD <- !Rem.Mask
patients <- insert.BraakCERAD(patients,BraakCERAD)

print('Alzheimers in NPDx1 is equal to BraakCERAD mask?')
print(identical(patients$NPDx1 == dz.exc,BraakCERAD))

write.csv(x = patients,file = paste(savedir,"patientSampleBraakCERAD.csv",sep=''))

# visualize representation of diseases
p <- ggplot(data = patients,aes(x=NPDx1,fill=NPDx1)) + geom_hline(yintercept = 100,color ='grey')+ 
  geom_bar() + theme_classic() +
  scale_y_continuous(expand= c(0,0))+
  ylab('Count') + xlab('Primary Histological diagnosis') + theme(text = element_text(size =8),
                                                             axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5)) +
  theme(legend.position = 'none')
p
ggsave(p,filename = paste(savedir,'SampleCharacteristics.pdf',sep=''),units= 'in',height = 3,width=4)

########################################
### Process gene, CSF, and MoCA data ###
########################################

cog <- read.csv('data/INDD_MoCA12119.csv',stringsAsFactors = F)
cog <- cog[cog$INDDID %in% patients$INDDID,]
cog <- cog[order(cog$INDDID),]
# remove any MOCAs that are all 0's or all NA
cog <- cog[cog$MoCATotal != 0 & !is.na(cog$MoCATotal),]
write.csv(x = cog,file = paste(savedir,"MoCA_processed.csv",sep=''))

Genes <- read.csv('data/INDD_Genetics12119.csv',stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% patients$INDDID,]
Genes <- Genes[order(Genes$INDDID),]
Genes$APOE[Genes$APOE == 'E2/E3  '] <- 'E2/E3'
Genes$APOE[Genes$APOE == 'E3/E4  '] <- 'E2/E4'
write.csv(x = Genes, file = paste(savedir,'Genetics_processed.csv',sep=''))

CSF.name <- 'Luminex'
CSF <- read.csv(paste('data/INDD_CSF',CSF.name,'12119.csv',sep=''),stringsAsFactors = F)
CSF <- CSF[CSF$INDDID %in% patients$INDDID,]
#specify which CSF features you care about
CSF.vars <- c('LuminexTTau','LuminexPTau','LuminexAbeta42')
#CSF.vars <- c('TotalTau','PhosphorylatedTau','Beta.amyloid42')
#CSF.vars <- c('ElisaTTau','ElisaPTau','ElisaAbeta42')

CSF <- CSF[!rowSums(is.na(CSF[,CSF.vars])),] # remove missing
CSF <- CSF[order(CSF$INDDID),]

CSF$CSFDate <- as.Date(CSF$CSFDate,format = '%m/%d/%Y')

write.csv(x=CSF, file = paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''))
