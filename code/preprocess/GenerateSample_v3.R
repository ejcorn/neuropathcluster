rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=TRUE)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/processfxns.R')

################################
### Load and Preprocess Data ###
################################

data.date <- '122219'
micro <- read.csv(file = paste0("data/INDD_Micro",data.date,".csv"),stringsAsFactors = FALSE)
patients <- read.csv(file = paste0("data/INDD_Patients",data.date,".csv"),stringsAsFactors = FALSE)
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

# process micro column names corresponding to region and type of pathology
# names are currently RegionType ... convert to Region_Type in order to unambiguously identify in future
# also shorten type names for plotting
list[pathItems.type,pathRegions.name] <- get.pathscore.names()
list[pathItems.type.short,pathRegions.name.short] <- get.pathscore.names(vers='short')
cnames <- colnames(microSample)
newnames <- rep(NA,length(cnames)) # hold new names
for(type.i in 1:length(pathItems.type)){
  type <- pathItems.type[[type.i]]
  type.short <- pathItems.type.short[[type.i]]
  type.idx <- grep(type,cnames) # where are items of that type located
  type.startstop <- lapply(cnames[type.idx], function(X) list(st.sp=str_locate_all(X,type)[[1]],nm = X)) # get indices of where type is located within string
  newnames.tmp <- sapply(type.startstop, function(X) paste0(substr(X$nm,1,X$st.sp[,'start']-1),'_',type.short)) # make new names 
  if(length(newnames.tmp)>0){colnames(microSample)[type.idx] <- newnames.tmp} # replace existing column names with new names
}

############################################################
### plot missing data in a region-by-feature type matrix ###
############################################################

missingMatrix <- is.na(microSample[,-1])
patients$AutopsyDate <- as.Date(patients$AutopsyDate,format = '%m/%d/%Y')
date.breaks.p <- seq.Date(as.Date('1990/01/01'),as.Date('2020/01/01'),by = '5 years')
df.plot <- data.frame(x=patients$AutopsyDate[patients$INDDID %in% microSample$INDDID],y=rowSums(missingMatrix))
p <- ggplot() + geom_vline(xintercept = as.Date('2007/01/01'),linetype='dashed',color='grey70')+
  geom_point(aes(x=x,y=y),alpha=0.5,stroke=0,size=0.75) + 
  scale_x_date(breaks=date.breaks.p,labels = year(date.breaks.p))+
  xlab('Autopsy Date') + ylab('Missing Features') + theme_classic() +
  theme(text=element_text(size=8),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
ggsave(p,filename = paste0(savedir,'MissingDataByAutopsyDate.pdf'),units= 'cm',height = 4.5,width=4.5)
save(df.plot, file=paste(savedir,'FigS1a_SourceData.RData',sep=''))

# remove pre 2007 cases
date.thresh <- as.Date('01/01/2007',format='%m/%d/%Y') 
INDDIDs.pre2007 <- patients$INDDID[patients$AutopsyDate < date.thresh]
microSample <- microSample[!microSample$INDDID %in% INDDIDs.pre2007,]
missingMatrix <- is.na(microSample[,-1])

missingMatrix.reshape <- region.by.item.matrix(missingMatrix)
clim <- c(0,max(100*missingMatrix.reshape,na.rm=T))
p <- imagesc(100*missingMatrix.reshape,cmap='Blues',clim = clim) + ggtitle('Missing Data (% Subjects)')+
  theme(text=element_text(size=6),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),plot.title = element_text(hjust=0.5,size=6)) + nice_cbar() 
ggsave(p,filename = paste0(savedir,'MissingDataRegionTypeAll.pdf'),units= 'cm',height = 4.5,width=5.3)
save(missingMatrix.reshape, file=paste(savedir,'FigS1b_SourceData.RData',sep=''))

# Look at the above plot. You will see that data is ubiquitously missing from a set of
# regions and types of features moreso than other features. They are:
# 'Ubiquitin','OC','OFC','MC','DG','Antibody', and 'LC'
# Additionally, 'CB_Angiopathy' has a lot of missing data
# These features will be excluded in the next step

featureTags <- c('Ubiquitin','OC','OFC','MC','DG','Antibody','LC','CB_Angiopathy')
retainedFeatureMask <- rowSums(sapply(featureTags, function(tag) grepl(tag,colnames(missingMatrix))))==0
missingMatrix.reshape <- region.by.item.matrix(missingMatrix[,retainedFeatureMask])
clim <- c(0,max(100*missingMatrix.reshape,na.rm=T))
p<-imagesc(100*missingMatrix.reshape,cmap='Blues',clim=clim) + ggtitle('Missing Data (% Subjects)')+
  theme(text=element_text(size=6),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),plot.title = element_text(hjust=0.5,size=6)) + nice_cbar() 
ggsave(p,filename = paste0(savedir,'MissingDataRegionTypeExcluded.pdf'),units= 'cm',height = 4.5,width=5.3)
save(missingMatrix.reshape, file=paste(savedir,'FigS1c_SourceData.RData',sep=''))

microSample <- microSample[,c(TRUE,retainedFeatureMask)] # add TRUE to keep INDDID column

#########################################################################
### look at characteristic of subject retention vs. feature retention ###
#########################################################################

DataRepresentation <- sapply(seq(0,1,length.out = 100), function(i) sum((rowSums(is.na(microSample)) < i*ncol(microSample))))
# Y axis shows percent of subjects with less than X% of features containing missing data
p1 <- ggplot() + geom_point(aes(x = seq(0,100,length.out = 100), y = 100*DataRepresentation / nrow(microSample))) +
  geom_vline(xintercept = 100*params$missing.thrsh.r,linetype='dashed',color='grey70') + 
  xlab('max. % of features missing') + ylab('% of subjects') + theme_classic() + ggtitle('Data Representation') +
  theme(plot.title = element_text(size = 8,hjust=0.5),text=element_text(size=8))
p1
ggsave(p1,filename = paste(savedir,'DataRepresentationSubjects.pdf',sep=''),units= 'cm',height = 4.5,width=4.5)
save(DataRepresentation, file=paste(savedir,'FigS1d_SourceData.RData',sep=''))

DataRepresentation <- sapply(seq(0,1,length.out = 100), function(i) 100*mean((colSums(is.na(microSample)) < i*nrow(microSample))))
# Y axis shows percent of features with less than X% of subjects containing missing data
p1 <- ggplot() + geom_point(aes(x = seq(0,100,length.out = 100), y = DataRepresentation)) +
  xlab('max. % of subjects missing') + ylab('% of features') + theme_classic() + ggtitle('Data Representation') +
  theme(plot.title = element_text(size =8,hjust=0.5),text=element_text(size=8))
p1
ggsave(p1,filename = paste(savedir,'DataRepresentationFeatures.pdf',sep=''),units= 'cm',height = 4.5,width=4.5)
save(DataRepresentation, file=paste(savedir,'FigS1b_SourceData.RData',sep=''))

###################################################
### Exclude subjects with a lot of missing data ###
###################################################

# given that we must retain enough features to broadly characterize copathology
# exclude subjects with more than X% of those features missing

missing.thrsh.r <- params$missing.thrsh.r
missing.Mask <- rowMeans(is.na(microSample[,-1])) < missing.thrsh.r
microSample <- microSample[missing.Mask,]

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

# save micro sample
write.csv(x = microSample,file = paste(savedir,"microSample.csv",sep=''))

# Narrow down clinical diagnostic labels

n.dx <- 2
ClinDx.all <- sapply(1:n.dx, function(i) paste('ClinicalDx',i,sep=''))
for(ClinDx.i in ClinDx.all){
  patients[,ClinDx.i][which(grepl("Alzheimer",patients[,ClinDx.i]))] <- "Alzheimer's disease"
  patients[,ClinDx.i][which(grepl("Corticobasal",patients[,ClinDx.i]))] <- "CBD"
  patients[,ClinDx.i][which(grepl("supranuclear",patients[,ClinDx.i]))] <- "PSP"
  patients[,ClinDx.i][which(grepl("Dementia with Lewy",patients[,ClinDx.i]))] <- "DLB"
  patients[,ClinDx.i][which(grepl("PPA",patients[,ClinDx.i]))] <- "PPA"
  patients[,ClinDx.i][which(grepl("bvFTD",patients[,ClinDx.i]))] <- "bvFTD"
  patients[,ClinDx.i][which(grepl("Parkinson\'s Disease with Dementia",patients[,ClinDx.i]))] <- "Parkinson's disease dementia"
  patients[,ClinDx.i][which(grepl("Parkinson\'s Disease (not demented)",patients[,ClinDx.i]) | grepl("Parkinsonism NOS",patients[,ClinDx.i]) | grepl("Parkinson\'s disease",patients[,ClinDx.i]))] <- "Parkinson's disease"
  patients[,ClinDx.i][which(grepl("Amyotrophic",patients[,ClinDx.i]))] <- "Amyotrophic Lateral Sclerosis"
  patients[,ClinDx.i][which(grepl("Mild",patients[,ClinDx.i]))] <- "MCI"
  patients[,ClinDx.i][which(grepl("Impaired, not MCI",patients[,ClinDx.i]))] <- "Generally Impaired"
  patients[,ClinDx.i][which(grepl("Multiple system",patients[,ClinDx.i]))] <- "Multiple system atrophy"
  patients[,ClinDx.i][which(grepl("ascular",patients[,ClinDx.i]))] <- "Cerebrovascular Disease"
  patients[,ClinDx.i][which(grepl("chizoph",patients[,ClinDx.i]) | grepl("Depression",patients[,ClinDx.i]) | grepl("psychiatric",patients[,ClinDx.i]))] <- "Psychiatric Illness"
  #patients[,ClinDx.i][which(grepl("Unremarkable adult",patients[,ClinDx.i]))] <- "Normal"
}

#Narrow down disease labels for primary - quaternary NPDx

n.dx <- 5
NPDx.all <- sapply(1:n.dx, function(i) paste('NPDx',i,sep=''))
for(NPdx.i in NPDx.all){
  patients[,NPdx.i][which(grepl("Alzheimer",patients[,NPdx.i]))] <- "Alzheimer's disease"
  patients[,NPdx.i][which(grepl("Corticobasal",patients[,NPdx.i]))] <- "CBD"
  patients[,NPdx.i][which(grepl("Chronic Traumatic Encephalopathy",patients[,NPdx.i]))] <- "CTE"  
  patients[,NPdx.i][which(grepl("supranuclear",patients[,NPdx.i]))] <- "PSP"
  
  patients[,NPdx.i][which(grepl("Lewy body disease",patients[,NPdx.i]))] <- "LBD"
  patients[,NPdx.i][which(grepl("Parkinson\'s disease, atypical",patients[,NPdx.i]))] <- "LBD"
  # per JR and Eddie Lee's emails on 1-21-20, group together all lewy body disease into one
  # patients[,NPdx.i][which(grepl("dementia with Lewy",patients[,NPdx.i]))] <- "LBD"
  # patients[,NPdx.i][which(grepl("Parkinson's disease dementia",patients[,NPdx.i]))] <- "LBD"
  # patients[,NPdx.i][which(grepl("Parkinson's disease",patients[,NPdx.i]))] <- "Parkinson's disease"
  
  patients[,NPdx.i][which(grepl("PPA",patients[,NPdx.i]))] <- "PPA"
  patients[,NPdx.i][which(grepl("amyloid angiopathy",patients[,NPdx.i]))] <- "CAA"
  patients[,NPdx.i][which(grepl("FTLD-TDP",patients[,NPdx.i]))] <- "FTLD-TDP"
  patients[,NPdx.i][which(grepl("Frontotemporal",patients[,NPdx.i]))] <- "FTLD-Other"
  patients[,NPdx.i][which(grepl("Limbic-predominant",patients[,NPdx.i]))] <- "LATE"
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
  patients[is.na(patients[,NPdx.i]),NPdx.i] <- '' # some missing values are NA, others are ''... make all '' for future logical comparisons
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

# perform alzheimer's staging per Montine et al. 2012 Table 3 -- also see ADNPCGrid.png in data/ folder
# march down grid from top then left to right
patients$ADStatus <- '' # make new column
patients$ADStatus[patients$ABeta == 0 & patients$CERAD == 0] <- 'None'
patients$ADStatus[patients$ABeta == 1 & patients$CERAD <= 1] <- 'Low'
patients$ADStatus[patients$ABeta == 1 & patients$CERAD >= 2 & patients$Braak03 <= 1] <- 'Low'
patients$ADStatus[patients$ABeta == 1 & patients$CERAD >= 2 & patients$Braak03 >= 1] <- 'Intermediate'
patients$ADStatus[patients$ABeta == 2 & patients$Braak03 <= 1] <- 'Low'
patients$ADStatus[patients$ABeta == 2 & patients$Braak03 >= 1] <- 'Intermediate'
patients$ADStatus[patients$ABeta == 3 & patients$CERAD <= 1 & patients$Braak03 <= 1] <- 'Low'
patients$ADStatus[patients$ABeta == 3 & patients$CERAD <= 1 & patients$Braak03 >= 1] <- 'Intermediate'
patients$ADStatus[patients$ABeta == 3 & patients$CERAD >= 2 & patients$Braak03 <= 1] <- 'Low'
patients$ADStatus[patients$ABeta == 3 & patients$CERAD >= 2 & patients$Braak03 == 2] <- 'Intermediate'
patients$ADStatus[patients$ABeta == 3 & patients$CERAD >= 2 & patients$Braak03 == 3] <- 'High'
patients$ADStatus[grepl('Definite',patients$NPDx1Likelihood) & grepl('Alzheimer\'s disease',patients$NPDx1)] <- 'High'

# in case ABeta is missing, do it based on Braak and CERAD only
patients$ADStatus[patients$ABeta == '' & patients$Braak03 == 0 & patients$CERAD == 0] <-'Low'
patients$ADStatus[patients$ABeta == '' & (patients$Braak03 < 2 | patients$CERAD < 2) & (patients$Braak03 != 0 | patients$CERAD != 0)] <-'Low'
patients$ADStatus[patients$ABeta == '' & patients$Braak03 >= 2 & patients$CERAD >= 2] <-'Intermediate'
patients$ADStatus[patients$ABeta == '' & patients$Braak03 == 3 & patients$CERAD == 3] <-'High'

# process DLB stages
patients$DLBType[patients$DLBType == 'Diffuse or Neocortical'] <- 'Neocortical'
patients$DLBType[patients$DLBType == 'Brainstem Predominant'] <- 'Brainstem'
patients$DLBType[patients$DLBType == 'Transitional or Limbic'] <- 'Limbic'
patients$DLBType[patients$DLBType == 'Amygdala Predominant'] <- 'Amygdala'
patients$DLBType[is.na(patients$DLBType)] <- '?'
patients$DLBType[patients$DLBType == ''] <- '?'
patients$DLBType[patients$DLBType == 'N/A'] <- '?'

# x<- grepl.npdx(patients,'LBD',5)
# patients$DLBType[x]
# MissingDLBType <- patients[x & !patients$DLBType %in% c('Neocortical','Brainstem','Amygdala','Limbic'),]
# write.csv(x=MissingDLBType,row.names = F,file = paste0(savedir,'MissingDLBType_EJC012420.csv'))
# JR.LBD <- read.csv('data/LBD cases 1-23-2020.csv',stringsAsFactors = F)
# sum(JR.LBD$INDDID %in% MissingDLBType$INDDID)

write.csv(x = patients,file = paste(savedir,"patientSample.csv",sep=''))

# Make a separate patient sample where NPDx1 contains Alzheimer's disease wherever a patient has Int/High ABC sstage

dz.exc <- 'Alzheimer\'s disease'
n.dx <- 5
Rem.Mask <- exclude.dz(patients,dz.exc,n.dx) # get mask for all patients with ABC stage int-high

ABC <- !Rem.Mask
# insert.ABC makes NPDX1 = 'Alzheimer\'s disease' if ABC stage is intermediate-high
# and all other diagnoses are shifted behind it in NPDx2-5
patients <- insert.ABC(patients,ABC)

print('Alzheimers in NPDx1 is equal to ABC mask?')
print(identical(patients$NPDx1 == dz.exc,ABC))
NPDx.all <- sapply(2:5,function(n) as.character(get(paste0('NPDx',n),patients)))
print(paste('AD in NPDx2-5:',sum(NPDx.all=='Alzheimer\'s disease',na.rm = T)))

write.csv(x = patients,file = paste(savedir,"patientSampleABC.csv",sep=''))

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

cog <- read.csv(paste0('data/INDD_MoCA',data.date,'.csv'),stringsAsFactors = F)
cog <- cog[cog$INDDID %in% patients$INDDID,]
cog <- cog[order(cog$INDDID),]
# remove any MOCAs that are all 0's or all NA
cog <- cog[cog$MoCATotal != 0 & !is.na(cog$MoCATotal),]
# remove any patients whose visuospatial scores is > 5 (one oddly had 995)
cog <- cog[-which(cog$VisuospatialTotal > 5),]
write.csv(x = cog,file = paste(savedir,"MoCA_processed.csv",sep=''))

Genes <- read.csv(paste0('data/INDD_Genetics',data.date,'.csv'),stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% patients$INDDID,]
Genes <- Genes[order(Genes$INDDID),]
Genes$APOE[Genes$APOE == 'E2/E3  '] <- 'E2/E3'
Genes$APOE[Genes$APOE == 'E3/E4  '] <- 'E2/E4'
write.csv(x = Genes, file = paste(savedir,'Genetics_processed.csv',sep=''))

CSF.names <- c('Luminex','Elisa')
for(CSF.name in CSF.names){
  CSF <- read.csv(paste0('data/INDD_CSF',CSF.name,data.date,'.csv'),stringsAsFactors = F)
  CSF <- CSF[CSF$INDDID %in% patients$INDDID,]
  CSF <- process.CSF(CSF,CSF.name)

  write.csv(x=CSF, file = paste(params$opdir,'processed/CSF',CSF.name,'_processed.csv',sep=''))
}