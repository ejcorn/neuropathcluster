# same as CharacterizeLouvainClusters.R but stratifies AD and LBD by disease stages

rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'analyzecluster/patientcharacteristics/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,-(1:2)] # Get rid of index column and INDDIDs
patientSample.orig <- read.csv(paste(params$opdir,'processed/patientSampleABC.csv',sep=''),stringsAsFactors = F)[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''),stringsAsFactors=F)[,2]

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

microSample <- remove.Disconnected.Subjects(microSample,DisconnectedSubjects)
patientSample.orig <- remove.Disconnected.Subjects(patientSample.orig,DisconnectedSubjects)

############################
### Cluster by diagnoses ###
############################

list[patientSample,dz.short]<- other.dz(patientSample.orig)

# add in graded colors for AD stages -- first get colors for each disease
pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
dz.colors <- pal(length(dz.short))
AD.Color <- dz.colors[which(dz.short == 'AD')] # then get AD color
LBD.Color <- dz.colors[which(dz.short == 'LBD')] # and get LBD color
# then make 3 colors for int high that get progressively lighter (because high, int is alphabetical and dz.short must match sort(unique(patientSample$NPDx1)))
AD.Stage.Colors <- c(darken(col = AD.Color,amount = 0.1,method='relative'),AD.Color,lighten(col = AD.Color,amount = 0.3,method='relative'))
#AD.Stage.Colors <- AD.Stage.Colors[1:length(unique(patie))]
# insert new colors
dz.colors <- c(dz.colors[0:(which(dz.short == 'AD')-1)],AD.Stage.Colors,dz.colors[(which(dz.short == 'AD')+1):length(dz.colors)])
# replace AD diagnoses with ADNC levels
patientSample$NPDx1[patientSample$NPDx1 %in% c('Alzheimer\'s disease')] <- paste('ADNC -',patientSample$ADStatus[patientSample$NPDx1 == 'Alzheimer\'s disease'])
# if ADNC - None shows up as primary dx... that means Alzheimer's disease listed but ABC criteria not met ... so AD listing probably a mistake AND probably no other pertinent diagnoses
# this only occurs in one patient and they are essentially normal except with very mild tau
patientSample$NPDx1[patientSample$NPDx1 == 'ADNC - None'] <- 'Unremarkable adult' 
list[patientSample,dz.short]<- other.dz(patientSample)
# make 3 colors for brainstem limbic neocortical that get darker (b/c alphabetical order, see above)
LBD.Stage.Colors <- c(lighten(col = LBD.Color,amount = 0.3,method='relative'),lighten(col = LBD.Color,amount = 0.2,method='relative'),LBD.Color,darken(col=LBD.Color,amount=0.1,method='relative'))
# insert new colors
dz.colors <- c(dz.colors[0:(which(dz.short == 'LBD')-1)],LBD.Stage.Colors,dz.colors[(which(dz.short == 'LBD')+1):length(dz.colors)])
# replace LBD diagnoses with LBD patterns
patientSample$NPDx1[patientSample$NPDx1 == 'LBD'] <- paste('LBD -',patientSample$DLBType[patientSample$NPDx1 == 'LBD'])

######################################
### Plot distribution of diagnoses ###
######################################

p1 <- ggplot(patientSample) + geom_bar(aes(x=NPDx1,fill=NPDx1)) + theme_classic()+
  scale_fill_manual(limits=sort(unique(patientSample$NPDx1)),values = dz.colors) + 
  theme(text=element_text(size=8),axis.text.x = element_text(size=7,angle=90,hjust=1,vjust=0.5),legend.position = 'none') +
  ylab('# of Patients') + xlab('') 

df.plt <- patientSample.orig[patientSample$NPDx1 == 'Other',]
df.plt$NPDx1[df.plt$NPDx1 == 'Alzheimer\'s disease' & df.plt$ADStatus == 'None'] <- 'Unremarkable adult'
p2 <- ggplot(df.plt)+ geom_bar(aes(x=NPDx1,fill=NPDx1)) + theme_classic()+
  scale_fill_brewer(limits=sort(unique(df.plt$NPDx1)),palette='Blues') + 
  ggtitle(paste('Other: n =',nrow(df.plt))) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(text=element_text(size=8),axis.text.x = element_text(size=8,angle=90,hjust=1,vjust=0.5),legend.position = 'none') +
  ylab('# of Patients') + xlab('')

df.plt <- patientSample.orig[patientSample$NPDx1 == 'Tau-Other',]
df.plt$NPDx1[df.plt$NPDx1 == 'Alzheimer\'s disease' & df.plt$ADStatus == 'None'] <- 'Unremarkable adult'
p3 <- ggplot(df.plt)+ geom_bar(aes(x=NPDx1,fill=NPDx1)) + theme_classic()+
  scale_fill_brewer(limits=sort(unique(df.plt$NPDx1)),palette='Blues') + 
  ggtitle(paste('Tau-Other: n =',nrow(df.plt))) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(text=element_text(size=8),axis.text.x = element_text(size=8,angle=90,hjust=1,vjust=0.5),legend.position = 'none') +
  ylab('# of Patients') + xlab('')

p.list <- list(p1,p2,p3)
p.all <- plot_grid(plotlist = p.list,nrow=1,rel_widths = c(2,2,1))
ggsave(filename = paste(savedir,'PatientDiagnosesOtherTauOther.pdf',sep=''),plot = p.all,
       height = 6,width=18,units='cm')

FigSXac <- lapply(p.list,function(X) X$data) # get source data
