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

############################################################
### Cellular path vs. molecular path at individual level ###
############################################################

Cellular <- c('Angiopathy','NeuronLoss','Gliosis')
cols <- c(Angiopathy='#E41A1C',NeuronLoss='#377EB8',Gliosis='#4DAF4A')
pathItems.type <- list("Tau","Thio","Antibody","TDP43","Syn")
# molecular path, averaged across region
ind.path <- as.data.frame(lapply(pathItems.type, function(ptype) data.frame(rowMeans(na.rm=T,microSample[,grep(ptype,colnames(microSample))]))))
colnames(ind.path) <- pathItems.type
load(file=paste(savedir,'CopathByPatient.RData',sep=''))

nboot <- 10

res.molcell <- res.comolcell <- list()

for(Cell in Cellular){
  # data frame of each type of cellular path, averaged across regions
  cell.path <- data.frame(y = rowMeans(na.rm = T,microSample[,grep(Cell,colnames(microSample))]))  
  df <- cbind(cell.path,ind.path)
  list[p.vals,ul.ll,m] <- boot.reg(df,nboot)
  res.molcell[[Cell]] <- list(p=p.vals,ul.ll,m)  

  # account for copathology as well
  df <- cbind(cell.path,ind.path,copath.by.patient)
  list[p.vals,ul.ll,m] <- boot.reg(df,nboot)
  res.comolcell[[Cell]] <- list(p=p.vals,ul.ll,m)
}

# FDR correct over all cellular features
p.vals.molcell <- lapply(res.molcell, function(X) X[['p']]) # get p-vals into list
p.vals.molcell <- list.fdr.correct(p.vals.molcell) # correct over all cellular features
p.vals.comolcell <- lapply(res.comolcell, function(X) X[['p']]) # get p-vals into list
p.vals.comolcell <- list.fdr.correct(p.vals.comolcell) # correct over all cellular features

# make plots with corrected p-values
p.molcell <- p.comolcell <- list()
for(Cell in Cellular){
  list[p.vals,ul.ll,m] <- res.molcell[[Cell]]
  p.molcell[[Cell]] <- plot.beta.bar(p.vals.molcell[[Cell]],ul.ll,m,Cell,cols[Cell]) 
  list[p.vals,ul.ll,m] <- res.comolcell[[Cell]]
  p.comolcell[[Cell]] <- plot.beta.bar(p.vals.comolcell[[Cell]],ul.ll,m,Cell,cols[Cell]) + theme(axis.text.x = element_text(angle=90,vjust=0.5))  
}

p.molcell <- plot_grid(plotlist=p.molcell,align='hv',nrow=1)
p.comolcell <- plot_grid(plotlist=p.comolcell,align='hv',nrow=1)
ggsave(plot = p.molcell,filename = paste(savedir,'MolecularVsCellular.pdf',sep=''),height = 4,width = 18,units = 'cm')
ggsave(plot = p.comolcell,filename = paste(savedir,'ComolecularVsCellular.pdf',sep=''),height = 7,width = 18,units = 'cm')
