get.CSF.vars <- function(CSF.type){
	# INPUTS:
	# CSF.type: Luminex or Elisa, type of CSF assay
	# 
	# OUTPUTS:
	# column names corresponding to a-beta, phosphotau and total tau for given type of CSF assay
	#specify which CSF features you care about
	if(CSF.type =='Luminex'){CSF.vars <- c('LuminexTTau','LuminexPTau','LuminexAbeta42')}
	if(CSF.type =='Elisa'){CSF.vars <- c('ElisaTTau','ElisaPTau','ElisaAbeta42')}
	#CSF.vars <- c('TotalTau','PhosphorylatedTau','Beta.amyloid42')
	
	return(CSF.vars)

}

process.CSF <- function(CSF,CSF.type){
	# INPUTS:
	# CSF: dataframe with CSF testing data
	# CSF.type: Luminex or Elisa, type of CSF assay
	#
	# OUTPUTS:
	# CSF: processed CSF data

	CSF.vars <- get.CSF.vars(CSF.type)

	CSF <- CSF[!rowSums(is.na(CSF[,CSF.vars])),] # remove missing
	CSF <- CSF[order(CSF$INDDID),]

	CSF$CSFDate <- as.Date(CSF$CSFDate,format = '%m/%d/%Y')
	return(CSF)

}

extract.CSF.sample <- function(CSF,CSF.type,n.sample='first'){
	# INPUTS:
	# CSF: processed CSF data (output of process.CSF() above)
	# CSF.type: Luminex or Elisa, type of CSF assay
	# n.sample: which CSF sample to use, first or last
	#
	# OUTPUTS:
	# CSF.sample: dataframe of CSF features from specified sample (first or last)

	CSF.vars <- get.CSF.vars(CSF.type)

	CSF.by.pt <- lapply(sort(unique(CSF$INDDID)), function(id) # order tests by date w/in subjs
	  CSF[CSF$INDDID==id,c('CSFDate',CSF.vars)])
	CSF.by.pt <- lapply(CSF.by.pt, function(X) X[order(X$CSFDate),])
	#CSF.sample <- do.call('rbind',lapply(CSF.by.pt, function(X) colMeans(X)))
	if(n.sample == 'first'){
		CSF.sample <- do.call('rbind',lapply(CSF.by.pt, function(X) X[1,])) # just use first available sample
	} else if(n.sample == 'last'){
		CSF.sample <- do.call('rbind',lapply(CSF.by.pt, function(X) X[nrow(X),])) # just use first available sample
	}
	# for some patients, could compute a feature based on change in CSF proteins
	#CSF.diff <- do.call('rbind',lapply(CSF.by.pt, function(X) X[nrow(X),] - X[1,]))
	CSF.sample <- cbind(data.frame(INDDID=sort(unique(CSF$INDDID))),CSF.sample)
	CSF.sample$CSFDate <- as.Date(CSF.sample$CSFDate)
	return(CSF.sample)
}

addnormal <- function(df,CSF.type,data.date,n.sample='first'){
	# INPUTS:
	# df: dataframe containing INDDIDs and CSF
	# CSF.type: Luminex or Elisa, type of CSF assay
	# data.date: date of data download
	# n.sample: which CSF sample to use, first or last
	#
	# OUTPUTS:
	# df: dataframe with normal CSF samples added
	
	patients <- read.csv(file = paste0("data/INDD_Patients",data.date,".csv"),stringsAsFactors = FALSE)
	dx.global <- read.csv(file = paste0("data/INDD_GlobalDiagnosis",data.date,".csv"),stringsAsFactors = FALSE)
	CSF <- read.csv(paste0('data/INDD_CSF',CSF.type,data.date,'.csv'),stringsAsFactors = F)
	nl.INDDIDs <- dx.global$INDDID[dx.global$GlobalDx == 'Normal' & dx.global$NPDx ==''] # people whos global dx listed as normal and no conflicting neuropath dx
	nl.INDDIDs <- unique(c(nl.INDDIDs,patients$INDDID[grep('Unremarkable',patients$NPDx1)],patients$INDDID[grep('Normal',patients$NPDx1)]))
	nl.INDDIDs <- nl.INDDIDs[!(nl.INDDIDs %in% df$INDDID) & nl.INDDIDs %in% CSF$INDDID] # remove duplicates
	CSF <- CSF[CSF$INDDID %in% nl.INDDIDs,]
	CSF <- process.CSF(CSF,CSF.type)
	CSF <- extract.CSF.sample(CSF,CSF.type,n.sample=n.sample)
	CSF <- CSF[!CSF$INDDID %in% df$INDDID,] # remove duplicates
	df <- rbind(df,CSF)
	return(list(df=df,nl.INDDIDs=nl.INDDIDs))
}
