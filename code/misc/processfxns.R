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

addnormal.CSF <- function(df,CSF.type,data.date,n.sample='first'){
	# INPUTS:
	# df: dataframe containing INDDIDs and CSF
	# CSF.type: Luminex or Elisa, type of CSF assay
	# data.date: date of data download
	# n.sample: which CSF sample to use, first or last
	#
	# OUTPUTS:
	# CSF: dataframe with normal CSF samples ONLY
	
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
	#df <- rbind(df,CSF) # concatenate CSF --- commented out so I can be more flexible in script
	return(list(CSF=CSF,nl.INDDIDs=nl.INDDIDs)) # return(list(df=df,nl.INDDIDs=nl.INDDIDs))
}


process.Gene <- function(Genes,GOIs,Alleles,remove.wt=FALSE){
	# INPUTS:
	# Genes: data frame with genetic data from INDD
	# GOIs: genes of interest to select
	# Alleles: list of lists for each gene containing allele names, assuming genotype listed as "A1/A2" and first allele listed in wild type
	# remove.wt: logical indicating whether to remove wild type allele from table
	#
	# OUTPUTS:
	# Allele.Tables: list of dataframes for each gene containing counts of each allele
	# Genes: Genes data frame with exclusion for missing data applied
	
	# always do this for APOE

	Genes$APOE[Genes$APOE == 'E2/E3  '] <- 'E2/E3'
	Genes$APOE[Genes$APOE == 'E3/E4  '] <- 'E2/E4'
	missing.mask <- rowSums(Genes[,GOIs] != '') == length(GOIs)
	Genes <- Genes[missing.mask,]  # remove missing
	Genes.df <- Genes[,GOIs] # select genes of interest
	
	Genotypes <- lapply(Genes.df, function(X) unique(X))
	# for each gene, get an Allele-by-Genotype table counting # of each allele per genotype
	n.alleles.per.genotype <- list()
	Allele.Tables <- list()
	for(g.i in names(Genotypes)){
	  n.alleles.per.genotype[[g.i]] <- 
	    sapply(Genotypes[[g.i]], function(gt.i)
	      sapply(Alleles[[g.i]], function(A) str_count(pattern = A,string = gt.i)))
	  Allele.Tables[[g.i]] <- t(sapply(Genes.df[,g.i], function(gt.i) n.alleles.per.genotype[[g.i]][,gt.i]))
	  rownames(Allele.Tables[[g.i]]) <- NULL
	 	# remove wild type allele [-1] to avoid redundancy, i.e. if you're not H1 you are definitely H2 so 
  		# including both as predictors is unnecessary
  		if(remove.wt){Allele.Tables[[g.i]] <- Allele.Tables[[g.i]][,Alleles[[g.i]][-1],drop=FALSE]}
  		Allele.Tables[[g.i]] <- as.data.frame(Allele.Tables[[g.i]])
	}
	return(list(Allele.Tables=Allele.Tables,Genes=Genes))

}

addnormal.Genes <- function(df,GOIs,Alleles,data.date,remove.wt=TRUE){
	# INPUTS:
	# df: dataframe containing INDDIDs and alleles
	# GOIs: genes of interest to select
	# Alleles: list of lists for each gene containing allele names, assuming genotype listed as "A1/A2" and first allele listed in wild type
	# data.date: date of data download
	# remove.wt: passed through to process.Gene() and specifies whether you want to remove the wildtype allele
	# 
	# OUTPUTS:
	# Allele.Tables: dataframe with Alleles from normall subjects ONLY
	# nl.INDDIDs: inddids of normal subjects
	
	patients <- read.csv(file = paste0("data/INDD_Patients",data.date,".csv"),stringsAsFactors = FALSE)
	dx.global <- read.csv(file = paste0("data/INDD_GlobalDiagnosis",data.date,".csv"),stringsAsFactors = FALSE)
	Genes <- read.csv(paste0('data/INDD_Genetics',data.date,'.csv'),stringsAsFactors = F)
	nl.INDDIDs <- dx.global$INDDID[dx.global$GlobalDx == 'Normal' & dx.global$NPDx ==''] # people whos global dx listed as normal and no conflicting neuropath dx
	nl.INDDIDs <- unique(c(nl.INDDIDs,patients$INDDID[grep('Unremarkable',patients$NPDx1)],patients$INDDID[grep('Normal',patients$NPDx1)]))
	nl.INDDIDs <- nl.INDDIDs[!(nl.INDDIDs %in% df$INDDID) & nl.INDDIDs %in% Genes$INDDID] # remove duplicates
	Genes <- Genes[Genes$INDDID %in% nl.INDDIDs,]
	list[Allele.Tables,Genes] <- process.Gene(Genes,GOIs,Alleles,remove.wt)	
	names(Allele.Tables) <- NULL # this is to prevent R from adding back in list item names to df columns
	# add IDs back and concatenate across genes
	Allele.Tables <- cbind(data.frame(INDDID=Genes$INDDID),do.call('cbind',Allele.Tables))
	Allele.Tables <- Allele.Tables[!Allele.Tables$INDDID %in% df$INDDID,] # remove duplicates
	nl.INDDIDs <- Allele.Tables$INDDID
	#df <- rbind(df,Genes)
	return(list(Allele.Tables=Allele.Tables,nl.INDDIDs=nl.INDDIDs)) #return(list(df=df,nl.INDDIDs=nl.INDDIDs))
}
