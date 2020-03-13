# opdir is the overarching output directory specific to a data matrix constructed a particular way
# i.e. missing thresholds for rows and columns, any other small change that can be specified in 'extralab'
params$opdir <- paste('neuropathcluster_R',params$missing.thrsh.r,'C',params$missing.thrsh.c,params$extralab,'/',sep='')
dir.create(params$opdir,recursive = T)

# resultsdir contains results of analyses for a particular choice of gamma
params$resultsdir <- paste(params$opdir,'results_G',params$gamma.opt,'/',sep='')
dir.create(params$resultsdir,recursive = T)

# sourcedata.dir contains numerical data underlying all figures and tables
params$sourcedata.dir <- paste0(params$opdir,'SourceData/')
dir.create(params$sourcedata.dir,recursive = T)
