###########################
### Set input variables ###
###########################

rm(list=ls())
homedir <- '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/'
setwd(homedir)

# list of prespecified variables
params <- list(missing.thrsh.r=0.2,
               missing.thrsh.c=0,
               extralab='_010320',
               dist.met='polychoric', # distance metric for modularity maximization
               gamma.opt=6,
               nreps_gammasweep=1000, # set number of reps for gamma sweep
               homedir=homedir,
               matlab.path='/Applications/MATLAB_R2019a.app/bin/matlab', # path to matlab binary
               BCT.path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT') # path to folder containing brain connectivity toolbox scripts (https://sites.google.com/site/bctnet/)

#####################
### Load packages ###
#####################

source('code/misc/packages.R')

#################################
### Define output directories ###
#################################

source('code/misc/directories.R')

# ~~~~~~~~~~ Run code to this point ~~~~~~~~~~~ #
# ~~~~~ then run any script independently ~~~~~ #

#######################
### Preprocess data ###
#######################

# go from raw INDD csvs to base data for analysis
source('code/preprocess/GenerateSample_v3.R') # Fig S1a-d

########################
### Cluster Patients ###
########################

source('code/kmeans/Imputation.R')
source('code/clustering/PrepDataCluster.R') # compute polychoric correlation matrix... time consuming
source('code/kmeans/pam.R')

# modularity maximization
if(params$gamma.opt < 2){
  source('code/clustering/runGammaSweepMATLAB_R.R')
  sampfrac <- 0.5
  source('code/clustering/runSplitReliabilityMATLAB_R.R')
  source('code/assesscluster/ProcessCluster.R') # make sure clusters are in same order every time
}

source('code/kmeans/PAMBootstrapReliabilityScrambleNorm.R') # bootstrap clustering
source('code/kmeans/PartitionSimilarity_PAMBootstrapReliabilityScrambleNorm.R') # compute similarity between bootstrapped partitions
source('code/kmeans/Plot_PAMBootstrapReliabilityScrambleNorm.R') # Figure S6a-h -- these three scripts are serially dependent

######################################################################
### Perform expl. factor analysis on polychoric correlation matrix ###
######################################################################

source('code/kmeans/EFA.R')
source('code/kmeans/PolychorFeatures.R')
source('code/kmeans/AnalyzeEFA.R') # Fig S9
source('code/kmeans/umap_cluster.R')

#########################
### Assess clustering ###
#########################

source('code/assesscluster/CharacterizeLouvainClustersDiseaseStages.R') # Fig 2d-f
source('code/assesscluster/BreakdownOtherDiagnoses.R') # Fig S2
source('code/assesscluster/ClinicalDxByCluster.R')
source('code/assesscluster/SubjectCorrMatLouvain.R') # Fig 2a-c
source('code/assesscluster/AgeMissingDataByCluster.R') # Fig S3a-g, Fig S5a-c, Fig S6a
source('code/assesscluster/EOADByCluster.R') # Fig S6b
source('code/assesscluster/ADLBDClusterivsClusterj.R') # Fig 3a-f
source('code/assesscluster/ClustersPathSpace.R')
source('code/assesscluster/PathMatrixClustered.R') # Fig S9a

#######################################
### Make surface plots of centroids ###
#######################################

source('code/plot_brains/prep_centroid_plots.R') # prep data for runPlotBrainsMATLAB_R.R
source('code/plot_brains/runPlotBrainsMATLAB_R.R') # Fig S8a

########################################
### Cognition, Genes, CSF by cluster ###
########################################

source('code/cogcsf/CogByCluster.R') # Fig 4a, Fig S4a
source('code/genes/ProcessAlleles.R') # process for next 2 scripts
source('code/genes/GenotypeProportionsByCluster.R') # Fig 7a-b
source('code/genes/AllelesByCluster.R')  # Fig 7c-g
source('code/cogcsf/CSFByCluster.R') # Fig 6a
source('code/cogcsf/CSFVsOnset.R') # Fig S4b
source('code/cogcsf/MMSEByCluster.R') # Fig S4c-d

source('code/cogcsf/predictMoCAMMSE.R') # Fig 5a-c computation
source('code/cogcsf/plot_predictMoCAMMSE.R') # Fig 5a-c plot

################################
### Exclude/isolate diseases ###
################################

dz.exc <- 'Alzheimer\'s disease'
exc.cl <- c(2,4,5)
n.dx <- 5
source('code/cogcsf/CogByClusterExcludeDisease.R') # Fig S10a
source('code/cogcsf/CSFByClusterExcludeDisease.R') # Fig S11a
source('code/genes/GenotypeProportionsByClusterExcludeDisease.R') # Fig S12a-b
source('code/genes/AllelesByClusterExcludeDisease.R') # Fig S12c-e

dz.exc <- 'PSP'
exc.cl <- 1
n.dx <- 5
source('code/genes/AllelesByClusterExcludeDisease.R') # Fig S13a-b

dz.iso <- 'Alzheimer\'s disease'
n.dx <- 5
source('code/cogcsf/CogByClusterIsolateDisease.R') # Fig S10b
source('code/cogcsf/CSFByClusterIsolateDisease.R') # Fig S11b

##############################
### Predict disease labels ###
##############################

# GLM -- comparing to overlapping traditional diagnoses
extralabs <- c('CSFOnlyAddNormalMMSE','CSFGeneAddNormalMMSE','GeneOnlyAddNormalMMSE') # Fig 6, Fig S14, Fig S16
for(extralab in extralabs){
  source('code/predictdisease/pd_prepdata_alldx.R') # construct data frame allowing for overlap in traditional dx
  source('code/predictdisease/pd_traintestglm.R')
  source('code/predictdisease/pd_plotperfglm.R')
  source('code/predictdisease/pd_plotfeatureweightsglm.R')
}

# random forest
extralabs <- c('CSFGeneAddNormalMMSE') # Fig S15
for(extralab in extralabs){
  source('code/predictdisease/pd_prepdata_alldx.R')
  source('code/predictdisease/pd_traintestrf.R')
  source('code/predictdisease/pd_plotperfrf.R')
  source('code/predictdisease/pd_plotfeatureweightsrf.R')
}

# visualize disease labels in 2D and 3D CSF space
source('code/predictdisease/CSFspace.R')

##############################################
### Predict disease labels exclude disease ###
##############################################

extralabs <- c('CSFOnlyAddNormalMMSE','CSFGeneAddNormalMMSE') # Fig S17a-b, Fig S17c-d
dz.exc <- 'Alzheimer\'s disease'
exc.cl <- c(2,4,5)
n.dx <- 5
for(extralab in extralabs){
  source('code/predictdisease/pd_traintestExcludeDisease.R')
  extralab <- paste(extralab,'Exclude',dz.exc,sep='')
  source('code/predictdisease/pd_plotperfExcludeDisease.R')
}

#############################################################
### Predict disease labels using weighted class balancing ###
#############################################################

extralab <- 'CSFOnlyAddNormalMMSE'
source('code/predictdisease/pd_traintestglm_WCB.R')
source('code/predictdisease/pd_plotperfglm_WCB.R')

extralab <- 'CSFGeneAddNormalMMSE'
source('code/predictdisease/pd_traintestrf_WCB.R')
source('code/predictdisease/pd_plotperfrf_WCB.R')