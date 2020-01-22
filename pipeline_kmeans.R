###########################
### Set input variables ###
###########################

rm(list=ls())
homedir <- '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/'
setwd(homedir)

# list of prespecified variables
params <- list(missing.thrsh.r=0.2,
               missing.thrsh.c=0,
               extralab='_122519',
               dist.met='spearman', # distance metric for modularity maximization
               gamma.opt=1.5,
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

#######################
### Preprocess data ###
#######################

# go from raw INDD csvs to base data for analysis
#source('code/preprocess/GenerateSample_v3.R')
#source('code/preprocess/demographics.R') # analyze demographics (age)

########################
### Cluster Patients ###
########################

#source('code/kmeans/kmeans.R')

#########################
### Assess clustering ###
#########################

# make sure clusters are in same order every time
source('code/assesscluster/CharacterizeLouvainClusters.R')
source('code/assesscluster/AgeMissingDataByCluster.R')

########################################
### Cognition, Genes, CSF by cluster ###
########################################

source('code/cogcsf/CogByCluster.R')
source('code/genes/ProcessAlleles.R')
source('code/genes/AlleleProportionsByCluster.R')
source('code/genes/AllelesByCluster.R')
source('code/cogcsf/CSFByCluster.R')
source('code/cogcsf/CSFByClusterOutliers.R') # ensure CSF analysis is robust to outliers

################################
### Exclude/isolate diseases ###
################################

dz.exc <- 'Alzheimer\'s disease'
exc.cl <- 2
n.dx <- 4
source('code/cogcsf/CogByClusterExcludeDisease.R')
source('code/cogcsf/CSFByClusterExcludeDisease.R')
source('code/cogcsf/CSFByClusterOutliersExcludeDisease.R') # ensure CSF analysis is robust to outliers
source('code/genes/AlleleProportionsByClusterExcludeDisease.R')
source('code/genes/AllelesByClusterExcludeDisease.R')

dz.exc <- 'PSP'
exc.cl <- 1
n.dx <- 4
source('code/genes/AllelesByClusterExcludeDisease.R')

dz.iso <- 'Alzheimer\'s disease'
n.dx <- 4
source('code/cogcsf/CogByClusterIsolateDisease.R')
source('code/cogcsf/CSFByClusterIsolateDisease.R')
source('code/genes/AllelesByClusterIsolateDisease.R')

# dz.excs <- c('Alzheimer\'s disease','PSP','Parkinson\'s disease')
# exc.cls <- c(2,1,4); names(exc.cls) <- dz.excs
# n.dx <- 4
# for(dz.exc in dz.excs){source('code/genes/AllelesByClusterExcludeDisease.R')}

##############################
### Predict disease labels ###
##############################

# GLM
extralabs <- c('CSFOnly')#,'GeneOnly','CSFGene')
for(extralab in extralabs){
  source('code/predictdisease/pd_prepdata.R')
  source('code/predictdisease/pd_traintestglm.R')
  source('code/predictdisease/pd_plotperfglm.R')
  source('code/predictdisease/pd_plotfeatureweightsglm.R')
}

# GLM -- comparing to overlapping traditional diagnoses
extralabs <- c('CSFOnly')
for(extralab in extralabs){
  source('code/predictdisease/pd_prepdata_alldx.R') # construct data frame allowing for overlap in traditional dx
  source('code/predictdisease/pd_traintestglm.R')
  source('code/predictdisease/pd_plotperfglm.R')
  source('code/predictdisease/pd_plotfeatureweightsglm.R')
}

# random forest
extralabs <- c('CSFGene')#,'CSFOnly','GeneOnly')
for(extralab in extralabs){
  source('code/predictdisease/pd_prepdata.R')
  source('code/predictdisease/pd_traintestrf.R')
  source('code/predictdisease/pd_plotperfrf.R')
  source('code/predictdisease/pd_plotfeatureweightsrf.R')
}

# visualize disease labels in 2D and 3D CSF space
source('code/predictdisease/CSFspace.R')

##############################################
### Predict disease labels exclude disease ###
##############################################

extralab <- c('CSFOnly')
dz.exc <- 'Alzheimer\'s disease'
exc.cl <- 2
n.dx <- 4
source('code/predictdisease/pd_traintestExcludeDisease.R')
extralab <- paste(extralab,'Exclude',dz.exc,sep='')
source('code/predictdisease/pd_plotperfExcludeDisease.R')
