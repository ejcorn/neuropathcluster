# neuropathcluster

Code to reproduce all analysis in Cornblath et al. 2019 ("Defining and predicting transdiagnostic categories of neurodegenerative disease").

## Requirements:
  - MATLAB R2017a or later
  - R 3.3.3 or later. Requisite packages are listed in code/misc/packages.R
  - Brain Connectivity Toolbox, MATLAB version: https://sites.google.com/site/bctnet/

## Directory structure

Master branch contains 4 major folders:
  - code/ contains folders, each of which contain scripts specific to certain analyses, i.e. ‘code/clustering’ contains code that prepares for and performs clustering, while ‘code/assesscluster’ contains codes that handles clustering output and assesses clusters.
  - data/ contains csv files with the neuropathology data, CSF data, MoCA, and genotyping data. These data are available upon request. The subfolder data/img contains files for cerebral surface visualizations.

## Input specification

The file ‘pipeline.R’ is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the. At the top of ‘pipeline.R,’ one must select their own home directory (path to the main directory containing the above 4 folders) and path to MATLAB binary, in addition to the following variables:
  - missing.thrsh.r: sets maximum proportion of missing data for a given subject, above which subjects are excluded. We used 0.75.
  - missing.thrsh.c: sets maximum proportion of missing subjects for a given feature, above which features are excluded. We used 1.0, so that every feature was included.
  - extralab: an extra label to identify a particular run of the pipeline. This can be any string. 
  - nreps_gammasweep: number of partitions to generate at each value of gamma, used to evaluate partition stability as a function of gamma
  - gamma.opt: choose gamma value for analyses of disease clusters (Fig. 4-8). The suitability of various gamma values can be identified by viewing the generated plot that outputs to resultsdir/optimcluster/MeanzRandbyGamma0to3nreps1000reps.pdf. We found a plateau of partition stability from gamma = 1.5-1.8, and present results for gamma = 1.7 in the paper.
  - BCT.path: local path to Brain Connectivity Toolbox, MATLAB version (found here: https://sites.google.com/site/bctnet/)

## Questions, suggestions, comments?

Please contact Eli Cornblath (Eli.Cornblath@pennmedicine.upenn.edu) with any questions regarding this project.
