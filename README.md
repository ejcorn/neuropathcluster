# neuropathcluster

Code to reproduce all analysis in Cornblath et al. 2020, *Defining and predicting transdiagnostic categories of neurodegenerative disease*, https://www.nature.com/articles/s41551-020-0593-y.

## Requirements:
  - R 3.3.3 or later. Requisite packages are listed in `code/misc/packages.R`
  - If you want to use modularity maximization for clustering (not used in publication):
    - MATLAB R2017a or later
    - Brain Connectivity Toolbox, MATLAB version: https://sites.google.com/site/bctnet/

## Directory structure

Master branch contains 4 major folders:
  - `code/` contains folders, each of which contain scripts specific to certain analyses, i.e. ‘code/clustering’ contains code that prepares for and performs clustering, while ‘code/assesscluster’ contains codes that handles clustering output and assesses clusters.
  - `data/` should contain csv files with the neuropathology data, CSF data, MoCA, and genotyping data. These data are available upon request (see paper). The subfolder data/img contains files for cerebral surface visualizations.

## Input specification

The file `pipeline_polychor.R` is located in the main directory. This file will coordinate the sequential execution of all scripts within the `code/` folder, generating all the figures and analyses in the published manuscript and supplement. At the top of `pipeline_polychor.R`, one must select their own home directory (path to the main directory containing the above 4 folders) and path to MATLAB binary (not necessary unless you want to try modularity maximization), in addition to the following variables:
  - missing.thrsh.r: sets maximum proportion of missing data for a given subject, above which subjects are excluded. We used 0.75.
  - missing.thrsh.c: sets maximum proportion of missing subjects for a given feature, above which features are excluded. We used 1.0, so that every feature was included.
  - extralab: an extra label to identify a particular run of the pipeline. This can be any string. 
  - nreps_gammasweep: number of partitions to generate at each value of gamma, used to evaluate partition stability as a function of gamma
  - dist.met: distance metric for clustering. Here polychoric correlation given semi-quantitative data, but could use spearman or pearson.
  - gamma.opt: select number of clusters to analyze (determined here by silhouette and AUC analysis), or gamma value if using modularity maximization.
  - BCT.path: local path to Brain Connectivity Toolbox, MATLAB version (found here: https://sites.google.com/site/bctnet/). Not needed unless using modularity maximization


## Questions, suggestions, comments?

Please contact Eli Cornblath (Eli D0T Cornblath @t pennmedicine dºt upenn do† edu) with any questions regarding this project.
