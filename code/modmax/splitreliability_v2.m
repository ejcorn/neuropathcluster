clearvars -except nreps gamma BCT_path sampfrac homedir opdir; close all; clc

% BCT_path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT';
% homedir = '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/';
% opdir = 'neuropathcluster_R0.75C1allpts_final/';
% nreps = 1000;
%%
savedir = [homedir,opdir,'optimcluster/']; 
mkdir(savedir);
addpath(genpath(BCT_path)); % add BCT functions to path
addpath(genpath([homedir,'code/matlab_functions'])) % add other ancillary matlab functions to path
load([homedir,opdir,'processed/pathDataForClustering.mat']);
X = X+.1; % shift back -- had to shift by non-integer to make R save as double for matlab

%% set parameters for subsampling
nobs = size(X,1);
nperms = 1000;
sampsize = round(sampfrac*nobs);
nreps = 50; % # of times to repeat community detection

partitions = nan(sampsize,nperms);
samples = nan(sampsize,nperms);
%%
parpool(4);

%% perform louvain clustering on subsamples
parfor i = 1:nperms
    disp(['Perm ',num2str(i)])
    sample = randperm(nobs,sampsize);
    % demean each measurement
    X_samp = (X(sample,:) - nanmean(X(sample,:),1))';    
    [W,p] = corr(X_samp,'type','spearman','rows','pairwise');        
    W = W.*(p < 0.05); % using p < 0.05 to maintain connected graph     
    C = zeros(sampsize,nreps);
    for R = 1:nreps
        Mi = community_louvain(W,gamma,[],'negative_sym');
        C(:,R) = Mi; 
    end
    % synchronize partitions
    C = multislice_pair_labeling(C);
    [C,~,~] = consensus_similarity(C');
    % store consensus partition
    partitions(:,i) = C;
    % store sample indices to compute centroids in R
    samples(:,i) = sample;
end

% synchronize partitions
partitions = multislice_pair_labeling(partitions);
%%
fname = [savedir,'subsamplePartitionsSF',num2str(sampfrac),'.mat'];
save(fname,'partitions','samples');

%%
delete(gcp);