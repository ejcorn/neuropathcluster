clearvars -except nreps BCT_path homedir opdir; close all; clc

% BCT_path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT';
% homedir = '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/';
% opdir = 'neuropathcluster_R0.75C1allpts/';
% nreps = 5;
%%
savedir = [homedir,opdir,'optimcluster/']; 
mkdir(savedir);
addpath(genpath(BCT_path)); % add BCT functions to path
addpath(genpath([homedir,'code/matlab_functions'])) % add other ancillary matlab functions to path
load([homedir,opdir,'processed/pathDataForClustering.mat']);
X = X+.1; % shift back -- had to shift by non-integer to make R save as double for matlab
X = X - nanmean(X,1); % demean each feature

%% compute correlation matrix
nobs = size(X,1);

[W_orig,p] = corr(X','type','spearman','rows','pairwise');

%% deal with nans in correlation matrix
% nans come from subject pairs where >0 subjects pairwise available features 
% are uniform, i.e. all 1's usually
W = W_orig;
% remove disconnected nodes, if any

DisconnectedSubjects = find(sum(isnan(W),1) == length(W));

W(DisconnectedSubjects,:) = []; p(DisconnectedSubjects,:) = [];
W(:,DisconnectedSubjects) = []; p(:,DisconnectedSubjects) = [];

%% ensure that diagonal is equal to 1
W(logical(eye(length(W)))) = 1;

%% threshold W
%p_fdr = triu(p) .* ~eye(length(p)); p_fdr = p_fdr(p_fdr > 0); % isolate all indiv. test p-values
%[~,thresh] = fdr_bh(p_fdr); % compute fdr threshold
%thrsh = 0.05 / (nobs*(nobs-1)/2); % bonferroni corrected threshold for signifance
%W = W.*(p < thrsh) % gives disconnected graph

%W = W.*(p < 0.05); % using p < 0.05 to maintain connected graph
% if not enough features available for a meaningful correlation
% between 2 patients you get a NaN. Set it to 0 to allow comm detection
%W(isnan(W)) = 0; 
[~,~,r]= dmperm(W) % confirm that graph is connected

save(fullfile(savedir,'subjectCorrMat.mat'),'W','DisconnectedSubjects')