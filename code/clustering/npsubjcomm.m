clear all; close all; clc
%%
addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/BCT'))
addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/NCT'))
addpath(genpath('~/Dropbox/Neurodegeneration/PathCogClinDx/code/'))
load ~/Dropbox/Neurodegeneration/PathCogClinDx/data/pathDataForClustering.mat
X = X+.1; % shift back -- had to shift by non-integer to make R save as double for matlab
X = X - nanmean(X,1); 

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

% set remaining nans to 0
% this isn't a bad assumption because the pathology is baseline very
% dissimilar, so you shouldn't draw an edge of any weight

%W(isnan(W)) = 0;

%% ensure that diagonal is equal to 1
W(logical(eye(length(W)))) = 1;

%% 
p_fdr = triu(p) .* ~eye(length(p)); p_fdr = p_fdr(p_fdr > 0); % isolate all indiv. test p-values
[~,thresh] = fdr_bh(p_fdr); % compute fdr threshold
thrsh = 0.05 / (nobs*(nobs-1)/2); % bonferroni corrected threshold for signifance
%W = W.*(p < thrsh) % gives disconnected graph
W = W.*(p < 0.05); % using p < 0.05 to maintain connected graph
[~,~,r]= dmperm(W) % confirm that graph is connected

%%
parpool(3)

%% use louvain consensus to group patients into clusters
gamma = 1.4;
% M0 contains initial community assignments based on amount of each
% type of pathology
nreps = 100;
C = zeros(nreps,nobs);
parfor i = 1:nreps
    disp(['Rep ',num2str(i)])
    [Mi,Q1] = community_louvain(W,gamma,[],'negative_sym');
    C(i,:) = Mi;
end

C = multislice_pair_labeling(C');
[M,~,~] = consensus_similarity(C');
M=M';
%%
[~,I] = sort(M,'ascend');
figure;imagesc(W(I,I)); colormap('redblue'); colorbar
save(['~/Dropbox/Neurodegeneration/PathCogClinDx/data/subjectClusterLouvain','.mat'],'M','W','DisconnectedSubjects')
max(M)

%%
load('~/Dropbox/Neurodegeneration/PathCogClinDx/data/subjectClusterLouvainold.mat','M')
M(DisconnectedSubjects) = [];
[~,I] = sort(M,'ascend');
figure;imagesc(W(I,I)); colormap('redblue'); colorbar

%% k-means
M = kmeansSpearman(X,4,'Replicates',20,'Distance','correlation');
[~,I] = sort(M);
figure;imagesc(W(I,I));
save('~/Dropbox/Neurodegeneration/PathCogClinDx/data/subjectClusterkMeans.mat','M')
