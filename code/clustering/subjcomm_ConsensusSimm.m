clearvars -except nreps BCT_path homedir opdir g; close all; clc

% BCT_path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT';
% homedir = '~/Dropbox/Neurodegeneration/PathCogClinDx/';
% opdir = 'neuropathcluster_R0.75C1allpts/';
% nreps = 1000;

addpath(genpath(BCT_path)); % add BCT functions to path
addpath(genpath([homedir,'code/matlab_functions'])) % add other ancillary matlab functions to path

savedir = [homedir,opdir,'optimcluster/']; 
mkdir(savedir);

gamma_rng = 0:0.1:3;

%% extract partition stability by gamma and find consensus partition
load(fullfile(savedir,'subjectCorrMat.mat'),'DisconnectedSubjects')
load([savedir,'LouvainSubjNPSweepGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'.mat']);

nobs = size(partitions{g},1);
disp('calculating z-rand scores')
zr_scores = zeros(nreps*(nreps-1)/2,1);
sim_mats = zeros(nreps,nreps,1);
partition_gamma = zeros(nobs,1);

disp(['Gamma = ',num2str(gamma_rng(g))]);
tic
[M,~,sim_mat] = consensus_similarity(partitions{g}'); %'
toc
partition_gamma(:,1) = M;
sim_mats(:,:,1) = sim_mat;
zr_unique = triu(sim_mat).*~eye(length(sim_mat));
zr_scores(:,1) = zr_unique(zr_unique ~= 0); % extract only upper triangle b/c matrix is symmetric

save([savedir,'ConsensusPartitionGamma',num2str(gamma_rng(g)),'NReps',num2str(nreps),'.mat'],'partition_gamma','gamma_rng','DisconnectedSubjects');
save(fullfile(savedir,['ZRandScoresGamma',num2str(gamma_rng(g)),'NReps',num2str(nreps),'.mat']),'zr_scores');
