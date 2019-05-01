clearvars -except nreps BCT_path homedir opdir; close all; clc

% BCT_path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT';
% homedir = '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/';
% opdir = 'neuropathcluster_R0.75C1allpts/';
% nreps = 1000;
%%
savedir = [homedir,opdir,'optimcluster/']; 
mkdir(savedir);
addpath(genpath(BCT_path)); % add BCT functions to path
addpath(genpath([homedir,'code/matlab_functions'])) % add other ancillary matlab functions to path
load([homedir,opdir,'processed/pathDataForClustering.mat']);
X = X+.1; % shift back -- had to shift by non-integer to make R save as double for matlab
X = X - nanmean(X,1); % demean each feature

nobs = size(X,1);

%% load correlation matrix
load(fullfile(savedir,'subjectCorrMat.mat'))
%%
parpool(4)

%% use louvain consensus to group patients into clusters
gamma_rng = 0:0.1:3;
n_gamma = length(gamma_rng);

partitions = cell(n_gamma,1);
modularity = cell(n_gamma,1);
tic
parfor g = 1:n_gamma  % loop through gamma values
    disp(['Gamma = ',num2str(gamma_rng(g))]);
    gamma = gamma_rng(g);
    C = zeros(nreps,nobs);
    Q = zeros(nreps,1);
    for i = 1:nreps
        disp(['Rep ',num2str(i)])
        [Mi,Qi] = community_louvain(W,gamma,[],'negative_asym');
        C(i,:) = Mi;
        Q(i) = Qi;
    end
    C = multislice_pair_labeling(C'); % synchronize
    partitions{g} = C; % store partitions
    modularity{g} = Q; % store modularity values
end
toc

save([savedir,'LouvainSubjNPSweepGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'.mat'],'partitions','modularity')
%% extract partition stability by gamma and find consensus partition

disp('calculating z-rand scores')
zr_scores = zeros(nreps*(nreps-1)/2,n_gamma);
sim_mats = zeros(nreps,nreps,n_gamma);
partitions_by_gamma = zeros(nobs,n_gamma);

parfor g = 1:n_gamma
    disp(['Gamma = ',num2str(gamma_rng(g))]);
    [M,~,sim_mat] = consensus_similarity(partitions{g}'); %'
    partitions_by_gamma(:,g) = M;
    sim_mats(:,:,g) = sim_mat;
    zr_unique = triu(sim_mat).*~eye(length(sim_mat));
    zr_scores(:,g) = zr_unique(zr_unique ~= 0); % extract only upper triangle b/c matrix is symmetric
end

partitions_by_gamma = multislice_pair_labeling(partitions_by_gamma);
save(fullfile(savedir,['subjectClusterLouvainPartitionsByGammaNReps',num2str(nreps),'.mat']),'partitions_by_gamma','gamma_rng','DisconnectedSubjects');
save(fullfile(savedir,['ZRandScoresLouvainSubjNPSweepGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'.mat']),'zr_scores');

%% visualize distribution of z-scored rand index at each level of gamma

% remove in final code
f=figure; [sp1,sp2] = subplot_ind2(n_gamma);
for g = 1:n_gamma
    subplot(sp1,sp2,g); hold on; axis square
    histogram(zr_scores(:,g));
    title(['\gamma = ',num2str(gamma_rng(g))]);
end
saveas(f,fullfile(savedir,['zRandHistogramsbyGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'reps.pdf']),'pdf');
save(fullfile(savedir,['ZRandScoresLouvainSubjNPSweepGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'.mat']),'zr_scores');
%% plot partition stability by gamma

zr_mean_by_gamma = mean(zr_scores,1);
[~,zrBestIdx] = max(zr_mean_by_gamma);
zrBestGamma = gamma_rng(zrBestIdx);
f=figure; hold on;
plot(gamma_rng,zr_mean_by_gamma,'Color','blue');
xticks(gamma_rng(1:2:end)); xticklabels(gamma_rng(1:2:end));
ylabel('Mean Z-Scored Rand Index'); xlabel('\gamma','FontWeight','bold');
mi = ylim(gca);
line([zrBestGamma zrBestGamma],[mi(1) max(zr_mean_by_gamma)],'LineStyle','--','Color','red');
prettifyEJC;
f.PaperUnits = 'inches';
f.PaperSize = [3 2];
f.PaperPosition = [0 0 3 2];
saveas(f,fullfile(savedir,['MeanzRandbyGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'reps.pdf']),'pdf');

zr_median_by_gamma = median(zr_scores,1);
[~,zrBestIdxMed] = max(zr_mean_by_gamma);
zrBestGammaMed = gamma_rng(zrBestIdxMed);
f=figure; hold on;
plot(gamma_rng,zr_median_by_gamma,'Color','blue');
xticks(gamma_rng(1:2:end)); xticklabels(gamma_rng(1:2:end));
ylabel('Median Z-Scored Rand Index'); xlabel('\gamma','FontWeight','bold');
mi = ylim(gca);
line([zrBestGammaMed zrBestGammaMed],[mi(1) max(zr_median_by_gamma)],'LineStyle','--','Color','red');
prettifyEJC;
f.PaperUnits = 'inches';
f.PaperSize = [3 2];
f.PaperPosition = [0 0 3 2];
saveas(f,fullfile(savedir,['MedianzRandbyGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'reps.pdf']),'pdf');


%% show number and size of clusters for each gamma value

k_max = max(max(partitions_by_gamma)); % find maximum number of clusters across all gamma
k_max = 10; % just focus on major clusters
cluster_size_by_gamma = zeros(k_max,n_gamma);
for g = 1:n_gamma
    P = partitions_by_gamma(:,g);
    for k_i = 1:k_max
        % count # of patients in each possible cluster
        cluster_size_by_gamma(k_i,g) = sum(P == k_i) / size(partitions_by_gamma,1); 
    end
end

f=figure; imagesc(flipud(cluster_size_by_gamma))
xticks(1:2:n_gamma); xticklabels(gamma_rng(1:2:end));
yticks(1:k_max); yticklabels(fliplr(1:k_max));
ylabel('Cluster'); xlabel('\gamma','FontWeight','bold');
h=colorbar; ylabel(h,'Proportion of subjects');
set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [4 2];
f.PaperPosition = [0 0 4 2];
saveas(f,fullfile(savedir,['NumberPerClusterbyGamma',num2str(min(gamma_rng)),'to',num2str(max(gamma_rng)),'nreps',num2str(nreps),'reps.pdf']),'pdf');

%%
delete(gcp)