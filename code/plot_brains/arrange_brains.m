clearvars -except k BCT_path homedir opdir resultsdir; close all; clc

k = 4;
BCT_path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT';
homedir = '~/Dropbox/Neurodegeneration/PathCogClinDx/';
resultsdir = 'neuropathcluster_R0.75C0.75/results_G2.1/';
%% Define constant variables

PathItems_Type = {'NeuronLoss','Gliosis','Angiopathy','Ubiquitin','Thio','TDP43','Tau','Syn','Antibody'};

%%
cd(homedir);
savedir = fullfile(homedir,resultsdir,'plot_brains'); 
addpath(genpath(BCT_path)); % add BCT functions to path
addpath(genpath([homedir,'code/plot_brains']))
addpath(genpath([homedir,'code/matlab_functions'])) % add other ancillary matlab functions to path

%% plot brain heatmaps for each cluster and each path type
% inputs to plotting functions
for k_i = 2
    savedir_k = fullfile(savedir,'brains',['Cluster',num2str(k_i)]);    
    for p_type = PathItems_Type %([5,6])
        % load and shrink subcortex/brainstem/cerebellum
        fname = fullfile(savedir_k,['SUB_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        sub_cdata = imread(fname);
        sub_cdata = crop_brain_image(sub_cdata);
        sub_cdata = imresize(sub_cdata,0.5);
        
        % load cortical surface
        fname = fullfile(savedir_k,['CRTX_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        crtx_cdata = imread(fname);
        crtx_cdata = crop_brain_image(crtx_cdata);
        
        % concatenate, title and save
        concat_img = img_vertcat_whitepad(crtx_cdata,sub_cdata);                
        f=figure;
        imshow(concat_img); title(p_type);
        f.PaperUnits = 'inches';
        f.PaperSize = [1 2];
        f.PaperPosition = [0 0 1 2];
        set(0,'CurrentFigure',f); set(gca,'FontSize',8);
        fname = fullfile(savedir_k,['FINAL_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        print(fname,'-dpng','-r400');
    end
end

