clearvars -except k BCT_path homedir opdir resultsdir; close all; clc

% k = 4;
% BCT_path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT';
% homedir = '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/';
% resultsdir = 'neuropathcluster_R0.75C1allpts_final/results_G1.7/';

%% Define constant variables

PathItems_Type = {'NeuronLoss','Gliosis','Angiopathy','Thio','TDP43','Tau','Syn'};

%%
cd(homedir);
savedir = fullfile(homedir,resultsdir,'plot_brains'); 
addpath(genpath(BCT_path)); % add BCT functions to path
addpath(genpath([homedir,'code/plot_brains']))
addpath(genpath([homedir,'code/matlab_functions'])) % add other ancillary matlab functions to path

%% plot brain heatmaps for each cluster and each path type
% inputs to plotting functions
set(0,'DefaultFigureVisible','off'); % don't show me 70+ figs
cmap = 'plasma'; % colormap
min_col = [199 185 218] / 255; grey_n = [0.75 0.75 0.75];
max_col = [45 45 129]/255;
cmap = bichrome_cmap(min_col,max_col);

alpha = 1; % alpha of brains
res = ['-r',num2str(600)];  % image resolution
for k_i = 1:k
    savedir_k = fullfile(savedir,'brains',['Cluster',num2str(k_i)]);
    mkdir(savedir_k); % separate folder for each cluster
    for p_type = PathItems_Type
        fname = ['Cluster ',num2str(k_i),'_',char(p_type),'.mat'];
        x = load(fullfile(savedir,'vals',fname));
        C_region_plot_t = x.C_region_plot_t;
        mask = ~isnan(C_region_plot_t(:,1));
        C_region_plot_t = C_region_plot_t(mask,:); % remove absent regions
        
        % plot subcortex, brainstem, cerebellum
        [f_ant,f_lat] = plot_subcortvol3(C_region_plot_t(:,2),C_region_plot_t(:,1),cmap);
        f_ant = figure_resize(f_ant,'inches',[1 1]);
        f_lat = figure_resize(f_lat,'inches',[1 1]);
        set(0, 'currentfigure', f_ant);
        print(fullfile(savedir_k,['ANT_SUB_Cluster',num2str(k_i),'_',char(p_type),'.png']),'-dpng',res);        
        close(f_ant)
        set(0, 'currentfigure', f_lat);
        print(fullfile(savedir_k,['LAT_SUB_Cluster',num2str(k_i),'_',char(p_type),'.png']),'-dpng',res);        
        close(f_lat)
        
        % plot cortical surface
        [f_lat,f_med] = plot_cortsurf_CNDRautopsy(C_region_plot_t(:,2),C_region_plot_t(:,1),cmap,alpha);
        f_lat = figure_resize(f_lat,'inches',[1 1]);        
        set(0, 'currentfigure', f_lat);
        print(fullfile(savedir_k,['LAT_CRTX_Cluster',num2str(k_i),'_',char(p_type),'.png']),'-dpng',res);        
        close(f_lat)
        f_med = figure_resize(f_med,'inches',[1 1]);
        set(0, 'currentfigure', f_med);
        print(fullfile(savedir_k,['MED_CRTX_Cluster',num2str(k_i),'_',char(p_type),'.png']),'-dpng',res);                
        close(f_med)
    end
end

%% combine cortex and subcortex into one plot

crtx_marg = 0.02; % percent margin for each cortical view
sub_marg = 0.02; % percent margin each subcortical view
res = ['-r',num2str(600)];  % image resolution
for k_i = 1:k
    savedir_k = fullfile(savedir,'brains',['Cluster',num2str(k_i)]);    
    for p_type = PathItems_Type
        % load and shrink subcortex/brainstem/cerebellum
        fname = fullfile(savedir_k,['ANT_SUB_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        ant_sub_cdata = imread(fname);
        ant_sub_cdata = crop_brain_image(ant_sub_cdata,sub_marg);
        ant_sub_cdata = imresize(ant_sub_cdata,0.6);
        fname = fullfile(savedir_k,['LAT_SUB_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        lat_sub_cdata = imread(fname);
        lat_sub_cdata = crop_brain_image(lat_sub_cdata,sub_marg);
        % lat image is smaller. scale to match downsize ant image
        ant_lat_scl = size(ant_sub_cdata,1)/size(lat_sub_cdata,1);
        lat_sub_cdata = imresize(lat_sub_cdata,ant_lat_scl);
                
        % load cortical surface
        fname = fullfile(savedir_k,['LAT_CRTX_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        lat_crtx_cdata = imread(fname);
        lat_crtx_cdata = crop_brain_image(lat_crtx_cdata,crtx_marg);
        fname = fullfile(savedir_k,['MED_CRTX_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        med_crtx_cdata = imread(fname);
        med_crtx_cdata = crop_brain_image(med_crtx_cdata,crtx_marg);
        
        % concatenate, crop again, title and save
        sub_cdata = [ant_sub_cdata,lat_sub_cdata];
        crtx_cdata = [lat_crtx_cdata,med_crtx_cdata];
        concat_img = img_vertcat_whitepad(crtx_cdata,sub_cdata);                
        concat_img = crop_brain_image(concat_img,0.01);
        f=figure;
        imshow(concat_img); %title(['Cluster ',num2str(k_i)]);
        f.PaperUnits = 'inches';
        f.PaperSize = [1.75 1];
        f.PaperPosition = [0 0 1.75 1];
        set(0,'CurrentFigure',f); set(gca,'FontSize',8);
        fname = fullfile(savedir_k,['FINAL_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        print(fname,'-dpng',res);
        close(f)
    end
end

%% 
PathItems_Type = {'Syn','Tau','TDP43','Thio','NeuronLoss','Angiopathy'};
n_types = length(PathItems_Type);
fig_cell = cell(n_types,k);
for k_i = 1:k
    savedir_k = fullfile(savedir,'brains',['Cluster',num2str(k_i)]);    
    im_cell = cell(n_types,1);
    for j = 1:n_types
        p_type = PathItems_Type{j};
        fname = fullfile(savedir_k,['FINAL_Cluster',num2str(k_i),'_',char(p_type),'.png']);
        im_cell{j} = crop_brain_image(imread(fname),0.03);
    end    
    fig_cell(1:n_types,k_i) = im_cell';
end
%%
f=figure;
imshow(cell2mat(fig_cell));
f.PaperUnits = 'centimeters';
f.PaperSize = [18 18];
f.PaperPosition = [0 0 18 18];
set(0, 'currentfigure', f);
print(fullfile(savedir,'brains','all_brains'),'-dpdf');    
