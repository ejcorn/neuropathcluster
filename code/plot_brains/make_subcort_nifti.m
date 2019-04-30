clear all; close all; clc
%% add relevant paths

addpath(genpath('/Applications/freesurfer/matlab'));
addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/brainmapping2'));
addpath(genpath('code'));
% Regions not included in Lausanne
% SN, LC, Medulla, Pons, Midbrain, Cerebellum
lausniftipath = 'data/img/ROIv_scale125_dilated.nii.gz';

%% make z-plane cutoffs in brainstem for midbrain, pons, medulla
% make midbrain 234, pons 235, medulla 236
nii = load_nii(lausniftipath);
nii_new = nii;

% view_nii(nii2); % visualize brainstem, identify upper and lower limit of
% pons crossing fibers
upper_pons_z = 27; 
lower_pons_z = 12;

% demarcate pons: below midbrain, above medulla
pons_z_mask = false(size(nii_new.img)); 
pons_z_mask(:,:,lower_pons_z:upper_pons_z) = true;
nii_new.img(nii_new.img == 234 & pons_z_mask) = 235;

% demarcate medulla: below pons
medulla_z_mask = false(size(nii_new.img)); 
medulla_z_mask(:,:,1:lower_pons_z) = true;
nii_new.img(nii_new.img == 234 & medulla_z_mask) = 236;

%% add cerebellum
cb = load_nii('data/img/Cerebellum-MNIflirt-maxprob-thr50-1mm.nii.gz');
cb.img(cb.img > 0) = 237; % binarize cerebellum
cb.img = int32(cb.img);
nii_new.img = nii_new.img + cb.img(1:2:end,1:2:end,1:2:end);
unique(nii_new.img)

%view_nii(nii_new)
%%
new_fname = 'data/img/ROIv_scale125_dilatedBS.nii.gz';
save_nii(nii_new,new_fname);

%% plot brains
cmap = 'plasma'; rank = true;
indices = [109:115,200:237]; 
indices = 116:237;
values = randperm(length(indices));
f = plot_subcortvol3(indices,values,cmap,rank)