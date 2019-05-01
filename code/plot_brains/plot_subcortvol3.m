function [f_ant,f_lat] = plot_subcortvol3(indices,betas,cmap)

%indices: Lausanne 234 labels of regions corresponding to betas
%betas: values to base color scale off of, corresponding to indices
%cmap: name of colormap to use

%f_ant = figure handle for anterior view of subcortex
%f_med = figure handle for lateral view of subcortex

%note that fliplr(subcortnifti) was used to ensure that the regions would
%be plotted on the right side, though it is unclear to me why they were not
%obeying the nifti coordinates in the first place


if (~exist('cmap','var'))
    cmap = 'parula';
end

%test data
%indices = C_region_plot_t(:,2); betas = C_region_plot_t(:,1);
%indices = [227:236]'; betas = randi(10)*rand(size(indices,1),1); cmap = 'plasma';
%betas(round(0.5*size(indices,1)):end) = -1*betas(round(0.5*size(indices,1)):end);
%betas = 1*abs(betas);

[~,ord] = sort(indices);  % reorder indices/betas in ascending order of indices
indices = indices(ord);
betas = betas(ord);
betas = reshape(betas,[],1); indices = reshape(indices,[],1); % convert inputs to row vectors

laus234nifti = load_nii('~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/data/img/ROIv_scale125_dilatedBS.nii.gz');
%% convert Nx1 betas to a Nx3 matrix of RGB values

betas(end+1) = 5;
betas = [1; betas];

C = colormap(cmap); close; L = size(C,1);
cmap_beta_idx = round(interp1(linspace(min(betas(:)),max(betas(:)),L),1:L,betas));
% be aware that this scaling method works worse on log/less linear data
%
cmap_beta_idx([1 length(betas)],:) = [];
betas([1 length(betas)],:) = [];

h_subcort234 = [109:112,114:115,227:230,232:237]; % Lausanne subcortical ROIs excluding accumbens area
sc = ismember(indices,h_subcort234);
indices = indices(sc); betas = betas(sc,:); %eliminate non-subcortical ROIs

% Duplicate subcortical structures bilaterally
% default input has only right-sided SC indices
right_side_sc = 109:115;
left_side_sc = 227:233; 
pres_mask = ismember(right_side_sc,indices); % find indices of right sided regions that are present
indices = [indices;left_side_sc(pres_mask)']; % add left side to indices
right_betas = betas(ismember(indices,right_side_sc)); % get right side betas
betas = [betas;right_betas]; % add corresponding elements to left side to rgbbetas

%%
cmin = 1; cmax = 5;
views = {[-1 0 0],[0 1 0]};
n_views = length(views);
s = size(indices,1);
z = size(laus234nifti.img); 


for j = 1:n_views
    if j == 1
        f_ant = figure;
        set(0,'CurrentFigure',f_ant);
    elseif j == 2
        f_lat = figure;
        set(0,'CurrentFigure',f_lat);
    end
    set(gca,'Visible','off');
    for i = 1:length(h_subcort234)
        sc_grey_roi_idx = h_subcort234(i);
        subcortnifti = zeros(z);
        ind = ismember(laus234nifti.img,sc_grey_roi_idx);
        subcortnifti(ind) = i;
        surf = isosurface(fliplr(subcortnifti),0); %fliplr accounts for LR swap from nifti file coordinates
        p = patch('Faces',surf.faces,'Vertices',surf.vertices,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none');
    end % gray background

    for i = 1:length(indices)
        sc_col_roi_idx = indices(i);
        subcortnifti = zeros(z);
        ind = ismember(laus234nifti.img,sc_col_roi_idx);
        subcortnifti(ind) = i;
        surf = isosurface(fliplr(subcortnifti),0);
        p = patch('Faces',surf.faces,'Vertices',surf.vertices,'EdgeColor','none','FaceVertexCData',betas(i)*ones(length(surf.vertices),1),'FaceColor','flat');
    end %one by one, plot subcortical ROIs

    daspect([1 1 1])
    view(views{j})
    %colorbar('southoutside')
    colormap(cmap);
    caxis([cmin cmax])
    
    %middle head on = 1 0 0, left lat = 0 1 0, right lat= 0 -1 0
    material dull
    camlight
    lighting gouraud
end
