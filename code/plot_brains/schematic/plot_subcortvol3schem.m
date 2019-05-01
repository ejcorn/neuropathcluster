function f = plot_subcortvol3schem(indices,cmap)

%indices: Lausanne 234 labels of regions corresponding to betas
%cmap: Nx3 matrix of rgb values, same length as indices, specifies region
%colors

%note that fliplr(subcortnifti) was used to ensure that the regions would
%be plotted on the right side, though it is unclear to me why they were not
%obeying the nifti coordinates in the first place

if (~exist('cmap','var'))
    disp('must provide color values')
    return
end

% test data
%cmap = roi_colors; indices = indices;

[~,ord] = sort(indices);  % reorder indices/betas in ascending order of indices
indices = indices(ord);
cmap = cmap(ord,:);
indices = reshape(indices,[],1); % convert inputs to row vectors

laus234nifti = load_nii('~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/data/img/ROIv_scale125_dilatedBS.nii.gz');
%% duplicate left and right

h_subcort234 = [109:112,114:115,227:230,232:237];
sc = ismember(indices,h_subcort234);
indices = indices(sc); cmap = cmap(sc,:); %eliminate non-subcortical ROIs

% Duplicate subcortical structures bilaterally
% default input has only right-sided SC indices
right_side_sc = 109:115;
left_side_sc = 227:233; 
pres_mask = ismember(right_side_sc,indices); % find indices of right sided regions that are present
indices = [indices;left_side_sc(pres_mask)']; % add left side to indices
right_cmap = cmap(ismember(indices,right_side_sc),:); % get right side betas
cmap = [cmap;right_cmap]; % add corresponding elements to left side to rgbbetas

%%

views = {[-1 0 0],[0 1 0]};
n_views = length(views);
s = size(indices,1);
z = size(laus234nifti.img); 

f = figure; view(3);colormap(cmap);
for j = 1:n_views
    subplot(1,n_views,j);
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
        p = patch('Faces',surf.faces,'Vertices',surf.vertices,'FaceColor',cmap(i,:),'EdgeColor','none');
    end %one by one, plot subcortical ROIs

    daspect([1 1 1])
    view(views{j})
    
    %middle head on = 1 0 0, left lat = 0 1 0, right lat= 0 -1 0
    material dull
    camlight
    lighting gouraud
    set(gca,'FontSize',8);
end
