function [f_lat,f_med] = plot_cortsurf_CNDRautopsy(indices,betas,cmap,alpha)

% plot cortical atlas from CNDR autopsy protocol

%indices = node indices; 1-6 cortical regions
%betas = values to scale color from
%cmap = colormap like parula or viridis
%alpha = transparency of brains

%f_lat = figure handle for lateral view of L hemisphere
%f_med = figure handle for medial view of L hemisphere

% remove subcortical regions
h_subcort234 = [109:115,227:237]; 
cort_mask = ~ismember(indices,h_subcort234);
indices = indices(cort_mask);
betas = betas(cort_mask);

lausannescale = 125; % no reason to use a larger scale
load data/img/human_regionNames.mat
lausnames = roinames{3};

% expand indices to CNDR autopsy parcellation
[indices,betas] = expand_indices_CNDR(indices,betas,lausnames);

if (~exist('full','var'))
    full = true;
end
% full sets whether or not to plot brains with nothing on them

if (~exist('cmap','var'))
    cmap = 'parula';
end

if (~exist('alpha','var') || isempty(alpha))
    alpha = 1;
end

%addpath('~/Documents/freesurfmatlab');
addpath(genpath('/Applications/freesurfer/matlab'));
%%
pth = 'data/img/laus_surf/';
[Lvertices,Lfaces]=freesurfer_read_surf([pth,'lh.pial']);
[Lv,LL,Lct] = read_annotation([pth,'lh.myaparc_',num2str(lausannescale),'.annot']);
Ls = size(Lct.struct_names,1);

for i = 1:Ls
    Lct.struct_names{i} = ['lh_',Lct.struct_names{i}];
end

% lind and rind will contain the ROI indices in lausnames, ordered as in
% ?ct.struct_names. i.e. the ROI indices/order you work with outside of
% annotation files, in the order they are found in the color table of the
% annotation file
% this is the key step in converting data ordered the usual lausanne way
% to surface coordinates
% this is not an eloquent way of explaining that fact
[~,lind] = ismember(Lct.struct_names,lausnames);    

%convert names to numbers
%make true the indices that match your betas
%replace the trues with the beta values

%%
%indices = randi(234,[25,1]); betas = randi(size(indices,1))*rand(size(indices,1),1);
%betas(round(0.5*size(indices,1)):end) = -1*betas(round(0.5*size(indices,1)):end);
%betas = -1*abs(betas); full = true;
%indices = 1:462; betas = zeros(462,1);

%% 
Lids = false(1,Ls);
Lids(find(ismember(lind,indices))) = true;
%Lids is a logical the size of Lct.table (contains node labels, reference
%for LL, which labels vertices). Lids tells you which nodes in indices are
%on the left brain.
k = 1; Linds = []; Lbetas = [];
for i = 1:length(indices)
    if ismember(indices(i),lind)
        Linds(k,1) = Lct.table(find(ismember(lind,indices(i))),5);
        Lbetas(k,1) = betas(i);
        k = k + 1;
    end
end

%% get faces for all regions not contained in indices in grey

%Lvertices(ismember(LL,Linds)); Rvertices(ismember(RL,Rinds)); %all vertices corresponding to areas to remove
Lremvertind = find(ismember(LL,Linds)); %indices of all vertices to remove
for i = 1:3
    Lx(:,i) = ismember(Lfaces(:,i),Lremvertind);    
end
% every entry of Lx/Rx that is a 1 denotes an index in Lfaces/Rfaces
% referencing a vertex for a region that needs to be deleted
 %rows of Lx/Rx with any 1s in them should be deleted because they belong to
 %remvertind (refer to vertices of ROIs in indices/betas)
 
Lfacesnull = Lfaces; Lfacesnull(any(Lx,2),:) = [];

%% get faces for all regions contained in indices on a color scale

Lremvertind = find(~ismember(LL,Linds)); %indices of all vertices to remove
%any vertex not corresponding to a member of Linds must be removed if you
%only want to plot regions in indices
for i = 1:3
    Lx(:,i) = ismember(Lfaces(:,i),Lremvertind);    
end
Lfaces(any(Lx,2),:) = [];

colormidpoint = max(betas) - 0.5*(max(betas) - min(betas)); 

%default: center scale around midpoint of betas, ie background value =
%midpoint of betas

[bool,pos] = detectsignchange(betas);

if ((size(betas,1) == 1) || bool) %if only one beta, make background value 0
    colormidpoint = 0; % if betas are all one sign, make background value 0
end

Llabels = repmat(colormidpoint, size(LL,1),1);


for i=1:length(Linds)
    Llabels(LL==Linds(i)) = Lbetas(i);
end

%% ensure that color scales are the same; not a great fix but?
% will replace 1 label with the lowest beta and 1 label with the highest beta 
% do this on each brain
% do it in the insula where it's tough to see
%
linsula = find(ismember(lausnames,'lh_insula_1')); linsula = Lct.table(find(ismember(lind,linsula)),5);

linsula = find(LL == linsula);
lminmax = [round(0.25*size(linsula,1)),round(0.755*size(linsula,1))]; % picking two arbitrary points in the L insula
linsula = linsula(lminmax);
clim = [min(betas);max(betas)];

if (length(betas) == 1)
    if betas > 0
        clim = [0; betas];
    end
    if betas < 0
        clim = [betas; 0];
    end
end

if bool
    if ~pos
        clim = [min(betas); -1*min(betas)];
    end
    if pos
        clim = [min(betas); max(betas)];
    end
end

Llabels(linsula) = clim;

%}

%clim = [min(betas);max(betas)];
%% round clim to n sig figs

%{
nsig = 2;
clim = [round(min(betas),nsig,'significant');round(max(betas),nsig,'significant')];

if (length(betas) == 1)
    if betas > 0
        clim = [0; round(betas,nsig,'significant')];
    end
    if betas < 0
        clim = [round(betas,nsig,'significant'); 0];
    end
end

if bool
    if ~pos
        clim = [round(min(betas),nsig,'significant'); -1*round(min(betas),nsig,'significant')];
    end
    if pos
        clim = [-1*round(max(betas),nsig,'significant'); round(max(betas),nsig,'significant')];
    end
end


Llabels(linsula) = clim; Rlabels(rinsula) = clim; 
%}

%% subplot version for viewing

cm = cmap;
clim = [1 5]; % manually set scale to 1-5, based on possible scores

f_lat = figure;
patch('Faces',Lfacesnull,'Vertices',Lvertices,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
patch('Faces',Lfaces,'Vertices',Lvertices,'FaceColor','flat','EdgeColor','none','FaceVertexCData',Llabels); colormap(cm);
view([-1 0 0]); daspect([1,1,1]); 
material dull
camlight; caxis(clim); axis off;
nogridtick(); set(gca,'fontsize',0.00000001); %set(gca,'Visible','off')

f_med = figure;
patch('Faces',Lfacesnull,'Vertices',Lvertices,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
patch('Faces',Lfaces,'Vertices',Lvertices,'FaceColor','flat','EdgeColor','none','FaceVertexCData',Llabels); colormap(cm);
view([1 0 0]); daspect([1,1,1]);
material dull
camlight; caxis(clim); axis off;
nogridtick(); set(gca,'fontsize',0.00000001); %set(gca,'Visible','off')

%}