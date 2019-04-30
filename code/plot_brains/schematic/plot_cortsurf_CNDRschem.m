function plot_cortsurf_CNDRschem(indices,betas,cmap,alpha,grid,newfig)

% plot cortical atlas from CNDR autopsy protocol

%indices = node indices; 1-6 cortical regions
%betas = values to scale color from
%cmap = colormap like parula or viridis
%alpha = transparency of brains
%grid = [0 0;0 1] to select which brain view(s) to plot,
%[L Lateral, R Lateral; L medial, R medial]
%newfig = T or F, do you want a new figure

%test data
%indices = [1:6,109:112,114:115,234:237]';
%betas = 1:length(indices);

% remove subcortical regions
h_subcort234 = [109:115,227:237]; 
cort_mask = ~ismember(indices,h_subcort234);
indices = indices(cort_mask);
betas = betas(cort_mask);
cmap = cmap(cort_mask,:);

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

if (~exist('grid','var') || isempty(grid))
    grid = [1 1; 0 0];
end

if (~exist('newfig','var') || isempty(newfig))
    newfig = false;
end

%addpath('~/Documents/freesurfmatlab');
addpath(genpath('/Applications/freesurfer/matlab'));
%%
pth = 'data/img/laus_surf/';
[Lvertices,Lfaces]=freesurfer_read_surf([pth,'lh.pial']);
[Lv,LL,Lct] = read_annotation([pth,'lh.myaparc_',num2str(lausannescale),'.annot']);
[Rvertices,Rfaces]=freesurfer_read_surf([pth,'rh.pial']);
[Rv,RL,Rct] = read_annotation([pth,'rh.myaparc_',num2str(lausannescale),'.annot']);

Ls = size(Lct.struct_names,1); Rs = size(Rct.struct_names,1);

for i = 1:Ls
    Lct.struct_names{i} = ['lh_',Lct.struct_names{i}];
end

for i = 1:Rs
    Rct.struct_names{i} = ['rh_',Rct.struct_names{i}];
end

% lind and rind will contain the ROI indices in lausnames, ordered as in
% ?ct.struct_names. i.e. the ROI indices/order you work with outside of
% annotation files, in the order they are found in the color table of the
% annotation file
% this is the key step in converting data ordered the usual lausanne way
% to surface coordinates
% this is not an eloquent way of explaining that fact
[~,lind] = ismember(Lct.struct_names,lausnames);    
[~,rind] = ismember(Rct.struct_names,lausnames);

%convert names to numbers
%make true the indices that match your betas
%replace the trues with the beta values

%%
%indices = randi(234,[25,1]); betas = randi(size(indices,1))*rand(size(indices,1),1);
%betas(round(0.5*size(indices,1)):end) = -1*betas(round(0.5*size(indices,1)):end);
%betas = -1*abs(betas); full = true;
%indices = 1:462; betas = zeros(462,1);

%% 
Lids = false(1,Ls); Rids = false(1,Rs);
Lids(find(ismember(lind,indices))) = true; Rids(find(ismember(rind,indices))) = true; 
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

k = 1; Rinds = []; Rbetas = [];
for i = 1:length(indices)
    if ismember(indices(i),rind)
        Rinds(k,1) = Rct.table(find(ismember(rind,indices(i))),5);
        Rbetas(k,1) = betas(i);
        k = k + 1;
    end
end

%% get faces for all regions not contained in indices in grey

%Lvertices(ismember(LL,Linds)); Rvertices(ismember(RL,Rinds)); %all vertices corresponding to areas to remove
Lremvertind = find(ismember(LL,Linds)); Rremvertind = find(ismember(RL,Rinds)); %indices of all vertices to remove
for i = 1:3
    Lx(:,i) = ismember(Lfaces(:,i),Lremvertind);
    Rx(:,i) = ismember(Rfaces(:,i),Rremvertind);
end
% every entry of Lx/Rx that is a 1 denotes an index in Lfaces/Rfaces
% referencing a vertex for a region that needs to be deleted
 %rows of Lx/Rx with any 1s in them should be deleted because they belong to
 %remvertind (refer to vertices of ROIs in indices/betas)
 
Lfacesnull = Lfaces; Lfacesnull(any(Lx,2),:) = [];
Rfacesnull = Rfaces; Rfacesnull(any(Rx,2),:) = [];

%% get faces for all regions contained in indices on a color scale

Lremvertind = find(~ismember(LL,Linds)); Rremvertind = find(~ismember(RL,Rinds)); %indices of all vertices to remove
%any vertex not corresponding to a member of Linds must be removed if you
%only want to plot regions in indices
for i = 1:3
    Lx(:,i) = ismember(Lfaces(:,i),Lremvertind);
    Rx(:,i) = ismember(Rfaces(:,i),Rremvertind);
end
Lfaces(any(Lx,2),:) = [];
Rfaces(any(Rx,2),:) = [];

colormidpoint = max(betas) - 0.5*(max(betas) - min(betas)); 

%default: center scale around midpoint of betas, ie background value =
%midpoint of betas

[bool,pos] = detectsignchange(betas);

if ((size(betas,1) == 1) || bool) %if only one beta, make background value 0
    colormidpoint = 0; % if betas are all one sign, make background value 0
end

Llabels = repmat(colormidpoint, size(LL,1),1); Rlabels = repmat(colormidpoint, size(RL,1),1);


for i=1:length(Linds)
    Llabels(LL==Linds(i)) = Lbetas(i);
end

for i=1:length(Rinds)
    Rlabels(RL==Rinds(i)) = Rbetas(i);
end

%patch('Faces',Lfaces,'Vertices',Lvertices,'FaceColor','flat','EdgeColor','none','FaceVertexCData',Llabels);
%patch('Faces',Rfaces,'Vertices',Rvertices,'FaceColor','flat','EdgeColor','none','FaceVertexCData',Rlabels);

%% ensure that color scales are the same; not a great fix but?
% will replace 1 label with the lowest beta and 1 label with the highest beta 
% do this on each brain
% do it in the insula where it's tough to see
%
linsula = find(ismember(lausnames,'lh_insula_1')); linsula = Lct.table(find(ismember(lind,linsula)),5);
rinsula = find(ismember(lausnames,'rh_insula_1')); rinsula = Rct.table(find(ismember(rind,rinsula)),5);

linsula = find(LL == linsula); rinsula = find(RL == rinsula);
lminmax = [round(0.25*size(linsula,1)),round(0.755*size(linsula,1))]; % picking two arbitrary points in the L insula
rminmax = [round(0.25*size(rinsula,1)),round(0.755*size(rinsula,1))]; % picking two arbitrary points in the R insula
linsula = linsula(lminmax); rinsula = rinsula(rminmax);
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

Llabels(linsula) = clim; Rlabels(rinsula) = clim; 

%}

clim = [min(betas);max(betas)];
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

%Rbool = (length(Rbetas) > 0);
%Lbool = (length(Lbetas) > 0); % don't plot brains with nothing 

cm = cmap;

if newfig
    figure; 
end

if grid(1,1)
    if sum(sum(grid)) > 1
        subplot(1,2,1);
    end
    patch('Faces',Lfacesnull,'Vertices',Lvertices,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
    patch('Faces',Lfaces,'Vertices',Lvertices,'FaceColor','flat','EdgeColor','none','FaceVertexCData',Llabels); colormap(cm);
    view([-1 0 0]); daspect([1,1,1]); 
    material dull
    camlight; caxis(clim); axis off;
    nogridtick(); set(gca,'fontsize',0.00000001); %set(gca,'Visible','off')
end

if grid(1,2)
    if sum(sum(grid)) > 1
        subplot(1,2,2); 
    end
    patch('Faces',Lfacesnull,'Vertices',Lvertices,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
    patch('Faces',Lfaces,'Vertices',Lvertices,'FaceColor','flat','EdgeColor','none','FaceVertexCData',Llabels); colormap(cm);
    view([1 0 0]); daspect([1,1,1]);
    material dull
    camlight; caxis(clim); axis off;
    nogridtick(); set(gca,'fontsize',0.00000001); %set(gca,'Visible','off')
end

%}