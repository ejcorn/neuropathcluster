clear all; close all; clc

BCT_path = '~/Dropbox/Cornblath_Bassett_Projects/code/BCT';
homedir = '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/code/plot_brains';
addpath(genpath(pwd));
%% Plot schematic indicating which brain regions were sampled

indices = [1:6,109:112,114:115,234:237]';
load(fullfile(homedir,'schematic/regioncolors.mat'))
f=plot_subcortvol3schem(indices,roi_colors);
f.PaperUnits = 'inches';
f.PaperSize = [1.75 1];
f.PaperPosition = [0 0 1.75 1];
set(0,'CurrentFigure',f); set(gca,'FontSize',8);
fname = fullfile(homedir,['schematic/subcortex.png']);
print(fname,'-dpng','-r400'); close(f);
 
f=figure;
plot_cortsurf_CNDRschem(indices,1:length(indices),roi_colors)
f.PaperUnits = 'inches';
f.PaperSize = [1.75 1];
f.PaperPosition = [0 0 1.75 1];
set(0,'CurrentFigure',f);
fname = fullfile(homedir,['schematic/cortex.png']);
print(fname,'-dpng','-r400');