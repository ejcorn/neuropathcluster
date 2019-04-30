function f = figure_resize(f,units,sz)

% f: figure handle
% units: character variable with either 'inches' or 'centimeters' or
% another valid unit
% sz; double, [width height]

f.PaperUnits = units;
f.PaperSize = sz;
f.PaperPosition = [0 0 sz];