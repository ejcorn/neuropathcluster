function x = nanzscore(x,d)
if ~exist('d','var')
    d=1; % z-score across rows by default
end
x = (x-nanmean(x,d))./nanstd(x,[],d);