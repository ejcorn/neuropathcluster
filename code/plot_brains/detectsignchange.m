function [bool, pos] = detectsignchange(betas)
% used to determine if range of betas for plotting on brainmaps cross zero
% if all betas are on one side of zero, then color scale doesn't reflect
% directionality
% thus, need to detect sets of betas with only positive or only negative
% values and return bool if true
% if that is the case, i will apply a colorscale that reflects the values
% over 0
% bool is true is there is no sign change, pos is true if all elements are
% positive
% bool will be false if betas have length of 1

s = length(betas);
bool = false; pos = false;
if (sum(betas > 0) == s)
    bool = true;
    pos = true;
end

if (sum(betas < 0) == s)
    bool = true;
end

if s == 1
    bool = false; 
end