function [exp_inds,exp_betas] = expand_indices_CNDR(indices,betas,lausnames)

% go from indices 1:6 of CNDR autopsied cortical regions
% expand to corresponding Lausanne 125 indices
% expand (i.e. copy) corresponding beta values

exp_inds = []; exp_betas = [];
for i = 1:length(indices)
    idx = indices(i);
    if idx == 1
        query = {'anteriorcingulate'};        
    elseif idx == 2
        query = {'occip','pericalc'};
    elseif idx == 3
        query = {'middletemporal','superiortemporal'};
    elseif idx == 4
        query = {'middlefrontal'};
    elseif idx == 5
        %query = {'inferiorparietal'};
        query = cellfun(@(x) sprintf('inferiorparietal_%d',x),num2cell([2:6]),'UniformOutput',false);
    elseif idx == 6
        query = {'entorhinal'};        
    end
    new_inds = find(contains(lausnames,query));
    % concatenate all parcels
    exp_inds = [exp_inds;new_inds];
    % repeat betas for each added parcel    
    exp_betas = [exp_betas;ones(length(new_inds),1)*betas(i)];
end
