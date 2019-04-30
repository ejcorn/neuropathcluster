function zr_scores = GET_PAIRWISE_UNIQUE_ZRAND(partitions)

% partitions should be pxn
% p = number of nodes
% n = number of partitions
% requires zrand.m from NCT
nreps = size(partitions,2);
nCombos = nreps*(nreps-1)/2;
combos = nchoosek(1:nreps,2);
combos1 = combos(:,1); combos2 = combos(:,2);
zr_scores = zeros(nCombos,1);
for i = 1:nCombos
    part1 = combos1(i); part2 = combos2(i);
    zr_scores(i) = zrand(double(partitions(:,part1)),double(partitions(:,part2)));
end
