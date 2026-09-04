function test_surform
% structure + exact equivalence with the legacy HYB copy
rng(1, 'twister');
X = randn(7, 3);
Z = bvar.util.surform(X);
assert(isequal(size(Z), [7, 21]), 'surform: wrong size');
assert(isequal(full(Z(2, 4:6)), X(2, :)), 'surform: wrong block placement');
assert(nnz(Z) == numel(X), 'surform: wrong sparsity');

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_jbes_hybtvp', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
assert(isequal(SURform(X), Z), 'surform: differs from legacy SURform');
end
