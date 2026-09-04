function test_tnormrnd
% bounds + exact draw-for-draw equivalence with the legacy copy under one seed
rng(4, 'twister');
t = bvar.util.tnormrnd(0.3, 2.0, -1, 1.5, 500);
assert(all(t >= -1 & t <= 1.5), 'tnormrnd: draws outside bounds');

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
rng(4, 'twister');
t1 = tnormrnd(0.3, 2.0, -1, 1.5, 500);       % legacy (path-shadowed name)
assert(isequal(t1, t), 'tnormrnd: differs from legacy under the same seed');
end
