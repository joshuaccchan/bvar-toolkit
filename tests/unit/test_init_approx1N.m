function test_init_approx1N
% deterministic function: exact equivalence with the legacy ml_varsv copy
rng(6, 'twister');
s2 = exp(randn(80, 1));
h1 = bvar.sv.init_approx1N(s2, -1.0, 0.95, 0.2);

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
h2 = getARh_approx1N(s2, -1.0, 0.95, 0.2);
assert(isequal(h2, h1), 'init_approx1N: differs from legacy getARh_approx1N');
end
