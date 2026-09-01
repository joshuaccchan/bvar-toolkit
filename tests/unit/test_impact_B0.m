function test_impact_B0
% deterministic function: exact equivalence with the legacy ml_varsv copy
% (chan2023_joe_mlvarsv prior_B0)
rng(13, 'twister');
n = 4; T = 60;
Y0 = randn(9, n);    % presample; only the last 4 rows enter
Y = randn(T, n);
kappa = 0.2;

[b1, V1] = bvt.priors.impact_B0(Y0, Y, kappa);

root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
[b2, V2] = prior_B0(Y0, Y, kappa);
assert(isequal(b2, b1), 'impact_B0: beta0 differs from legacy prior_B0');
assert(isequal(V2, V1), 'impact_B0: Vbeta differs from legacy prior_B0');
end
