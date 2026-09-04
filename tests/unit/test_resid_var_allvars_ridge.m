function test_resid_var_allvars_ridge
% deterministic function: exact equivalence with the legacy HYB get_resid_var_v2
rng(42, 'twister');
n = 5;
Y0 = randn(9, n);
Y  = randn(60, n);
s_core = bvar.priors.resid_var_allvars_ridge(Y0, Y);

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_jbes_hybtvp', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
s_leg = get_resid_var_v2(Y0, Y);
assert(isequal(s_leg, s_core), ...
    'resid_var_allvars_ridge: differs from legacy get_resid_var_v2');

% never-merge guard: this is NOT the univariate AR(4) variant
assert(~isequal(s_core, bvar.priors.resid_var_ar4(Y0, Y)), ...
    'resid_var_allvars_ridge: unexpectedly equal to resid_var_ar4 - variants merged?');
end
