function test_prior_niw_opthyper_ncp
% deterministic: bvar.priors.niw variant 'opthyper_ncp' reproduces
% AD_OptHyper's 5-kappa prior_NCP exactly - non-integer lag-decay exponent
% kappa(2), IW hyperparameters kappa(4)/kappa(5), sparse-diagonal VA0
rng(45, 'twister');
n = 3; p = 4; Tt = 64;
Y0 = randn(6, n); Yt = randn(Tt, n);
kappa = [0.04, 1.7, 100, 3, 0.8];
[A1, V1, nu1, S1] = bvar.priors.niw(p, kappa, Y0, Yt, 'opthyper_ncp');

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'cjz2019_ad_opthyper', 'legacy');
addpath(leg); c = onCleanup(@() rmpath(leg));
[A2, V2, nu2, S2] = prior_NCP(p, kappa, Y0, Yt);
assert(isequal(A2, A1), 'niw: A0 differs from AD_OptHyper prior_NCP');
assert(isequal(V2, V1), 'niw: VA0 differs from AD_OptHyper prior_NCP');
assert(issparse(V1) && issparse(V2), 'niw: VA0 must be sparse as in the legacy copy');
assert(isequal(nu2, nu1), 'niw: nu0 differs from AD_OptHyper prior_NCP');
assert(isequal(S2, S1), 'niw: S0 differs from AD_OptHyper prior_NCP');
end
