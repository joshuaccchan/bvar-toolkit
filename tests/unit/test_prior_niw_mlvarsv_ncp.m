function test_prior_niw_mlvarsv_ncp
% deterministic: bvar.priors.niw variant 'mlvarsv_ncp' reproduces ml_varsv's
% prior_NCP exactly, including U_hat; p = 2 pins the hard-coded 4-row
% presample branch
rng(44, 'twister');
n = 4; p = 2; Tt = 58;
Y0 = randn(7, n); Yt = randn(Tt, n);
c1 = 0.04; c2 = 100;
[A1, V1, nu1, S1, U1] = bvar.priors.niw(p, [c1 c2], Y0, Yt, 'mlvarsv_ncp');

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
[A2, V2, nu2, S2, U2] = prior_NCP(p, c1, c2, Y0, Yt);
assert(isequal(A2, A1), 'niw: A0 differs from ml_varsv prior_NCP');
assert(isequal(V2, V1), 'niw: VA0 differs from ml_varsv prior_NCP');
assert(isequal(nu2, nu1), 'niw: nu0 differs from ml_varsv prior_NCP');
assert(isequal(S2, S1), 'niw: S0 differs from ml_varsv prior_NCP');
assert(isequal(U2, U1), 'niw: U_hat differs from ml_varsv prior_NCP');
end
