function test_prior_niw_largebvar_nc
% deterministic: bvar.priors.niw variant 'largebvar_nc' reproduces large_BVAR's
% prior_NC exactly (p = 4 as in all legacy callers)
rng(43, 'twister');
n = 3; p = 4; Tt = 62;
Y0 = randn(6, n); Yt = randn(Tt, n);
c1 = 0.2^2; c2 = 100;
[A1, V1, nu1, S1] = bvar.priors.niw(p, [c1 c2], Y0, Yt, 'largebvar_nc');

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2020_springer_largebvar', 'legacy');
addpath(leg); c = onCleanup(@() rmpath(leg));
[A2, V2, nu2, S2] = prior_NC(p, c1, c2, Y0, Yt);
assert(isequal(A2, A1), 'niw: A0 differs from large_BVAR prior_NC');
assert(isequal(V2, V1), 'niw: VA0 differs from large_BVAR prior_NC');
assert(isequal(nu2, nu1), 'niw: nu0 differs from large_BVAR prior_NC');
assert(isequal(S2, S1), 'niw: S0 differs from large_BVAR prior_NC');
end
