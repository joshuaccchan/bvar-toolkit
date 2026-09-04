function test_prior_minn_largebvar
% deterministic: bvar.priors.minn with n0pre = p reproduces large_BVAR's
% prior_Minn exactly (p must be 4 - the AR(4) design is only conformable then,
% and every legacy caller uses p = 4)
rng(41, 'twister');
n = 3; p = 4; Tt = 60;
Y0 = randn(6, n); Yt = randn(Tt, n);
c1 = 0.2^2; c2 = 0.05^2; c3 = 100;
[b1, V1, S1] = bvar.priors.minn(p, c1, c2, c3, Y0, Yt, p);

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2020_springer_largebvar', 'legacy');
addpath(leg); c = onCleanup(@() rmpath(leg));
[b2, V2, S2] = prior_Minn(p, c1, c2, c3, Y0, Yt);
assert(isequal(b2, b1), 'minn: beta differs from large_BVAR prior_Minn');
assert(isequal(V2, V1), 'minn: V differs from large_BVAR prior_Minn');
assert(isequal(S2, S1), 'minn: Sig_hat differs from large_BVAR prior_Minn');
end
