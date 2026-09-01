function test_prior_minn_mlvarsv
% deterministic: bvt.priors.minn with n0pre = 4 reproduces ml_varsv's
% prior_Minn exactly, including the U_hat output; p = 2 pins the hard-coded
% 4-row presample branch (differs from n0pre = p there)
rng(42, 'twister');
n = 4; p = 2; Tt = 55;
Y0 = randn(7, n); Yt = randn(Tt, n);
c1 = 0.04; c2 = 0.01; c3 = 100;
[b1, V1, S1, U1] = bvt.priors.minn(p, c1, c2, c3, Y0, Yt, 4);

root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
[b2, V2, S2, U2] = prior_Minn(p, c1, c2, c3, Y0, Yt);
assert(isequal(b2, b1), 'minn: beta differs from ml_varsv prior_Minn');
assert(isequal(V2, V1), 'minn: V differs from ml_varsv prior_Minn');
assert(isequal(S2, S1), 'minn: Sig_hat differs from ml_varsv prior_Minn');
assert(isequal(U2, U1), 'minn: U_hat differs from ml_varsv prior_Minn');

% and at p = 4 the two documented settings coincide
p4 = 4;
[b3, V3, S3] = bvt.priors.minn(p4, c1, c2, c3, Y0, Yt, p4);
[b4, V4, S4] = prior_Minn(p4, c1, c2, c3, Y0, Yt);
assert(isequal(b4, b3) && isequal(V4, V3) && isequal(S4, S3), ...
    'minn: n0pre = p and n0pre = 4 disagree at p = 4');
end
