function test_prior_niw_kron_script
% deterministic: bvt.priors.niw variant 'kron_script' with kappa = [.2^2 100]
% reproduces the BVAR_code workspace script construct_prior_A exactly
% (the script hard-codes c1 = .2^2, c2 = 100 and computes only A0/VA0/sig2;
% its callers all set S0 = eye(n), nu0 = n+3 just before running it)
rng(46, 'twister');
n = 3; p = 4; T = 61;
Y0 = randn(4, n); shortY = randn(T, n);
k = n*p + 1;
[A1, V1, nu1, S1] = bvt.priors.niw(p, [.2^2 100], Y0, shortY, 'kron_script');

root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan2020_jbes_kronecker', 'legacy');
addpath(leg); c = onCleanup(@() rmpath(leg));
% pre-declare the variables the script assigns so the parser binds them as
% variables in this workspace, then run the script exactly as the legacy
% mains do (with k, n, T, p, Y0, shortY already set)
A0 = []; VA0 = []; sig2 = []; %#ok<NASGU>
construct_prior_A;
assert(isequal(A0, A1), 'niw: A0 differs from construct_prior_A');
assert(isequal(VA0, V1), 'niw: VA0 differs from construct_prior_A');
assert(isequal(nu1, n+3), 'niw: nu0 must be n+3 as set by the legacy callers');
assert(isequal(S1, eye(n)), 'niw: S0 must be eye(n) as set by the legacy callers');
end
