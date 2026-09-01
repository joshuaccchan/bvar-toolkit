function test_acp_stru
% deterministic function: exact equivalence with the legacy R1 copy
% (chan2022_qe_acp prior_ACP_stru), both idx_ns branches
rng(11, 'twister');
n = 5; p = 3;
kappa = [0.04, 0.0016, 1, 100];
sig2 = exp(randn(n, 1));
idx_ns = [2 4];

p_ns = bvt.priors.acp_stru(n, p, kappa, sig2, idx_ns);
p_st = bvt.priors.acp_stru(n, p, kappa, sig2);   % nargin==4 branch, idx_ns = []

root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan2022_qe_acp', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
q_ns = prior_ACP_stru(n, p, kappa, sig2, idx_ns);
q_st = prior_ACP_stru(n, p, kappa, sig2);
assert(isequal(q_ns, p_ns), 'acp_stru: differs from legacy prior_ACP_stru (idx_ns set)');
assert(isequal(q_st, p_st), 'acp_stru: differs from legacy prior_ACP_stru (idx_ns default)');
end
