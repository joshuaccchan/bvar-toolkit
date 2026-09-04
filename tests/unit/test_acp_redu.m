function test_acp_redu
% deterministic function: exact equivalence with the legacy R1 copy
% (chan2022_qe_acp prior_ACP_redu), both idx_ns branches and the
% asymmetric-kappa grid pattern used in main_ACP_jointden
rng(12, 'twister');
n = 5; p = 3;
kappa = [0.04, 0.0016, 1, 100];
kappa_asym = [0.1, 0.0009, 1, 100];
sig2 = exp(randn(n, 1));
idx_ns = [1 3 5];

p_ns = bvar.priors.acp_redu(n, p, kappa, sig2, idx_ns);
p_st = bvar.priors.acp_redu(n, p, kappa, sig2);   % nargin==4 branch, idx_ns = []
p_as = bvar.priors.acp_redu(n, p, kappa_asym, sig2, idx_ns);

root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2022_qe_acp', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
q_ns = prior_ACP_redu(n, p, kappa, sig2, idx_ns);
q_st = prior_ACP_redu(n, p, kappa, sig2);
q_as = prior_ACP_redu(n, p, kappa_asym, sig2, idx_ns);
assert(isequal(q_ns, p_ns), 'acp_redu: differs from legacy prior_ACP_redu (idx_ns set)');
assert(isequal(q_st, p_st), 'acp_redu: differs from legacy prior_ACP_redu (idx_ns default)');
assert(isequal(q_as, p_as), 'acp_redu: differs from legacy prior_ACP_redu (asymmetric kappa)');
end
