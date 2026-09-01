function test_acp_2019wp_relationship
% documentation test: the 2019wp standalone prior_ACPi (chan2019wp_acp) is a
% strict earlier subset of the R1 local prior_ACPi inside bvt.priors.acp_stru -
% same mi/Vi/nui/Si values as is_ns = false, but Vi packaged as a sparse
% diagonal matrix. It is NOT canonicalized; this test pins the relationship.
rng(14, 'twister');
n = 4; p = 2;
kappa = [0.04, 0.0016, 1, 100];
sig2 = exp(randn(n, 1));
prior = bvt.priors.acp_stru(n, p, kappa, sig2);   % idx_ns = [] -> is_ns false

root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan2019wp_acp', 'legacy');
addpath(leg); c = onCleanup(@() rmpath(leg));
kb = n*p + 1;
cnt = 0;
for ii = 1:n
    [mi, Vi, nui, Si] = prior_ACPi(n, p, ii, kappa, sig2);
    assert(issparse(Vi) && isequal(size(Vi), [kb+ii-1, kb+ii-1]), ...
        'acp 2019wp: expected sparse diagonal Vi of size ki x ki');
    v = full(diag(Vi));
    assert(isequal(mi(1:kb), prior.beta0(:, ii)) && ...
        isequal(mi(kb+1:end), prior.alp0(cnt+1:cnt+ii-1)), ...
        'acp 2019wp: prior means differ from R1 values');
    assert(isequal(v(1:kb), prior.Vbeta(:, ii)) && ...
        isequal(v(kb+1:end), prior.Valp(cnt+1:cnt+ii-1)), ...
        'acp 2019wp: prior variances differ from R1 values');
    assert(isequal(nui, prior.nu(ii)) && isequal(Si, prior.S(ii)), ...
        'acp 2019wp: nu/S differ from R1 values');
    cnt = cnt + ii - 1;
end
end
