function test_sv_params_mlvarsv
% bvar.sv.sv_params with phi_bnd = .998 reproduces the ml_varsv sample_SVpara
% draw-for-draw, at all three shapes that package calls it with. Settles the
% step-4 record, which had claimed the ml_varsv copy was structurally
% irreproducible: its "no n+r column split, full-vector mu/phi/sig2 indexing"
% differences are no-ops whenever numel(mu) == size(h,2), and every ml_varsv
% call site satisfies that (VAR_CSV.m 61: h is T x 1, mu = 0; VAR_ARSV_redu.m 84
% and VAR_ARSVO_redu.m 91: h is T x n, mu is n x 1; VAR_FSV.m 82: h is
% T x (n+r), mu is (n+r) x 1 - so r = 0 inside sv_params there too). The phi
% truncation bound (.998 here vs the OISV default .999) is the only difference
% that can bite; the last block below shows it biting.
%
% NOT claimed: r > 0. With numel(mu) < size(h,2) the two bodies genuinely
% differ - sv_params demeans only the first n columns, sample_SVpara demeans
% all of them.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg)); %#ok<NASGU>

% chan_koop_yu2024_jbes_oisv ships a same-name, numerically different copy
resolved = which('sample_SVpara');
assert(strncmpi(resolved, leg, numel(leg)), ...
    'sample_SVpara must resolve from the ml_varsv legacy copy, got %s', resolved);

T = 60; nsweep = 12;
shapes = {'VAR-CSV (h Tx1, mu = 0)', 'VAR-SV (h Txn, mu nx1)', 'VAR-FSV (h Tx(n+r), mu (n+r)x1)'};
for ks = 1:3
    switch ks
        case 1      % VAR_CSV.m 61: scalar hyperparameters, mu passed as 0
            m = 1; mu0 = 0;
            Hyper = struct('nuh',3, 'Sh',.1*2, 'phi0',.98, 'Vphi',.05^2, ...
                'mu0',0, 'Vmu',100);
        case 2      % main_varsv.m 89-91 with n = 5
            m = 5; mu0 = log(.4 + (1:m)'/10);
            Hyper = struct('nuh',3*ones(m,1), 'Sh',.1*2*ones(m,1), ...
                'phi0',.98*ones(m,1), 'Vphi',.05^2*ones(m,1), ...
                'mu0',zeros(m,1), 'Vmu',100*ones(m,1));
        case 3      % main_varsv.m 102-104 with n = 5, r = 2
            m = 7; mu0 = log(.4 + (1:m)'/10);
            Hyper = struct('nuh',3*ones(m,1), 'Sh',.1*2*ones(m,1), ...
                'phi0',.98*ones(m,1), 'Vphi',.05^2*ones(m,1), ...
                'mu0',zeros(m,1), 'Vmu',100*ones(m,1));
    end
    phi0 = .9*ones(m,1); if ks == 1, phi0 = .9; end

    % fixed h paths, drawn once so nothing but the sampler consumes rng below
    rng(600+ks, 'twister');
    H = cell(nsweep,1);
    for kk = 1:nsweep
        H{kk} = cumsum(.2*randn(T,m))/5;
    end
    assert(numel(mu0) == size(H{1},2), 'shape %d is not the r = 0 case', ks);

    rng(701, 'twister'); mu = mu0; phi = phi0; got = cell(nsweep,4);
    for kk = 1:nsweep
        [mu,phi,sig2,fl] = bvar.sv.sv_params(H{kk}, mu, phi, Hyper, .998);
        got(kk,:) = {mu,phi,sig2,fl};
    end
    sCore = rng;

    rng(701, 'twister'); mu = mu0; phi = phi0; want = cell(nsweep,4);
    for kk = 1:nsweep
        [mu,phi,sig2,fl] = sample_SVpara(H{kk}, mu, phi, Hyper);
        want(kk,:) = {mu,phi,sig2,fl};
    end
    sLeg = rng;

    assert(isequal(got, want), 'sv_params(.998) differs from ml_varsv sample_SVpara: %s', shapes{ks});
    assert(isequal(sCore.State, sLeg.State), 'rng call sequence differs: %s', shapes{ks});
    if ks > 1
        assert(any(cellfun(@(x) any(x==1), want(:,4))), ...
            'degenerate test, no phi candidate ever accepted: %s', shapes{ks});
    end
end

% the bound is what separates the two copies: force a candidate into
% [.998,.999) and show .999 accepts it where .998 (and the legacy) reject
rng(5, 'twister');
h2 = cumsum(.05*randn(40,1));
HyperT = struct('nuh',3, 'Sh',.1*2, 'phi0',.9985, 'Vphi',1e-8, 'mu0',0, 'Vmu',100);
rng(99,'twister'); [~,pLeg]  = sample_SVpara(h2, 0, .5, HyperT);
rng(99,'twister'); [~,p998]  = bvar.sv.sv_params(h2, 0, .5, HyperT, .998);
rng(99,'twister'); [~,p999]  = bvar.sv.sv_params(h2, 0, .5, HyperT);   % OISV default
assert(isequal(pLeg, p998), 'phi_bnd = .998 must match the legacy bound');
assert(~isequal(pLeg, p999), ...
    'degenerate teeth: no candidate landed in [.998,.999), so the bound was never tested');
end
