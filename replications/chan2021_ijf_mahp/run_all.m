% chan2021_ijf_mahp/run_all - functionized estimation pipeline of Chan (2021,
% IJF): large BVAR with SV under the Minnesota-type normal-gamma ('MNG'),
% normal-gamma ('NG'), or data-based Minnesota ('Minn') prior.
%
%   out = run_all(model, nsim, burnin, seed)
%
%   model  - 'MNG' | 'NG' | 'Minn' (default 'MNG', the legacy main_BVAR.m
%            line 17 default model = 1)
%   nsim   - posterior draws kept (default preset nsim_default = 10000)
%   burnin - burn-in sweeps (default preset burnin_default = 1000)
%   seed   - optional; when nonempty, rng(seed,'twister') is set FIRST, before
%            any draw. When omitted/empty the ambient rng state is used as-is.
%
% Functionized 2026-09-01 (step 5, MAHP flagship). Reproduces the legacy
% pipeline main_BVAR.m -> BVAR_MNG.m / BVAR_NG.m / BVAR_Minn.m draw-for-draw
% bitwise (verified by tests/unit/test_mahp_equivalence.m at small nsim) with
% ONE deliberate divergence: the legacy scripts re-seed the global stream from
% the wall clock (randn('seed',sum(clock*100)); rand('seed',sum(clock*1000)) -
% BVAR_MNG.m line 32, BVAR_NG.m line 31, BVAR_Minn.m line 25), which makes
% as-shipped runs irreproducible AND switches MATLAB to the legacy v4/v5
% generators; run_all drops that line so the caller controls seeding via
% `seed` (or the ambient state) on the modern twister stream. The legacy
% clock-seed sits AFTER the psi gamrnd chain-init draws, so under the dropped
% line the whole run (init draws included) consumes one coherent stream.
% The trailing legacy histogram figures are not reproduced.
%
% All constants come from preset.m in this folder (each field cites its legacy
% source line); the data file is read from legacy/ READ-ONLY; the Gibbs blocks
% are the extracted core functions
%   bvt.samplers.eq_gauss      (per-equation coefficient draw),
%   bvt.sv.ksc_rw_h0           (KSC auxiliary-mixture SV draw, legacy SVRW),
%   bvt.samplers.gig_shrinkage (kappa/psi GIG ladder, variants mng/ng/minn),
%   bvt.samplers.nu_psi_ng     (normal-gamma shape MH step, legacy sample_nu_psi),
% called in the legacy order: theta -> h -> kappa(/psi) -> nu_psi -> h0 -> Sigh.
% The h0 and Sigh draws (identical 3-and-2-line blocks in all three legacy
% scripts) remain inline below.
%
% Output struct: posterior means beta_hat, alp_hat, h_hat, kappa_hat (and
% nu_psi_hat for MNG/NG), the raw stores store_beta/store_alp/store_h/
% store_Sigh/store_kappa (and store_nu_psi, count_nu_psi for MNG/NG), plus
% sig2, dimensions, settings and the preset used.
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3), 1212-1226.

function out = run_all(model, nsim, burnin, seed)
thisdir = fileparts(mfilename('fullpath'));

    % make bvt.* and gigrnd resolvable when called standalone
if isempty(which('bvt.priors.vtheta')) || isempty(which('gigrnd'))
    root = fileparts(fileparts(thisdir));
    addpath(fullfile(root, 'core'));
    addpath(fullfile(root, 'third_party'));
end

    % constants: preset.m in THIS folder (cd guard pins name resolution)
od = cd(thisdir);
guard = onCleanup(@() cd(od));
pr = preset();
clear guard

if nargin < 1 || isempty(model),  model  = 'MNG';               end
if nargin < 2 || isempty(nsim),   nsim   = pr.nsim_default;     end
if nargin < 3 || isempty(burnin), burnin = pr.burnin_default;   end
if nargin < 4, seed = []; end
names = {'MNG', 'NG', 'Minn'};
im = find(strcmpi(model, names), 1);
if isempty(im)
    error('run_all:model', 'model must be ''MNG'', ''NG'' or ''Minn''');
end
model = names{im};
if ~isempty(seed)
    rng(seed, 'twister');
end

    % data (legacy folder, read-only) and design [main_BVAR.m lines 21-35]
raw = load(fullfile(thisdir, 'legacy', pr.data_file));
data = raw(1:pr.data_rows, :);
Y0 = data(1:pr.n0, pr.var_id);
Y  = data(pr.n0+1:end, pr.var_id);
[T, n] = size(Y);
p = pr.p;
[~, Z] = bvt.util.build_lags([Y0(end-p+1:end, :); Y], p);   % identical to the legacy inline construction
k_beta = n^2*p + n;                                         % main_BVAR.m line 36
k_alp  = n*(n-1)/2;                                         % main_BVAR.m line 37
sig2 = bvt.priors.resid_var_ar4(Y0, Y);                     % legacy get_resid_var, main_BVAR.m line 38
[C, idx_kappa1, idx_kappa2] = bvt.priors.minnesota_C(n, p, sig2);   % legacy get_C, main_BVAR.m line 39

    % priors [main_BVAR.m lines 42-46, via preset]
ah = pr.ah; Vh = pr.Vh;
nuh0 = pr.nuh0; Sh0 = pr.Sh0;
c01 = pr.c01; c02 = pr.c02;
lam0_nu_psi = pr.lam0_nu_psi;

    % initialize storage [BVAR_MNG.m lines 9-14 / BVAR_NG.m lines 8-13 / BVAR_Minn.m lines 16-20]
is_ng_type = any(strcmp(model, {'MNG', 'NG'}));
if strcmp(model, 'NG')
    store_kappa = zeros(nsim, 1);
else
    store_kappa = zeros(nsim, 2);
end
if is_ng_type
    store_nu_psi = zeros(nsim, 1);
end
store_beta = zeros(nsim, k_beta);
store_alp = zeros(nsim, k_alp);
store_Sigh = zeros(nsim, n);
store_h = zeros(nsim, T, n);

    % initialize the Markov chain
    % [BVAR_MNG.m lines 17-28 / BVAR_NG.m lines 16-27 / BVAR_Minn.m lines 9-14]
switch model
    case 'MNG'
        nu_psi = pr.mng.nu_psi_init;
        kappa = pr.mng.kappa_init;
    case 'NG'
        nu_psi = pr.ng.nu_psi_init;
        kappa = pr.ng.kappa_init;
    case 'Minn'
        kappa = pr.minn.kappa_init;
end
h0 = log(sig2);
h = repmat(h0', T, 1);
Sigh = Sh0;
if is_ng_type
        % first rng consumption of the run, exactly as in the legacy scripts
    psi_kappa1 = gamrnd(nu_psi, 2/nu_psi, n*p, 1);
    psi_kappa2 = gamrnd(nu_psi, 2/nu_psi, (n-1)*n*p, 1);
    Psi = ones(k_beta, 1);
    Psi(idx_kappa1) = psi_kappa1;
    Psi(idx_kappa2) = psi_kappa2;
    count_nu_psi = 0;
end

    % MCMC starts here (legacy clock-seed line deliberately dropped - see header)
fprintf('Starting MCMC for BVAR-%s.... \n', model);
for isim = 1:nsim + burnin
        % sample alp and beta
    switch model
        case 'MNG'      % BVAR_MNG.m lines 38-39 (incl. the *2 rescaling)
            [Valp, Vbeta] = bvt.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C.*Psi, sig2);
            Valp = pr.mng.vtheta_scale*Valp; Vbeta = pr.mng.vtheta_scale*Vbeta;
        case 'NG'       % BVAR_NG.m line 37 (Psi, not C.*Psi; no rescaling)
            [Valp, Vbeta] = bvt.priors.vtheta(idx_kappa1, idx_kappa2, [kappa, kappa, pr.ng.kappa3, pr.ng.kappa4], Psi, sig2);
        case 'Minn'     % BVAR_Minn.m line 30
            [Valp, Vbeta] = bvt.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C, sig2);
    end
    [beta, alp, U] = bvt.samplers.eq_gauss(Y, Z, h, Valp, Vbeta);

        % sample h [BVAR_MNG.m lines 62-65]
    for ij = 1:n
        Ystar = log(U(:, ij).^2 + pr.sv_offset);
        h(:, ij) = bvt.sv.ksc_rw_h0(Ystar, h(:, ij), Sigh(ij), h0(ij));
    end

        % sample kappa (and, for MNG/NG, psi then nu_psi)
    switch model
        case 'MNG'      % BVAR_MNG.m lines 68-87
            [kappa, psi_kappa1, psi_kappa2] = bvt.samplers.gig_shrinkage('mng', ...
                beta, idx_kappa1, idx_kappa2, C, kappa, psi_kappa1, psi_kappa2, ...
                nu_psi, c01, c02, n, p, pr.psi_floor);
            Psi(idx_kappa1) = psi_kappa1;
            Psi(idx_kappa2) = psi_kappa2;
            [nu_psi, flag] = bvt.samplers.nu_psi_ng(psi_kappa1, psi_kappa2, nu_psi, lam0_nu_psi);
            count_nu_psi = count_nu_psi + flag;
        case 'NG'       % BVAR_NG.m lines 66-84
            [kappa, psi_kappa1, psi_kappa2] = bvt.samplers.gig_shrinkage('ng', ...
                beta, idx_kappa1, idx_kappa2, C, kappa, psi_kappa1, psi_kappa2, ...
                nu_psi, c01, c02, n, p, pr.psi_floor);
            Psi(idx_kappa1) = psi_kappa1;
            Psi(idx_kappa2) = psi_kappa2;
            [nu_psi, flag] = bvt.samplers.nu_psi_ng(psi_kappa1, psi_kappa2, nu_psi, lam0_nu_psi);
            count_nu_psi = count_nu_psi + flag;
        case 'Minn'     % BVAR_Minn.m lines 59-62
            kappa = bvt.samplers.gig_shrinkage('minn', ...
                beta, idx_kappa1, idx_kappa2, C, kappa, [], [], [], c01, c02, n, p, []);
    end

        % sample h0 [BVAR_MNG.m lines 90-92; identical in BVAR_NG/BVAR_Minn]
    Kh0 = sparse(1:n, 1:n, 1./Sigh + 1./Vh);
    h0_hat = Kh0\(ah./Vh + h(1, :)'./Sigh);
    h0 = h0_hat + chol(Kh0, 'lower')'\randn(n, 1);

        % sample Sigh [BVAR_MNG.m lines 95-96; identical in BVAR_NG/BVAR_Minn]
    e = h - [h0'; h(1:T-1, :)];
    Sigh = 1./gamrnd(nuh0 + T/2, 1./(Sh0 + sum(e.^2)'/2));

    if isim > burnin
        isave = isim - burnin;
        switch model
            case 'MNG'      % BVAR_MNG.m lines 100-101
                store_kappa(isave, :) = kappa(1:2)';
                store_nu_psi(isave, :) = nu_psi;
            case 'NG'       % BVAR_NG.m lines 97-98
                store_kappa(isave, :) = kappa;
                store_nu_psi(isave, :) = nu_psi;
            case 'Minn'     % BVAR_Minn.m line 75
                store_kappa(isave, :) = kappa(1:2);
        end
        store_Sigh(isave, :) = Sigh';
        store_beta(isave, :) = beta';
        store_alp(isave, :) = alp';
        store_h(isave, :, :) = h;
    end

    if (mod(isim, 5000) == 0)
        disp([num2str(isim) ' loops... '])
    end
end

    % posterior means [BVAR_MNG.m lines 116-119; identical in BVAR_NG/BVAR_Minn]
h_hat = squeeze(mean(store_h));
alp_hat = mean(store_alp)';
beta_hat = mean(store_beta)';
kappa_hat = mean(store_kappa)';

out = struct();
out.model = model;
out.nsim = nsim;
out.burnin = burnin;
out.seed = seed;
out.T = T; out.n = n; out.p = p;
out.sig2 = sig2;
out.beta_hat = beta_hat;
out.alp_hat = alp_hat;
out.h_hat = h_hat;
out.kappa_hat = kappa_hat;
out.store_beta = store_beta;
out.store_alp = store_alp;
out.store_h = store_h;
out.store_Sigh = store_Sigh;
out.store_kappa = store_kappa;
if is_ng_type
    out.nu_psi_hat = mean(store_nu_psi);
    out.store_nu_psi = store_nu_psi;
    out.count_nu_psi = count_nu_psi;
end
out.preset = pr;
end
