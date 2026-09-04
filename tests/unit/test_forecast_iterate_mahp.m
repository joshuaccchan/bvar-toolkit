function test_forecast_iterate_mahp
% seeded draw-for-draw equivalence of the functionized MAHP FORECAST pipeline
% (bvar.priors.* + bvar.samplers.* + bvar.sv.ksc_rw_h0 for the estimation stage,
% bvar.forecast.iterate('mahp_sv',...) for the per-draw forecast stage) with the
% legacy workspace script forecast_BVAR_MNG.m, primed exactly as
% chan2021_ijf_mahp main_forecasting.m primes it for one vintage t and run at
% small nsim from a tempdir copy. Asserts isequal on every forecast quantity
% (tmpyhat1, tmpyhat4, all MCMC stores, kappa_hat, kappaCI) and on the
% terminal rng state (same randn/rand/gamrnd call sequence).
%
% ZERO PATCHES: unlike the estimation scripts (BVAR_MNG.m etc., whose sole
% clock-seed line test_mahp_equivalence removes), the forecast scripts carry
% NO clock-seed line - verified below - so the tempdir copy runs byte-verbatim
% and the harness controls seeding via rng(seed,'twister') before dispatch.
% main_forecasting.m's per-vintage setup (lines 18-46 and 65-74) is replicated
% in run_legacy below, which is also where the hard-coded nsim = 20000 /
% burnin = 100 are overridden (the script reads them from the workspace).
%
% Two vintages: t = 91 (the legacy T0; both horizons evaluated) and t = T-2
% (the tt==4 guard is OFF: tmpyhat4 must stay all-zero on both sides while the
% simulation still consumes the same rng draws).
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2021_ijf_mahp', 'legacy');

tmpdir = tempname; mkdir(tmpdir);
ctmp = onCleanup(@() cleanup_tmp(tmpdir));   %#ok<NASGU> % rmpath BEFORE rmdir
files = {'forecast_BVAR_MNG.m', 'SVRW.m', 'getVtheta.m', 'get_C.m', ...
    'get_resid_var.m', 'gigrnd.m', 'sample_nu_psi.m'};
for kf = 1:numel(files)
    copyfile(fullfile(leg, files{kf}), fullfile(tmpdir, files{kf}));
end
assert(isempty(strfind(fileread(fullfile(tmpdir, 'forecast_BVAR_MNG.m')), 'randn(''seed''')), ...
    'unexpected clock-seed line in forecast_BVAR_MNG.m - the zero-patch premise is broken'); %#ok<STREMP>
addpath(tmpdir);

    % shared data/design (main_forecasting.m lines 21-33), legacy csv read-only
p = 4;
raw = load(fullfile(leg, 'macrodata_Q_2018Q4.csv'));
data = raw(1:238, :);
var_id = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,160,161,245];
Y0 = data(1:8, var_id);
Y = data(9:end, var_id);
T = size(Y, 1);

    % {vintage t, nsim, burnin}
cases = {91, 20, 10; T-2, 6, 3};
seed = 20260901;
for kc = 1:size(cases, 1)
    [t, nsim, burnin] = cases{kc, :};
    L = run_legacy(tmpdir, t, nsim, burnin, seed, p, Y0, Y);
    R = run_core(t, nsim, burnin, seed, p, Y0, Y);
    lbl = sprintf('vintage t=%d', t);
    assert(isequal(L.store_beta,  R.store_beta),  '%s: store_beta differs',  lbl);
    assert(isequal(L.store_alp,   R.store_alp),   '%s: store_alp differs',   lbl);
    assert(isequal(L.store_h_T,   R.store_h_T),   '%s: store_h_T differs',   lbl);
    assert(isequal(L.store_Sigh,  R.store_Sigh),  '%s: store_Sigh differs',  lbl);
    assert(isequal(L.store_kappa, R.store_kappa), '%s: store_kappa differs', lbl);
    assert(isequal(L.kappa_hat, R.kappa_hat) && isequal(L.kappaCI, R.kappaCI), ...
        '%s: kappa_hat/kappaCI differ', lbl);
    assert(isequal(L.tmpyhat1, R.tmpyhat1), '%s: tmpyhat1 differs', lbl);
    assert(isequal(L.tmpyhat4, R.tmpyhat4), '%s: tmpyhat4 differs', lbl);
    assert(isequal(L.rngstate, R.rngstate), '%s: rng call sequence differs', lbl);
    if t > T-4
        assert(~any(L.tmpyhat4(:)) && ~any(R.tmpyhat4(:)), ...
            '%s: tmpyhat4 expected all-zero (tt==4 guard off)', lbl);
    end
end
end

function cleanup_tmp(tmpdir)
if any(strcmpi(strsplit(path, pathsep), tmpdir))
    rmpath(tmpdir);
end
if exist(tmpdir, 'dir')
    rmdir(tmpdir, 's');
end
end

function out = run_legacy(tmpdir, t, nsim, burnin, seed, p, Y0, Y) %#ok<INUSD,INUSL> % nsim/burnin/p/Y0/Y are read by the dispatched script from this workspace
% replicate main_forecasting.m's per-vintage setup (lines 33-46, 65-74) with
% the small nsim/burnin already assigned, then dispatch the byte-verbatim
% forecast script exactly as main_forecasting.m line 77 does.
[T, n] = size(Y);                                       %#ok<ASGLU>
ah = zeros(n, 1); Vh = 10*ones(n, 1);                   %#ok<NASGU>
nuh0 = 5*ones(n, 1); Sh0 = .01*ones(n, 1).*(nuh0-1);    %#ok<NASGU>
c01 = [1, 1/.04];                                       %#ok<NASGU>
c02 = [1, 1/.04^2];                                     %#ok<NASGU>
lam0_nu_psi = 1;                                        %#ok<NASGU>
Yt = Y(1:t, :);
Tt = size(Yt, 1);                                       %#ok<NASGU>
tmpYt = [Y0(end-p+1:end, :); Yt];
Zt = zeros(t, n*p);
for ii = 1:p
    Zt(:, (ii-1)*n+1:ii*n) = tmpYt(p-ii+1:end-ii, :);
end
Zt = [ones(t, 1) Zt];                                   %#ok<NASGU>

    % pre-declare script-assigned variables read back below (parser binding)
tmpyhat1 = []; tmpyhat4 = [];
store_beta = []; store_alp = []; store_h_T = []; store_Sigh = []; store_kappa = [];
kappa_hat = []; kappaCI = [];

resolved = which('forecast_BVAR_MNG');
assert(strncmpi(resolved, tmpdir, numel(tmpdir)), ...
    'forecast_BVAR_MNG must resolve from the tempdir copy, got %s', resolved);

rng(seed, 'twister');
forecast_BVAR_MNG;
s = rng;

out = struct('tmpyhat1', tmpyhat1, 'tmpyhat4', tmpyhat4, ...
    'store_beta', store_beta, 'store_alp', store_alp, 'store_h_T', store_h_T, ...
    'store_Sigh', store_Sigh, 'store_kappa', store_kappa, ...
    'kappa_hat', kappa_hat, 'kappaCI', kappaCI, 'rngstate', s.State);
end

function out = run_core(t, nsim, burnin, seed, p, Y0, Y)
% the functionized path: estimation stage from the extracted core blocks in
% the legacy order (constants verbatim from forecast_BVAR_MNG.m, incl. its
% estimation-vs-forecast divergences recorded in preset.m pr.forecast: psi
% floor 1e-16, halved psi_kappa1 init scale, doubled Valp/Vbeta), forecast
% stage via bvar.forecast.iterate('mahp_sv',...) once per kept draw.
[T, n] = size(Y);
ah = zeros(n, 1); Vh = 10*ones(n, 1);
nuh0 = 5*ones(n, 1); Sh0 = .01*ones(n, 1).*(nuh0-1);
c01 = [1, 1/.04];
c02 = [1, 1/.04^2];
lam0_nu_psi = 1;
Yt = Y(1:t, :);
Tt = size(Yt, 1);
[~, Zt] = bvar.util.build_lags([Y0(end-p+1:end, :); Yt], p);   % identical to the legacy inline construction

rng(seed, 'twister');
    % ---- estimation stage [forecast_BVAR_MNG.m lines 8-109] ----
tmpyhat1 = zeros(nsim, 2*n+1);
tmpyhat4 = zeros(nsim, 2*n+1);
k_beta = n^2*p+n;
k_alp = n*(n-1)/2;                                      %#ok<NASGU>
sig2 = bvar.priors.resid_var_ar4(Y0, Yt);                % legacy get_resid_var
[C, idx_kappa1, idx_kappa2] = bvar.priors.minnesota_C(n, p, sig2);   % legacy get_C

store_kappa = zeros(nsim, 2);
store_alp = zeros(nsim, n*(n-1)/2);
store_beta = zeros(nsim, k_beta);
store_h_T = zeros(nsim, n);
store_Sigh = zeros(nsim, n);

nu_psi = .5;
kappa = [.4, .001, 1, 100];
h0 = log(sig2);
h = repmat(h0', Tt, 1);
Sigh = Sh0;
psi_kappa1 = gamrnd(nu_psi, 2/nu_psi/2, n*p, 1);        % forecast_BVAR_MNG.m line 31: HALVED init scale
psi_kappa2 = gamrnd(nu_psi, 2/nu_psi, (n-1)*n*p, 1);
Psi = ones(k_beta, 1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

for isim = 1:nsim + burnin
        % sample alp and beta [lines 41-61]
    [Valp, Vbeta] = bvar.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C.*Psi, sig2);
    Valp = 2*Valp; Vbeta = 2*Vbeta;
    [beta, alp, U] = bvar.samplers.eq_gauss(Yt, Zt, h, Valp, Vbeta);
        % sample h [lines 64-67]
    for ij = 1:n
        Ystar = log(U(:, ij).^2 + .0001);
        h(:, ij) = bvar.sv.ksc_rw_h0(Ystar, h(:, ij), Sigh(ij), h0(ij));
    end
        % sample kappa1/kappa2 and psi [lines 70-85; forecast psi floor 1e-16]
    [kappa, psi_kappa1, psi_kappa2] = bvar.samplers.gig_shrinkage('mng', ...
        beta, idx_kappa1, idx_kappa2, C, kappa, psi_kappa1, psi_kappa2, ...
        nu_psi, c01, c02, n, p, 1e-16);
    Psi(idx_kappa1) = psi_kappa1;
    Psi(idx_kappa2) = psi_kappa2;
        % sample nu_psi [line 88; one-output call, flag discarded]
    nu_psi = bvar.samplers.nu_psi_ng(psi_kappa1, psi_kappa2, nu_psi, lam0_nu_psi);
        % sample h0 [lines 91-93]
    Kh0 = sparse(1:n, 1:n, 1./Sigh + 1./Vh);
    h0_hat = Kh0\(ah./Vh + h(1, :)'./Sigh);
    h0 = h0_hat + chol(Kh0, 'lower')'\randn(n, 1);
        % sample Sigh [lines 96-97]
    e = h - [h0'; h(1:Tt-1, :)];
    Sigh = 1./gamrnd(nuh0+Tt/2, 1./(Sh0 + sum(e.^2)'/2));

    if isim > burnin
        isave = isim - burnin;
        store_kappa(isave, :) = kappa(1:2);
        store_beta(isave, :) = beta';
        store_alp(isave, :) = alp';
        store_h_T(isave, :) = h(end, :)';
        store_Sigh(isave, :) = Sigh';
    end
end
kappa_hat = mean(store_kappa)';
kappaCI = quantile(store_kappa, [0.25 .975]);

    % ---- forecast stage [lines 112-147] via bvar.forecast.iterate ----
cfg = struct('Yt', Yt, 'Y', Y, 'p', p, 't', t, 'T', T);
for isim = 1:nsim
    draw = struct('alp', store_alp(isim, :)', 'beta', store_beta(isim, :)', ...
        'h_T', store_h_T(isim, :)', 'Sigh', store_Sigh(isim, :)');
    fcr = bvar.forecast.iterate('mahp_sv', draw, cfg);
    tmpyhat1(isim, :) = fcr(1, :);
    tmpyhat4(isim, :) = fcr(2, :);
end
s = rng;

out = struct('tmpyhat1', tmpyhat1, 'tmpyhat4', tmpyhat4, ...
    'store_beta', store_beta, 'store_alp', store_alp, 'store_h_T', store_h_T, ...
    'store_Sigh', store_Sigh, 'store_kappa', store_kappa, ...
    'kappa_hat', kappa_hat, 'kappaCI', kappaCI, 'rngstate', s.State);
end
