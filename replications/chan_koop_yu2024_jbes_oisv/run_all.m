% chan_koop_yu2024_jbes_oisv/run_all - functionized full-sample estimation
% pipeline of Chan, Koop and Yu (2024, JBES): large order-invariant ('OI') or
% Cholesky/triangularized ('CS') Bayesian SVAR with stochastic volatility and
% the Minnesota-type horseshoe prior, at the default or reversed variable
% ordering - the four configurations run by legacy main_SVAR_fullsample.m
% (CS1, CS2, OI1, OI2).
%
%   out = run_all(model, flip, nsim, burnin, seed)
%
%   model  - 'OI' (legacy model 1, SVARSV_MH.m) | 'CS' (legacy model 2,
%            CS_MH.m); default 'OI' (the commented legacy default
%            func_main_SVAR_v2.m line 11: model = 1)
%   flip   - 0: default ordering | 1: reversed ordering (legacy rev_option,
%            func_main_SVAR_v2.m lines 19-21); default 0
%   nsim   - posterior draws kept (default preset nsim_default = 30000)
%   burnin - burn-in sweeps (default preset burnin_default = 5000)
%   seed   - optional; when nonempty, rng(seed,'twister') is set FIRST, before
%            any draw. When omitted/empty the ambient rng state is used as-is.
%
% Functionized 2026-09-02 (step 7, OISV family pass). Reproduces the legacy
% pipeline main_SVAR_fullsample.m -> func_main_SVAR_v2.m -> SVARSV_MH.m /
% CS_MH.m draw-for-draw bitwise (verified by tests/unit/test_oisv_equivalence.m
% at small nsim) with ONE deliberate divergence: legacy SVARSV_MH.m line 24
% re-seeds the global stream from the wall clock (randn('seed',sum(clock*100));
% rand('seed',sum(clock*1000))), which makes as-shipped OI runs irreproducible
% AND switches MATLAB to the legacy v4/v5 generators; run_all drops that line
% so the caller controls seeding via `seed` (or the ambient state) on the
% modern twister stream. CS_MH.m carries the same line only COMMENTED OUT
% (line 49) - the legacy CS run already consumes the ambient stream. In both
% samplers every rng draw sits after that point, so the whole run (chain-init
% draws included) consumes one coherent stream. The legacy wall-clock timing
% displays (start_time/etime) are not reproduced.
%
% All constants come from preset.m in this folder (each field cites its legacy
% source line); the data file is read from legacy/ READ-ONLY; the Gibbs blocks
% are the extracted core functions
%   bvt.structural.b0_row_sampler   (OI row-wise B0 rotation draw, using
%                                    bvt.util.anormrnd for the first coordinate),
%   bvt.samplers.eq_svar_oi         (OI per-equation VAR coefficient draw),
%   bvt.samplers.eq_tri_cs          (CS per-equation VAR coefficient draw),
%   bvt.samplers.alp_tri_cs         (CS free impact-element draw),
%   bvt.samplers.horseshoe_kappa_psi (horseshoe psi/z_psi/kappa/z_kappa block,
%                                    shared verbatim by both legacy samplers),
% with the pre-existing core reused: bvt.priors.resid_var_ar4 (legacy
% get_resid_var), bvt.priors.minnesota_C (legacy get_C), bvt.priors.vtheta
% whose Vbeta output reproduces legacy getVbeta exactly (same three assignment
% lines on the same inputs; the Valp output is discarded - NaN under the OI
% kappa(3) = NaN, never read), bvt.sv.ksc_ar1_mean (legacy sample_SV),
% bvt.sv.sv0_params (legacy sample_SV0para, OI) and bvt.sv.sv_params (legacy
% sample_SVpara, CS) at their OISV-canonical default phi bounds (.99 / .999),
% bvt.util.build_lags (the inline lag construction), bvt.util.vec, and
% bvt.structural.construct_Sigt (legacy utility/construct_Sigt.m = the private
% duplicate inside func_main_SVAR_v2.m) for the posterior covariance paths.
% Blocks called in the legacy order - OI: B0 -> alpha -> h -> (phi,sig2) ->
% horseshoe; CS: B -> alp -> h -> (mu,phi,sig2) -> horseshoe.
%
% Output struct: the six func_main_SVAR_v2 outputs (Mean_kappa, Mean_Impact,
% Median_kappa, Median_Impact, Sig_mean, Sig_median - Sig_median is zeros(T,n,n)
% exactly as in the legacy func, lines 46/62), the raw stores (OI: store_kappa/
% store_hpara/store_alpha/store_B0/store_h; CS: store_beta/store_alp/store_h/
% store_hpara/store_kappa), count_phi, the legacy script-tail posterior
% summaries, plus dimensions, settings and the preset used.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function out = run_all(model, flip, nsim, burnin, seed)
thisdir = fileparts(mfilename('fullpath'));

    % make bvt.* resolvable when called standalone
if isempty(which('bvt.structural.b0_row_sampler'))
    root = fileparts(fileparts(thisdir));
    addpath(fullfile(root, 'core'));
end

    % constants: preset.m in THIS folder (cd guard pins name resolution)
od = cd(thisdir);
guard = onCleanup(@() cd(od));
pr = preset();
clear guard

if nargin < 1 || isempty(model),  model  = 'OI';                end
if nargin < 2 || isempty(flip),   flip   = 0;                   end
if nargin < 3 || isempty(nsim),   nsim   = pr.nsim_default;     end
if nargin < 4 || isempty(burnin), burnin = pr.burnin_default;   end
if nargin < 5, seed = []; end
names = {'OI', 'CS'};
im = find(strcmpi(model, names), 1);
if isempty(im)
    error('run_all:model', 'model must be ''OI'' or ''CS''');
end
model = names{im};
if ~isempty(seed)
    rng(seed, 'twister');
end

    % data (legacy folder, read-only) and design [func_main_SVAR_v2.m lines 8-31]
p = pr.p;
data = load(fullfile(thisdir, 'legacy', pr.data_file));
var_id = pr.var_id;
if flip                                 % rev_option == 1 [func_main_SVAR_v2.m lines 19-21]
    var_id = var_id(end:-1:1);          % = flip(var_id); spelled out because the argument shadows flip()
end
Y0 = data(1:pr.n0, var_id);             % save the first 24 obs as the initial conditions
Y  = data(pr.n0+1:end, var_id);
[T, n] = size(Y);
[~, X] = bvt.util.build_lags([Y0(end-p+1:end, :); Y], p);   % identical to the legacy inline construction

    % MCMC [func_main_SVAR_v2.m lines 33-63]
switch model
    case 'OI'
        res = mcmc_oi(Y, X, Y0, T, n, p, nsim, burnin, pr);
            % post-processing [func_main_SVAR_v2.m lines 36-46]
        Mean_kappa = res.kappa_mean; Median_kappa = median(res.store_kappa)';
        Mean_Impact = reshape(res.B0_mean, n, n)';
        B0_median = median(res.store_B0)'; Median_Impact = reshape(B0_median, n, n)';

        store_Sig = zeros(T, n, n);
        for isim = 1:nsim
            B0 = reshape(res.store_B0(isim, :), n, n)'; % stacked by rows
            h = squeeze(res.store_h(isim, :, :));
            store_Sig = store_Sig + bvt.structural.construct_Sigt(h, B0);
        end
        Sig_mean = store_Sig/nsim; Sig_median = zeros(T, n, n);
    case 'CS'
        res = mcmc_cs(Y, X, Y0, T, n, p, nsim, burnin, pr);
            % post-processing [func_main_SVAR_v2.m lines 49-62]
        Mean_kappa = res.kappa_mean; Median_kappa = median(res.store_kappa)';
        A_id = nonzeros(tril(reshape(1:n^2, n, n), -1)'); A = eye(n); A(A_id) = res.alp_mean;
        Mean_Impact = A;
        A_id = nonzeros(tril(reshape(1:n^2, n, n), -1)'); A = eye(n); alp_median = median(res.store_alp)'; A(A_id) = alp_median;
        Median_Impact = A;

        store_Sig = zeros(T, n, n);
        for isim = 1:nsim
            h = squeeze(res.store_h(isim, :, :));
            alp = squeeze(res.store_alp(isim, :))';
            A = eye(n); A(A_id) = alp;
            store_Sig = store_Sig + bvt.structural.construct_Sigt(h, A);
        end
        Sig_mean = store_Sig/nsim; Sig_median = zeros(T, n, n);
end

out = res;
out.model = model;
out.flip = flip;
out.nsim = nsim;
out.burnin = burnin;
out.seed = seed;
out.T = T; out.n = n; out.p = p;
out.var_id = var_id;
out.Mean_kappa = Mean_kappa;
out.Median_kappa = Median_kappa;
out.Mean_Impact = Mean_Impact;
out.Median_Impact = Median_Impact;
out.Sig_mean = Sig_mean;
out.Sig_median = Sig_median;
out.preset = pr;
end

% -------------------------------------------------------------------------
function res = mcmc_oi(Y, X, Y0, T, n, p, nsim, burnin, pr)
% SVARSV_MH.m functionized line-for-line (the clock-seed line 24 dropped - see
% the run_all header). Draw order per sweep: B0 -> alpha -> h -> (phi,sig2) ->
% horseshoe block.
k = 1+n*p;                                              % SVARSV_MH.m line 3
k_alpha = n*k;                                          % line 4
sig2 = bvt.priors.resid_var_ar4(Y0, Y);                 % line 5 (legacy get_resid_var)
[C, idx_kappa1, idx_kappa2] = bvt.priors.minnesota_C(n, p, sig2);   % line 6 (legacy get_C)

    % priors [SVARSV_MH.m lines 9-11, via preset]
Hyper.nuh = pr.oi.nuh; Hyper.Sh = pr.oi.Sh;
Hyper.phi0 = pr.oi.phi0; Hyper.Vphi = pr.oi.Vphi;
Hyper.B0 = pr.oi.B0_prior_mean; Hyper.VB0 = pr.oi.VB0;

    % initialization for storage [lines 14-19]
store_kappa = zeros(nsim, 2);
store_alpha = zeros(nsim, k_alpha);
store_B0 = zeros(nsim, n^2);
store_h = zeros(nsim, T, n);
store_hpara = zeros(nsim, 2*n);
count_phi = zeros(n, 1);

    % (legacy clock-seed line 24 deliberately dropped)
disp('Starting MCMC for BSVAR.... ');

    % initialize the Markov chain [lines 28-44]
kappa = pr.oi.kappa_init;   % [.1,.1,NaN,100]: kappa(1) own lag; kappa(2) other lag; kappa(4) intercepts
A = (X'*X + pr.ls_ridge*speye(k))\(X'*Y);
U = Y-X*A;
alpha = A(:);                                           %#ok<NASGU> % line 31 (recomputed each sweep)
Sig_hat = U'*U/T;
B0 = diag(1./sqrt(diag(Sig_hat)));
h = repmat(log(diag(Sig_hat))', T, 1);
sig2 = 1./gamrnd(Hyper.nuh, 1./Hyper.Sh);
phi = min(Hyper.phi0 + sqrt(Hyper.Vphi).*randn(n, 1), pr.phi_init_bnd);
E = zeros(T, n);                                        %#ok<PREALL> % line 37 (recomputed each sweep)
z_psi1 = 1./gamrnd(.5, 1, n*p, 1);
z_psi2 = 1./gamrnd(.5, 1, (n-1)*n*p, 1);
z_kappa = 1./gamrnd(.5, 1, 2, 1);
psi_kappa1 = 1./gamrnd(.5, z_psi1);
psi_kappa2 = 1./gamrnd(.5, z_psi2);
Psi = ones(k*n, 1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

for isim = 1:nsim + burnin
        % sample B0 [lines 48-72 -> bvt.structural.b0_row_sampler]
    U = Y-X*A;
    B0 = bvt.structural.b0_row_sampler(U, h, B0, Hyper.B0, Hyper.VB0);

        % sample alpha [lines 75-88 -> bvt.samplers.eq_svar_oi;
        % legacy getVbeta = the Vbeta output of bvt.priors.vtheta]
    [~, tmpdV] = bvt.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C.*Psi, sig2);
    A = bvt.samplers.eq_svar_oi(Y, X, B0, h, A, tmpdV);
    alpha = A(:);

        % sample h [lines 90-95]
    E = (Y - X*A)*B0';
    for ii = 1:n
        ystar = log(E(:, ii).^2 + pr.sv_offset);
        h(:, ii) = bvt.sv.ksc_ar1_mean(ystar, h(:, ii), pr.oi.h_mean_in_sv, phi(ii), sig2(ii));
    end

        % sample phi and sig2 [lines 97-99; default phi bound = OISV .99]
    [phi, sig2, flag_phi] = bvt.sv.sv0_params(h, phi, Hyper);
    count_phi = count_phi + flag_phi;

        % sample psi, z_psi, kappa1/kappa2, z_kappa
        % [lines 101-120 -> bvt.samplers.horseshoe_kappa_psi, theta = alpha]
    [psi_kappa1, psi_kappa2, z_psi1, z_psi2, kappa, z_kappa] = ...
        bvt.samplers.horseshoe_kappa_psi(alpha, idx_kappa1, idx_kappa2, C, kappa, z_psi1, z_psi2, z_kappa);
    Psi(idx_kappa1) = psi_kappa1;
    Psi(idx_kappa2) = psi_kappa2;

    if isim > burnin                                    % [lines 122-129]
        isave = isim-burnin;
        store_kappa(isave, :) = kappa(1:2);
        store_hpara(isave, :) = [phi', sig2'];
        store_alpha(isave, :) =  alpha';
        store_B0(isave, :) =  reshape(B0', n^2, 1); % stacked by rows
        store_h(isave, :, :) =  h;
    end

    if (mod(isim, pr.progress_every) == 0)              % [lines 131-133]
        disp([num2str(isim) ' loops... '])
    end
end

    % posterior summaries [SVARSV_MH.m lines 139-144]
res = struct();
res.h_mean = squeeze(mean(store_h));
res.B0_mean = mean(store_B0)'; res.B0_median = median(store_B0)';
res.alpha_mean = mean(store_alpha)';
res.hpara_mean = mean(store_hpara)';
res.kappa_mean = mean(store_kappa)'; res.kappa_median = median(store_kappa)';
res.A_mean = reshape(res.alpha_mean, k, n)';
res.store_kappa = store_kappa;
res.store_hpara = store_hpara;
res.store_alpha = store_alpha;
res.store_B0 = store_B0;
res.store_h = store_h;
res.count_phi = count_phi;
end

% -------------------------------------------------------------------------
function res = mcmc_cs(Y, X, Y0, T, n, p, nsim, burnin, pr)
% CS_MH.m functionized line-for-line (its clock-seed line 49 is already
% commented out in the legacy script - nothing to drop). Draw order per sweep:
% B -> alp -> h -> (mu,phi,sig2) -> horseshoe block.
k_alp = n*(n-1)/2;                                      % CS_MH.m line 3
k_beta = n^2*p + n;                                     % line 4
k = k_beta/n;                                           % line 5

    % priors [lines 8-16, via preset]
sig2 = bvt.priors.resid_var_ar4(Y0, Y);                 % line 8 (legacy get_resid_var)
kappa = pr.cs.kappa_init;   % [.1,.1,1,100]: kappa(1) own lag; kappa(2) other lag; kappa(4) intercepts; kappa(3) alpha
[C, idx_kappa1, idx_kappa2] = bvt.priors.minnesota_C(n, p, sig2);   % line 10 (legacy get_C)
Hyper.beta0 = pr.cs.beta0;
Hyper.alp0 = pr.cs.alp0; Hyper.Valp = pr.cs.Valp;
Hyper.mu0 = pr.cs.mu0; Hyper.Vmu = pr.cs.Vmu;
Hyper.nuh = pr.cs.nuh; Hyper.Sh = pr.cs.Sh;
Hyper.phi0 = pr.cs.phi0; Hyper.Vphi = pr.cs.Vphi;

    % initialize storage [lines 19-24]
store_alp = zeros(nsim, k_alp);
store_beta = zeros(nsim, k_beta);
store_h = zeros(nsim, T, n);
store_hpara = zeros(nsim, 3*n);
store_kappa = zeros(nsim, 2);
count_phi = zeros(n, 1);

A_id = nonzeros(tril(reshape(1:n^2, n, n), -1)');       % line 26
A = eye(n);                                             % line 27

    % initialize the Markov chain [lines 30-46]
sig2 = 1./gamrnd(Hyper.nuh, 1./Hyper.Sh);
phi = Hyper.phi0; phi = min(Hyper.phi0 + sqrt(Hyper.Vphi).*randn(n, 1), pr.phi_init_bnd); %#ok<NASGU> % dead first assignment kept verbatim [line 31]
XX = X'*X; B = ((XX + pr.ls_ridge*speye(k))\(X'*Y))';
XB = X*B'; U = Y-XB; Sig_hat = U'*U/T;
mu = zeros(n, 1);
for ii = 1:n
    s2i = U(:, ii).^2; mu(ii) = mean(log(s2i));
end
h = repmat(log(diag(Sig_hat))', T, 1);
z_psi1 = 1./gamrnd(.5, 1, n*p, 1);
z_psi2 = 1./gamrnd(.5, 1, (n-1)*n*p, 1);
z_kappa = 1./gamrnd(.5, 1, 2, 1);
psi_kappa1 = 1./gamrnd(.5, z_psi1);
psi_kappa2 = 1./gamrnd(.5, z_psi2);
Psi = ones(k*n, 1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

for isim = 1:nsim + burnin
        % sample B [lines 53-73 -> bvt.samplers.eq_tri_cs;
        % legacy getVbeta = the Vbeta output of bvt.priors.vtheta]
    [~, tmpdV] = bvt.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C.*Psi, sig2);
    [B, XB] = bvt.samplers.eq_tri_cs(Y, X, XB, B, A, h, tmpdV, Hyper.beta0);
    beta = reshape(B', k_beta, 1);

        % sample alp [lines 76-88 -> bvt.samplers.alp_tri_cs]
    E = Y - XB;
    alp = bvt.samplers.alp_tri_cs(E, h, Hyper.Valp);
    A(A_id) = alp;

        % sample h [lines 90-95]
    AE = E*sparse(A');
    for ii = 1:n
        ystar = log(AE(:, ii).^2 + pr.sv_offset);
        h(:, ii) = bvt.sv.ksc_ar1_mean(ystar, h(:, ii), mu(ii), phi(ii), sig2(ii));
    end

        % sample mu, phi and sig2 [lines 97-99; default phi bound = OISV .999]
    [mu, phi, sig2, flag_phi] = bvt.sv.sv_params(h, mu, phi, Hyper);
    count_phi = count_phi + flag_phi;

        % sample psi, z_psi, kappa1/kappa2, z_kappa
        % [lines 101-120 -> bvt.samplers.horseshoe_kappa_psi, theta = beta]
    [psi_kappa1, psi_kappa2, z_psi1, z_psi2, kappa, z_kappa] = ...
        bvt.samplers.horseshoe_kappa_psi(beta, idx_kappa1, idx_kappa2, C, kappa, z_psi1, z_psi2, z_kappa);
    Psi(idx_kappa1) = psi_kappa1;
    Psi(idx_kappa2) = psi_kappa2;

    if isim > burnin                                    % [lines 122-129]
        isave = isim-burnin;
        store_beta(isave, :) = reshape(B', k_beta, 1);
        store_alp(isave, :) = alp';
        store_h(isave, :, :) = h;
        store_hpara(isave, :) = [mu', phi', sig2'];
        store_kappa(isave, :) = kappa(1:2);
    end

    if (mod(isim, pr.progress_every) == 0)              % [lines 131-133]
        disp([num2str(isim) ' loops... '])
    end
end

    % posterior summaries [CS_MH.m lines 138-142]
res = struct();
res.beta_mean = mean(store_beta)'; res.Beta_mean = reshape(res.beta_mean, k, n)';
res.alp_mean = mean(store_alp)'; res.alp_median = median(store_alp)';
res.h_mean = squeeze(mean(store_h));
res.hpara_mean = mean(store_hpara)';
res.kappa_mean = mean(store_kappa)'; res.kappa_median = median(store_kappa)';
res.store_beta = store_beta;
res.store_alp = store_alp;
res.store_h = store_h;
res.store_hpara = store_hpara;
res.store_kappa = store_kappa;
res.count_phi = count_phi;
end
