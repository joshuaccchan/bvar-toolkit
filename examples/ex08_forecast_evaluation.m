%% ex08 - Recursive forecasting and predictive evaluation
%
% BOOK: Chapter 14, Large VARs with Stochastic Volatility, in Bayesian
% Macroeconometrics: Methods and Applications (Chapman & Hall/CRC, forthcoming).
%
% This example carries a model through to a forecast comparison: estimate on
% data up to a vintage, forecast one and four quarters ahead, score the forecast
% against the outturn, move the vintage forward and repeat. bvar.forecast holds
% the two pieces that do the work, iterate for one draw and tables for the
% accumulation and the summary.
%
% THE MODEL. The large BVAR with stochastic volatility of Chan (2021), written
% in STRUCTURAL form: A0 y_t = b + B_1 y_{t-1} + ... + B_p y_{t-p} + eps_t, with
% A0 unit lower triangular and eps_it ~ N(0, exp(h_it)). Equation i therefore
% regresses y_i on the lags AND on the contemporaneous y_1..y_{i-1}, which is
% why Xi = [Z -Y(:,1:i-1)] appears inside bvar.samplers.eq_gauss. Two blocks of
% coefficients follow: beta, the n^2 p + n lag coefficients, and alp, the
% n(n-1)/2 free elements of A0.
%
% THE PRIOR. The Minnesota-type normal-gamma prior, model 1 of the paper, which
% is normal-gamma in the exact sense: lag coefficient j has
%
%     beta_j | psi_j ~ N(0, kappa * C_j * psi_j),   psi_j ~ Gamma(nu_psi, 2/nu_psi),
%
% so integrating out psi_j leaves heavier tails and more mass at zero than a
% normal. Minnesota-TYPE refers to two things the paper's plain normal-gamma
% (model 2) does without: C, from bvar.priors.minnesota_C, carrying the decay in
% lag order and the residual-variance ratio across equations, and a global
% shrinkage split in two, kappa1 on own lags and kappa2 on cross lags, with Gamma
% priors of mean .04 and .0016. kappa1, kappa2, every psi_j and nu_psi are all
% estimated, where a standard Minnesota prior fixes kappa and sets psi_j to one.
% The generalized inverse Gaussian is the CONDITIONAL POSTERIOR of psi_j and of
% the kappas, not part of the prior; bvar.samplers.gig_shrinkage draws it.
%
% SETTINGS. Twelve variables, twelve vintages and 200 draws after a burn-in of
% 100. The published exercise uses 23 variables, about 139 vintages and 20000
% draws; see replications/chan2021_ijf_mahp. Twelve vintages is far too few for
% the RMSFEs below to rank anything, and they are printed to show the machinery.
%
% DATA. Read-only from replications/chan2021_ijf_mahp/legacy/
% macrodata_Q_2018Q4.csv, the quarterly US panel of Chan (2021), first twelve of
% the paper's twenty-three series.
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for Large
% Bayesian VARs, International Journal of Forecasting, 37(3): 1212-1226.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

repo = fileparts(fileparts(mfilename('fullpath')));
fprintf('\n=== ex08: recursive forecasting and predictive evaluation ===\n');

%% ------------------------------------------------------------------
%  1. Data and sizes
%  ------------------------------------------------------------------
csv = fullfile(repo, 'replications', 'chan2021_ijf_mahp', 'legacy', ...
    'macrodata_Q_2018Q4.csv');
raw  = readmatrix(csv);                       % never written to; legacy is frozen
data = raw(1:238, :);
var_id = [1,2,22,23,35,37,57,58,59,76,81,83];  % first 12 of the paper's 23
p = 4;
Y0 = data(1:8, var_id);
Y  = data(9:end, var_id);
[T, n] = size(Y);

k_beta = n^2*p + n;
k_alp  = n*(n-1)/2;
fprintf('\nT = %d quarters, n = %d variables, p = %d lags\n', T, n, p);
fprintf('%d lag coefficients and %d free elements of A0: %d in all, from %d observations\n', ...
    k_beta, k_alp, k_beta+k_alp, T*n);
fprintf('The shrinkage prior below is what makes that ratio workable.\n');

%% ------------------------------------------------------------------
%  2. Prior constants, verbatim from the paper's forecasting driver
%  ------------------------------------------------------------------
ah = zeros(n,1);  Vh = 10*ones(n,1);
nuh0 = 5*ones(n,1);  Sh0 = .01*ones(n,1).*(nuh0-1);
c01 = [1, 1/.04];        % prior for kappa1, the own-lag shrinkage
c02 = [1, 1/.04^2];      % prior for kappa2, the cross-lag shrinkage
lam0_nu_psi = 1;

nsim = 200; burnin = 100;
nvintage = 12;
t_first = T - nvintage;                  % vintages t_first .. T-1
fprintf('\n%d vintages, %d draws after %d burn-in, at each vintage\n', ...
    nvintage, nsim, burnin);

%% ------------------------------------------------------------------
%  3. The recursive exercise
%  ------------------------------------------------------------------
yhat1 = zeros(nvintage, 3*n+1);
yhat4 = zeros(0, 3*n+1);
kappa_path = zeros(nvintage, 2);
rng(20260908, 'twister')
fprintf('\nvintage   t    kappa1      kappa2\n');

t_run = tic;
for iv = 1:nvintage
    t = t_first + iv - 1;
    [store_alp, store_beta, store_h_T, store_Sigh, kappa_hat] = ...
        estimate_vintage(Y0, Y, t, p, n, nsim, burnin, ah, Vh, nuh0, Sh0, ...
                         c01, c02, lam0_nu_psi);
    kappa_path(iv,:) = kappa_hat';
    fprintf('%5d %5d %10.5f %11.7f\n', iv, t, kappa_hat(1), kappa_hat(2));

        % one forecast per kept draw, then one accumulation row per horizon
    cfg = struct('Yt', Y(1:t,:), 'Y', Y, 'p', p, 't', t, 'T', T);
    tmp1 = zeros(nsim, 2*n+1);  tmp4 = zeros(nsim, 2*n+1);
    for isim = 1:nsim
        draw = struct('alp', store_alp(isim,:)', 'beta', store_beta(isim,:)', ...
            'h_T', store_h_T(isim,:)', 'Sigh', store_Sigh(isim,:)');
        fcr = bvar.forecast.iterate('mahp_sv', draw, cfg);
        tmp1(isim,:) = fcr(1,:);
        tmp4(isim,:) = fcr(2,:);
    end
    yhat1(iv,:) = bvar.forecast.tables('accum_row', tmp1, Y(t+1,:));
    if t <= T-4                                   % the legacy storage guard
        yhat4(end+1,:) = bvar.forecast.tables('accum_row', tmp4, Y(t+4,:)); %#ok<SAGROW>
    end
end
elapsed = toc(t_run);
fprintf('\n%d vintages in %.0f s\n', nvintage, elapsed);

%% ------------------------------------------------------------------
%  4. What iterate returns, and what tables makes of it
%  ------------------------------------------------------------------
fprintf(['\niterate returns a 2 x (2n+1) matrix per draw: row 1 for h = 1 and\n' ...
         'row 2 for h = 4, each holding n point forecasts, n per-variable log\n' ...
         'predictive likelihoods, and the joint log predictive likelihood over\n' ...
         'all %d variables. tables(''accum_row'', ...) collapses the %d draws of\n' ...
         'one vintage into a single row: the posterior mean of the point\n' ...
         'forecasts, and the log of the MEAN predictive density rather than the\n' ...
         'mean of the logs, computed with a max-subtraction for stability.\n'], n, nsim);

S = bvar.forecast.tables('mahp', yhat1, yhat4);
vnames = compose("var%d", 1:n);

fprintf('\n%s\n', repmat('-', 1, 58));
fprintf('%-10s %12s %12s %10s %10s\n', 'variable', 'RMSFE h=1', 'RMSFE h=4', 'ALPL h=1', 'ALPL h=4');
fprintf('%s\n', repmat('-', 1, 58));
for i = 1:n
    fprintf('%-10s %12.4f %12.4f %10.3f %10.3f\n', vnames(i), ...
        S.RMSFE_1(i), S.RMSFE_4(i), S.aveprelike_1(i), S.aveprelike_4(i));
end
fprintf('%-10s %12s %12s %10.3f %10.3f\n', 'joint', '', '', ...
    S.aveprelike_1(end), S.aveprelike_4(end));
fprintf('%s\n', repmat('-', 1, 58));
fprintf(['tables trims the evaluation window the way the paper does: the h = 1\n' ...
         'rows start at the fourth vintage so that both horizons are scored over\n' ...
         'the same period. With %d vintages that leaves %d rows at h = 1 and %d at\n' ...
         'h = 4, which is why these numbers rank nothing.\n'], ...
         nvintage, size(yhat1,1)-3, size(yhat4,1));

%% ------------------------------------------------------------------
%  5. The prior adapting
%  ------------------------------------------------------------------
figure('Name','ex08: shrinkage learned at each vintage','Position',[100 100 760 320]);
subplot(1,2,1); plot(t_first:(t_first+nvintage-1), kappa_path(:,1), 'k-o', 'LineWidth', 1.1);
title('kappa_1  (own lags)'); xlabel('vintage t'); grid on
subplot(1,2,2); plot(t_first:(t_first+nvintage-1), kappa_path(:,2), 'k-o', 'LineWidth', 1.1);
title('kappa_2  (cross lags)'); xlabel('vintage t'); grid on

fprintf(['\nThe two panels show the posterior mean of each shrinkage\n' ...
         'hyperparameter at each vintage. A Minnesota prior fixes both in\n' ...
         'advance; here they are estimated, and they move as the sample grows.\n' ...
         'kappa2 sits about thirty times below kappa1 throughout, so cross-lag\n' ...
         'coefficients are shrunk much harder than own lags.\n']);

fprintf('\nex08 done. For the published exercise use\n');
fprintf('  replications/chan2021_ijf_mahp/legacy/main_forecasting.m\n');
fprintf('at its own settings: 23 variables, nsim = 20000, T0 = 91.\n');

%% ------------------------------------------------------------------
function [store_alp, store_beta, store_h_T, store_Sigh, kappa_hat] = ...
    estimate_vintage(Y0, Y, t, p, n, nsim, burnin, ah, Vh, nuh0, Sh0, ...
                     c01, c02, lam0_nu_psi)
% One Gibbs run on data through vintage t. Block order and constants follow
% chan2021_ijf_mahp/legacy/forecast_BVAR_MNG.m, including the two settings that
% differ from the estimation script: the psi floor is 1e-16 and the psi_kappa1
% initial scale is halved. Both are recorded in that package's preset.m.
Yt = Y(1:t, :);  Tt = size(Yt, 1);
[~, Zt] = bvar.util.build_lags([Y0(end-p+1:end,:); Yt], p);
k_beta = n^2*p + n;

sig2 = bvar.priors.resid_var_ar4(Y0, Yt);
[C, idx_kappa1, idx_kappa2] = bvar.priors.minnesota_C(n, p, sig2);

nu_psi = .5;
kappa = [.4, .001, 1, 100];
h0 = log(sig2);  h = repmat(h0', Tt, 1);  Sigh = Sh0;
psi_kappa1 = gamrnd(nu_psi, 2/nu_psi/2, n*p, 1);
psi_kappa2 = gamrnd(nu_psi, 2/nu_psi, (n-1)*n*p, 1);
Psi = ones(k_beta, 1);
Psi(idx_kappa1) = psi_kappa1;  Psi(idx_kappa2) = psi_kappa2;

store_kappa = zeros(nsim, 2);
store_beta  = zeros(nsim, k_beta);
store_alp   = zeros(nsim, n*(n-1)/2);
store_h_T   = zeros(nsim, n);
store_Sigh  = zeros(nsim, n);

for isim = 1:nsim + burnin
    [Valp, Vbeta] = bvar.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C.*Psi, sig2);
    Valp = 2*Valp;  Vbeta = 2*Vbeta;      % the forecast driver's doubling
    [beta, alp, U] = bvar.samplers.eq_gauss(Yt, Zt, h, Valp, Vbeta);

    for ij = 1:n
        h(:,ij) = bvar.sv.ksc_rw_h0(log(U(:,ij).^2 + .0001), h(:,ij), Sigh(ij), h0(ij));
    end

    [kappa, psi_kappa1, psi_kappa2] = bvar.samplers.gig_shrinkage('mng', ...
        beta, idx_kappa1, idx_kappa2, C, kappa, psi_kappa1, psi_kappa2, ...
        nu_psi, c01, c02, n, p, 1e-16);
    Psi(idx_kappa1) = psi_kappa1;  Psi(idx_kappa2) = psi_kappa2;

    nu_psi = bvar.samplers.nu_psi_ng(psi_kappa1, psi_kappa2, nu_psi, lam0_nu_psi);

    Kh0 = sparse(1:n, 1:n, 1./Sigh + 1./Vh);
    h0 = Kh0\(ah./Vh + h(1,:)'./Sigh) + chol(Kh0,'lower')'\randn(n,1);
    e = h - [h0'; h(1:Tt-1,:)];
    Sigh = 1./gamrnd(nuh0 + Tt/2, 1./(Sh0 + sum(e.^2)'/2));

    if isim > burnin
        isave = isim - burnin;
        store_kappa(isave,:) = kappa(1:2);
        store_beta(isave,:)  = beta';
        store_alp(isave,:)   = alp';
        store_h_T(isave,:)   = h(end,:)';
        store_Sigh(isave,:)  = Sigh';
    end
end
kappa_hat = mean(store_kappa)';
end
