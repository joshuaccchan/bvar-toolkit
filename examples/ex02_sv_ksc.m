%% ex02 - Univariate stochastic volatility by the KSC auxiliary mixture
%
% THE MODEL (random-walk stochastic volatility, known initial log-volatility):
%
%       y_t = exp(h_t/2) * eps_t,     eps_t ~ N(0,1)
%       h_t = h_{t-1} + u_t,          u_t   ~ N(0, sigh2)
%       h_1 ~ N(h0, sigh2)
%
% h_t is the log variance at time t. The model is nonlinear in the states, so
% the ex01 machinery does not apply directly. The Kim, Shephard and Chib
% (1998) trick makes it apply anyway:
%
%   1. square and log the data:  ystar_t = log(y_t^2) = h_t + log(eps_t^2).
%      This is LINEAR in h_t. The catch is that log(eps_t^2) is log-chi^2(1),
%      not normal.
%   2. approximate that log-chi^2(1) density by a 7-component mixture of
%      normals (the KSC constants, hard-coded in the sampler).
%   3. conditional on which mixture component applies at each t, the model IS
%      linear and Gaussian - so draw the whole h path with the ex01 precision
%      sampler, and draw the component indicators from a 7-point distribution.
%
% bvar.sv.ksc_rw_h0 does steps 2 and 3 in one call. Its whole body is 20 lines:
% a T x 7 table of mixture weights, then exactly the banded-precision draw of
% ex01. It consumes one rand(T,1) and one randn(T,1) per call.
%
% SIBLINGS IN core/+bvar/+sv/ (same idea, different state equation):
%   bvar.sv.ksc_rw_h0      random walk, KNOWN initial value h_1 ~ N(h0,sigh2)   <- here
%   bvar.sv.ksc_rw_diffuse random walk, diffuse start h_1 ~ N(0,Vh); returns S too
%   bvar.sv.ksc_ar1_mean   stationary AR(1) with a mean, h_t - mu = phi(h_{t-1}-mu) + u
%   bvar.sv.csv_armh       common (scalar) stochastic volatility for an n-variate
%                         VAR, drawn by an accept-reject Metropolis-Hastings step
%   bvar.sv.sv_params /
%   bvar.sv.sv0_params     draw (mu, phi, sigh2) for the AR(1) specifications
%
% WHAT TO LOOK AT: the posterior mean of h should track the simulated truth
% closely (correlation well above 0.9), the posterior of sigh2 should cover the
% true value, and the LAST section shows a trap that bites in practice - the
% offset inside log(y^2 + c) is not scale-free.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

rng(20260903, 'twister')
fprintf('\n=== ex02: univariate SV via the KSC 7-component mixture ===\n');

%% ------------------------------------------------------------------
%  1. Simulate an SV series
%  ------------------------------------------------------------------
T          = 400;
sigh2_true = 0.05;        % variance of the log-volatility innovations
h0_true    = log(1.0);    % initial log variance (so the series starts at sd 1)

h_true    = zeros(T,1);
h_true(1) = h0_true + sqrt(sigh2_true)*randn;
for t = 2:T
    h_true(t) = h_true(t-1) + sqrt(sigh2_true)*randn;
end
y = exp(h_true/2).*randn(T,1);

fprintf('\nsimulated T = %d, true sigh2 = %.3f, true h0 = %.3f\n', T, sigh2_true, h0_true);
fprintf('  volatility sd exp(h/2) ranges over [%.2f, %.2f]\n', ...
    min(exp(h_true/2)), max(exp(h_true/2)));

%% ------------------------------------------------------------------
%  2. Priors and chain initialization
%     sigh2 ~ IG(nuh0, Sh0)   with prior mean Sh0/(nuh0-1) = 0.01
%     h0    ~ N(ah0, Vh0)     a loose normal, as in the MAHP replication
%  ------------------------------------------------------------------
nuh0 = 5;  Sh0 = 0.01*(nuh0-1);       % preset values used throughout the toolkit
ah0  = 0;  Vh0 = 10;

sv_offset = 1e-4;    % the c in log(y^2 + c); guards log(0). See section 5.

nsim = 800; burnin = 200;
store_h     = zeros(nsim, T);
store_sigh2 = zeros(nsim, 1);
store_h0    = zeros(nsim, 1);

h     = zeros(T,1);   % chain init: flat log-volatility
sigh2 = 0.1;
h0    = 0;

ystar = log(y.^2 + sv_offset);        % the transformed data, computed ONCE

%% ------------------------------------------------------------------
%  3. The Gibbs sampler: three blocks, one line each
%  ------------------------------------------------------------------
tic
for isim = 1:nsim + burnin

        % (a) the whole log-volatility path, KSC mixture + precision sampler
    h = bvar.sv.ksc_rw_h0(ystar, h, sigh2, h0);

        % (b) sigh2 | h, h0        (conjugate inverse gamma)
    e     = h - [h0; h(1:T-1)];
    sigh2 = 1/gamrnd(nuh0 + T/2, 1/(Sh0 + sum(e.^2)/2));

        % (c) h0 | h, sigh2        (Gaussian; this is the MAHP h0 block for n = 1)
    Kh0    = 1/sigh2 + 1/Vh0;
    h0_hat = Kh0\(ah0/Vh0 + h(1)/sigh2);
    h0     = h0_hat + 1/sqrt(Kh0)*randn;

    if isim > burnin
        isave = isim - burnin;
        store_h(isave,:)     = h';
        store_sigh2(isave,:) = sigh2;
        store_h0(isave,:)    = h0;
    end
end
elapsed = toc;

%% ------------------------------------------------------------------
%  4. What came out
%  ------------------------------------------------------------------
h_mean = mean(store_h)';
h_ci   = quantile(store_h, [.05 .95])';

fprintf('\n%d sweeps (%d kept) in %.2f seconds (%.2f ms per sweep)\n', ...
    nsim+burnin, nsim, elapsed, 1e3*elapsed/(nsim+burnin));
fprintf('\nposterior summaries:\n');
fprintf('  %-8s %10s %10s %20s\n', 'param', 'true', 'post mean', '90% CI');
fprintf('  %-8s %10.3f %10.3f   [%7.3f, %7.3f]\n', 'sigh2', sigh2_true, ...
    mean(store_sigh2), quantile(store_sigh2,.05), quantile(store_sigh2,.95));
fprintf('  %-8s %10.3f %10.3f   [%7.3f, %7.3f]\n', 'h0', h0_true, ...
    mean(store_h0), quantile(store_h0,.05), quantile(store_h0,.95));

fprintf('\nrecovery of the log-volatility PATH:\n');
fprintf('  correlation(h_mean, h_true)       : %.3f\n', corr(h_mean, h_true));
fprintf('  RMSE(h_mean, h_true)              : %.3f\n', sqrt(mean((h_mean-h_true).^2)));
fprintf('  RMSE of the naive proxy log(y^2)  : %.3f  (a single observation is\n', ...
    sqrt(mean((log(y.^2+sv_offset)+1.2704-h_true).^2)));
fprintf('                                       a very noisy variance estimate)\n');
fprintf('  share of periods with h_true in the 90%% band: %.2f\n', ...
    mean(h_true >= h_ci(:,1) & h_true <= h_ci(:,2)));

%% ------------------------------------------------------------------
%  5. A practical trap: the offset in log(y^2 + c) is NOT scale-free
%
%     Multiplying y by a constant a shifts the true log-variance path by the
%     constant 2*log(a) and changes NOTHING else - the shape of h, its
%     innovation variance sigh2, everything else is invariant. That is true
%     of the model. It is NOT true of the sampler, because the offset c in
%     ystar = log(y^2 + c) has fixed units. Once y^2 falls to the order of c,
%     the transformed data stop responding to y at all and the estimated
%     path is dragged towards a constant.
%
%     Below: the SAME series divided by 1, 100, 1000 and 10000 - think of a
%     growth rate expressed in percent (1.0), as a raw log difference (0.01),
%     and as a monthly raw log difference (0.001). Each estimate is shifted
%     back by 2*log(a) so all four are directly comparable to the truth.
%
%     Rule of thumb used throughout this toolkit: scale your data so that
%     typical y^2 is orders of magnitude larger than c - macro series in
%     PERCENT, not in raw log differences.
%  ------------------------------------------------------------------
scales   = [1 100 1000 10000];
nsim_s   = 400; burnin_s = 100;
h_mean_s = zeros(T, numel(scales));
fprintf('\nscale sensitivity of the log(y^2 + %g) offset:\n', sv_offset);
fprintf('  %-14s %12s %14s %12s %12s\n', 'y divided by', 'median y^2', ...
    'corr with h', 'sd of h_mean', 'RMSE vs h');
for is = 1:numel(scales)
    a  = scales(is);
    ys = y/a;
    ystar_s = log(ys.^2 + sv_offset);
    hs = zeros(T,1); sigh2_s = 0.1; h0_s = 0; S = zeros(nsim_s, T);
    for isim = 1:nsim_s + burnin_s
        hs      = bvar.sv.ksc_rw_h0(ystar_s, hs, sigh2_s, h0_s);
        e       = hs - [h0_s; hs(1:T-1)];
        sigh2_s = 1/gamrnd(nuh0 + T/2, 1/(Sh0 + sum(e.^2)/2));
        Kh0     = 1/sigh2_s + 1/Vh0;
        h0_s    = Kh0\(ah0/Vh0 + hs(1)/sigh2_s) + 1/sqrt(Kh0)*randn;
        if isim > burnin_s, S(isim-burnin_s,:) = hs'; end
    end
        % shift back onto the original scale: dividing y by a moves the true
        % log variance by -2*log(a).
    h_mean_s(:,is) = mean(S)' + 2*log(a);
    fprintf('  %-14d %12.2e %14.3f %12.3f %12.3f\n', a, median(ys.^2), ...
        corr(h_mean_s(:,is), h_true), std(h_mean_s(:,is)), ...
        sqrt(mean((h_mean_s(:,is)-h_true).^2)));
end
fprintf('  true path, for reference:  sd %.3f\n', std(h_true));
fprintf('  -> as median y^2 approaches c = %g, the estimated path flattens\n', sv_offset);
fprintf('     (sd falls to %.0f%% of the well-scaled answer) and the RMSE blows up,\n', ...
    100*std(h_mean_s(:,end))/std(h_mean_s(:,1)));
fprintf('     even though the underlying series is literally the same data.\n');

%% ------------------------------------------------------------------
%  6. Figures (non-blocking)
%  ------------------------------------------------------------------
figure('Name','ex02 stochastic volatility');
subplot(3,1,1)
plot(1:T, y, 'Color', [.4 .4 .4]); box off
title('simulated SV series y_t')

subplot(3,1,2)
plot(1:T, h_true, 'k', 'LineWidth', 1.2); hold on
plot(1:T, h_mean, 'r', 'LineWidth', 1.2);
plot(1:T, h_ci, 'r:'); hold off; box off
title('log-volatility h_t: truth (black), posterior mean and 90% band (red)')

subplot(3,1,3)
plot(1:T, h_true, 'k', 'LineWidth', 1.2); hold on
plot(1:T, h_mean_s, 'LineWidth', 1); hold off; box off
title('same data at four scales, shifted back: the log(y^2 + c) offset flattens h')
legend([{'truth'}, arrayfun(@(a) sprintf('y/%d', a), scales, 'uni', 0)], ...
    'Location','best'); legend boxoff

drawnow

fprintf('\nex02 done. Next: ex03_minnesota_bvar (a BVAR with a Minnesota prior).\n');
