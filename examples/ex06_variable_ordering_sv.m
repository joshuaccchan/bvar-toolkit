%% ex06 - Variable ordering in a VAR with stochastic volatility
%
% BOOK: Chapter 13, Section 13.1.3, and Chapter 14, Section 14.3, in Bayesian
% Macroeconometrics: Methods and Applications (Chapman & Hall/CRC, 2027).
%
% Two models of the reduced-form covariance, each estimated in the published
% variable order and in the reverse order:
%
%   Cholesky SV, as in ex04:  Sigma_t = A^{-1} diag(exp(h_t)) A^{-T},
%                             A unit lower triangular
%   Order-invariant SV:       Sigma_t = B0^{-1} diag(exp(h_t)) B0^{-T},
%                             B0 unrestricted, log-volatilities zero-mean
%
% Each model is also run a second time in the published order under another
% seed, which measures the Monte Carlo error the ordering gap is read against.
% Both models are estimated by bvar.models.var_sv.
%
% DATA. replications/chan_koop_yu2024_jbes_oisv/legacy/FRED_MD_20vars.csv,
% read-only: monthly FRED-MD, 1959:03 to 2019:12, columns 4, 6, 12 and 13 -
% industrial production, the unemployment rate, PCE inflation and the federal
% funds rate, which the paper ranks first among its 20 variables. IP and PCE
% inflation are multiplied by 100 to convert them to percent. Thirteen lags, the
% first 24 months as initial conditions, and the prior constants of the
% package's preset.m. The chains are 1000 draws after 200 burn-in, against the
% published 30,000 after 5,000.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.
% Cogley, T. and Sargent, T.J. (2005). Drifts and Volatilities: Monetary
% Policies and Outcomes in the Post WWII US, Review of Economic Dynamics,
% 8(2): 262-302.
% McCracken, M.W. and Ng, S. (2016). FRED-MD: A Monthly Database for
% Macroeconomic Research, Journal of Business and Economic Statistics,
% 34(4): 574-589.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

repo = fileparts(fileparts(mfilename('fullpath')));
fprintf('\n=== ex06: variable ordering in a VAR with stochastic volatility ===\n');

%% ------------------------------------------------------------------
%  1. Data, the two orderings, and the settings
%  ------------------------------------------------------------------
pkg = fullfile(repo, 'replications', 'chan_koop_yu2024_jbes_oisv');
od = cd(pkg); guard = onCleanup(@() cd(od));     % preset.m resolves in its own folder
pr = preset();
clear guard

data = load(fullfile(pkg, 'legacy', pr.data_file));   % never written to; legacy is frozen
data(:, [4 12]) = 100*data(:, [4 12]);                % IP and PCE inflation in percent
var_id = pr.var_id1;                                  % [4 6 12 13], the paper's first four
vnames = ["IP" "unemployment" "PCE inflation" "fed funds"];
p = pr.p;
n0 = pr.n0;

nsim = 1000; burnin = 200;
seeds = [20260915, 31415926];

n = numel(var_id);
fprintf('\nvariables, in the published order: %s\n', strjoin(cellstr(vnames), ', '));
fprintf('%d monthly observations from 1959:03, %d as initial conditions, p = %d lags\n', ...
    size(data,1), n0, p);
fprintf('%d draws after %d burn-in, two seeds\n', nsim, burnin);

%% ------------------------------------------------------------------
%  2. Six chains: two models x two orderings, plus a second seed for each
%     model in the published order
%  ------------------------------------------------------------------
cfg = struct('model', {'OI','OI','OI','CS','CS','CS'}, ...
             'rev',   {false, true, false, false, true, false}, ...
             'seed',  {seeds(1), seeds(1), seeds(2), seeds(1), seeds(1), seeds(2)}, ...
             'label', {'OI order 1','OI order 2','OI order 1, seed 2', ...
                       'CS order 1','CS order 2','CS order 1, seed 2'});

res = cell(1, numel(cfg));
for ic = 1:numel(cfg)
    ord = var_id;
    if cfg(ic).rev, ord = ord(end:-1:1); end
    Y0 = data(1:n0, ord);
    Y  = data(n0+1:end, ord);

    fprintf('\n%-20s ', cfg(ic).label);
    t0 = tic;
    res{ic} = bvar.models.var_sv(Y0, Y, p, 'model', cfg(ic).model, ...
        'nsim', nsim, 'burnin', burnin, 'seed', cfg(ic).seed);
    if cfg(ic).rev                       % back to the published order for comparison
        res{ic}.Sig_mean = res{ic}.Sig_mean(:, end:-1:1, end:-1:1);
    end
    fprintf('%5.0f s\n', toc(t0));
end

%% ------------------------------------------------------------------
%  3. Every element of Sigma_t, against the two-seed yardstick
%  ------------------------------------------------------------------
%  Variances are compared as a share of their level. Correlations are compared
%  in correlation points, because a correlation near zero can change sign while
%  moving very little.
vpath = @(S,i) S(:,i,i);
cpath = @(S,i,j) S(:,i,j)./sqrt(S(:,i,i).*S(:,j,j));
relgap = @(a,b) 100*mean(abs(a-b)) / mean(abs(a)+abs(b))*2;
absgap = @(a,b) mean(abs(a-b));
col = [4 5; 4 6; 1 2; 1 3];                    % CS orderings, CS seeds, OI orderings, OI seeds

fprintf('\ndifference between two posterior-mean paths\n');
fprintf('%-34s %11s %11s %11s %11s\n', 'variances, as a share of their level', ...
    'CS orders', 'CS seeds', 'OI orders', 'OI seeds');
for ii = 1:n
    fprintf('  %-32s', sprintf('var(%s)', vnames(ii)));
    for ic = 1:4
        fprintf(' %10.1f%%', relgap(vpath(res{col(ic,1)}.Sig_mean, ii), vpath(res{col(ic,2)}.Sig_mean, ii)));
    end
    fprintf('\n');
end

prs = nchoosek(1:n, 2);
gap = zeros(size(prs,1), 4);
fprintf('%-34s\n', 'correlations, in correlation points');
for ip = 1:size(prs,1)
    i = prs(ip,1); j = prs(ip,2);
    for ic = 1:4
        gap(ip,ic) = absgap(cpath(res{col(ic,1)}.Sig_mean, i, j), cpath(res{col(ic,2)}.Sig_mean, i, j));
    end
    fprintf('  %-32s %11.3f %11.3f %11.3f %11.3f\n', ...
        sprintf('corr(%s, %s)', vnames(i), vnames(j)), gap(ip,:));
end

[~, ip] = max(gap(:,1));                        % where the ordering moves CS most
i = prs(ip,1); j = prs(ip,2);
fprintf(['\nreversing the ordering moves the CS correlation between %s and %s by %.3f on\n' ...
    'average, against %.3f between two seeds of the same CS run. The same correlation\n' ...
    'under OI moves %.3f with the ordering and %.3f between seeds.\n'], ...
    vnames(i), vnames(j), gap(ip,1), gap(ip,2), gap(ip,3), gap(ip,4));

%% ------------------------------------------------------------------
%  4. The two panels of Figures 8 and 9, on four variables
%  ------------------------------------------------------------------
T = size(res{1}.Sig_mean, 1);
dates = datetime(1961,3,1) + calmonths(0:T-1);
lbl = {cfg([4 5 1 2]).label};
figure('Name', 'ex06: two orderings, two models');
subplot(2,1,1); hold on; box off
for ic = [4 5 1 2], plot(dates, vpath(res{ic}.Sig_mean, 4), 'LineWidth', 1); end
title('variance of the federal funds rate'); legend(lbl, 'Location', 'northeast'); hold off
subplot(2,1,2); hold on; box off
for ic = [4 5 1 2], plot(dates, cpath(res{ic}.Sig_mean, i, j), 'LineWidth', 1); end
title(sprintf('correlation between %s and %s', vnames(i), vnames(j))); hold off
drawnow

fprintf('\nThe order-invariant model is drawn under two orderings that are the same model,\n');
fprintf('so its two paths differ only by Monte Carlo error. The Cholesky model is drawn\n');
fprintf('under two orderings that are two different models.\n');
