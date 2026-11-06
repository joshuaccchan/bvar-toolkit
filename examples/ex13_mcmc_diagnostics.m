%% ex13 - Diagnostics for MCMC output: how well do the VAR-SV samplers mix?
%
% BOOK: Chapter 6, Section 6.5, in Bayesian Macroeconometrics: Methods and
% Applications (Chapman & Hall/CRC, 2027).
%
% The two samplers of ex06, the Cholesky and the order-invariant VAR with
% stochastic volatility, are run by bvar.models.var_sv with their parameter
% draws kept, and bvar.diag measures how well each chain mixes. All three
% diagnostics use the long-run variance of a chain, its spectral density at
% frequency zero:
%
%   inefficiency factor   IF = 1 + 2*sum_l rho(l), the number of draws that
%                         carry the information of one independent draw; R/IF
%                         is the effective sample size
%   Monte Carlo standard  sqrt(IF*var/R), the simulation error of a posterior
%   error                 mean
%   Geweke's Z            the mean of the first 10% of the chain against the mean
%                         of its second half, over their standard error;
%                         approximately N(0,1) if the chain is stationary
%
% DATA. As in ex06: replications/chan_koop_yu2024_jbes_oisv/legacy/
% FRED_MD_20vars.csv, read-only, monthly 1959:03 to 2019:12, industrial
% production, the unemployment rate, PCE inflation and the federal funds rate, 13
% lags, the first 24 months as initial conditions. IP and PCE inflation are
% multiplied by 100 to convert them to percent. Each chain keeps 3000 draws
% after 500 burn-in.
%
% See:
% Geweke, J. (1992). Evaluating the Accuracy of Sampling-Based Approaches to the
% Calculation of Posterior Moments. In: J.M. Bernardo, J.O. Berger, A.P. Dawid and
% A.F.M. Smith (Eds), Bayesian Statistics 4, 169-193, Oxford University Press.
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

repo = fileparts(fileparts(mfilename('fullpath')));
fprintf('\n=== ex13: diagnostics for the output of the VAR-SV samplers ===\n');

%% ------------------------------------------------------------------
%  1. The two samplers, with their parameter draws kept
%  ------------------------------------------------------------------
pkg = fullfile(repo, 'replications', 'chan_koop_yu2024_jbes_oisv');
od = cd(pkg); guard = onCleanup(@() cd(od));     % preset.m resolves in its own folder
pr = preset();
clear guard

data = load(fullfile(pkg, 'legacy', pr.data_file));   % never written to; legacy is frozen
data(:, [4 12]) = 100*data(:, [4 12]);                % IP and PCE inflation in percent, as in ex06
var_id = pr.var_id1;                                  % [4 6 12 13], as in ex06
vnames = ["IP" "unemployment" "PCE inflation" "fed funds"];
p = pr.p; n0 = pr.n0; n = numel(var_id);
Y0 = data(1:n0, var_id);
Y  = data(n0+1:end, var_id);
k = 1 + n*p;

nsim = 3000; burnin = 500; seed = 20260919;
models = {'CS', 'OI'};
mlabel = {'Cholesky SV', 'order-invariant SV'};
res = cell(1, 2);
fprintf('\n%d draws after %d burn-in\n', nsim, burnin);
for im = 1:2
    t0 = tic;
    res{im} = bvar.models.var_sv(Y0, Y, p, 'model', models{im}, 'nsim', nsim, ...
        'burnin', burnin, 'seed', seed, 'draws', true);
    fprintf('%-20s %4.0f s\n', mlabel{im}, toc(t0));
end

%% ------------------------------------------------------------------
%  2. How far the autocorrelations reach, and the truncation lag L
%  ------------------------------------------------------------------
%  The long-run variance sums the autocovariances up to lag L, weighted by
%  1 - l/(L+1). Once L passes the lag at which the autocorrelations die out, the
%  inefficiency factor stops changing; before that it understates the
%  dependence.
D = res{2}.draws;                                     % the order-invariant chain
own = @(i) (i-1)*k + 1 + i;                           % own first lag of equation i in A(:)
IFd = bvar.diag.inefficiency_factor(D.impact(:, 1:n+1:end), 200);
[~, is] = max(IFd);                                   % the slowest diagonal element of B0
three = [D.A(:, own(1)), D.kappa(:, 2), D.impact(:, (is-1)*n + is)];
tlabel = {'own first lag, IP equation', 'cross-lag shrinkage kappa_2', ...
    sprintf('B0(%d,%d), %s equation', is, is, vnames(is))};
Ls = [20 50 100 200];
fprintf('\ninefficiency factor by truncation lag L, order-invariant model\n');
fprintf('%-32s', ''); for L = Ls, fprintf('%8s', sprintf('L = %d', L)); end; fprintf('\n');
for j = 1:3
    fprintf('  %-30s', tlabel{j});
    for L = Ls, fprintf('%8.1f', bvar.diag.inefficiency_factor(three(:,j), L)); end
    fprintf('\n');
end

lags = 0:300;
acf = @(x) arrayfun(@(l) ((x(1+l:end) - mean(x))'*(x(1:end-l) - mean(x)))/sum((x - mean(x)).^2), lags);
figure('Name', 'ex13: autocorrelations of the draws');
hold on; box off
for j = 1:3, plot(lags, acf(three(:,j)), 'LineWidth', 1); end
yline(0, ':');
xlabel('lag'); title('autocorrelation of the draws, order-invariant model');
legend(tlabel, 'Location', 'northeast');
hold off
drawnow

%% ------------------------------------------------------------------
%  3. Mixing by block, and the Monte Carlo error of the posterior means
%  ------------------------------------------------------------------
L = 200;
fprintf('\ninefficiency factors at L = %d, median (largest) within each block\n', L);
fprintf('%-34s %20s %20s\n', '', mlabel{:});
blk = {'shrinkage kappa', 'kappa'; 'SV persistence phi', 'phi'; ...
       'SV innovation variance sig2', 'sig2'; 'log-volatility mean mu', 'mu'; ...
       'impact matrix', 'impact'; 'VAR coefficients', 'A'};
worst = cell(1, 2);
for b = 1:size(blk, 1)
    fprintf('  %-32s', blk{b,1});
    for im = 1:2
        Dm = res{im}.draws;
        if ~isfield(Dm, blk{b,2})
            fprintf('%20s', '-');
            continue
        end
        X = Dm.(blk{b,2});
        if strcmp(blk{b,2}, 'impact') && strcmp(models{im}, 'CS')
            X = X(:, tril(true(n), -1));             % the free elements below the diagonal
        end
        IF = bvar.diag.inefficiency_factor(X, L);
        fprintf('%20s', sprintf('%.1f (%.1f)', median(IF), max(IF)));
        rel = bvar.diag.mcse(X, L)./std(X);         % simulation error / posterior spread
        [r, jr] = max(rel);
        if isempty(worst{im}) || r > worst{im}{1}
            worst{im} = {r, param_name(blk{b,2}, jr, n, k, models{im}, vnames), nsim/IF(jr)};
        end
    end
    fprintf('\n');
end
for im = 1:2
    fprintf(['%s: the largest Monte Carlo standard error of a posterior mean is %.2f\n' ...
        '  posterior standard deviations (%s, about %.0f effective draws)\n'], ...
        mlabel{im}, worst{im}{:});
end

    % where the order-invariant chain is slow: the diagonal of B0 and phi
rho = corr(D.phi(:, is), D.impact(:, (is-1)*n + is));
fprintf(['\nOrder-invariant model: the diagonal of B0 has inefficiency factors %s, and\n' ...
    'across draws B0(%d,%d) and the persistence of the %s log-volatility have\n' ...
    'correlation %.2f. Scaling a row of B0 and shifting the level of that equation''s\n' ...
    'log-volatility leave Sigma_t unchanged. The prior on B0 and the zero mean of the\n' ...
    'log-volatility pin the scale, the latter less tightly the closer phi is to one.\n'], ...
    mat2str(round(IFd, 1)), is, is, vnames(is), rho);

%% ------------------------------------------------------------------
%  4. Geweke's Z for the scalar parameters
%  ------------------------------------------------------------------
%  The default truncation lag, floor(4*(n/100)^(2/9)) for a segment of n draws,
%  suits weakly dependent draws. A chain whose autocorrelations last hundreds of
%  lags needs a longer one.
[~, ~, info] = bvar.diag.geweke(res{1}.draws.kappa(:, 1));
fprintf(['\nGeweke''s Z with the default lags (%d for the first %d draws, %d for the\n' ...
    'last %d) and with L = %d\n'], info.LA, numel(info.A), info.LB, numel(info.B), L);
fprintf('%-34s %20s %20s\n', '', mlabel{:});
fprintf('%-34s %20s %20s\n', '', sprintf('default  L = %d', L), sprintf('default  L = %d', L));
both = {{}, {}};                                  % reject with both lags
for b = 1:3
    for i = 1:size(res{1}.draws.(blk{b,2}), 2)
        nm = param_name(blk{b,2}, i, n, k, 'OI', vnames);
        fprintf('  %-32s', nm);
        for im = 1:2
            x = res{im}.draws.(blk{b,2})(:, i);
            z = [bvar.diag.geweke(x), bvar.diag.geweke(x, [], [], L)];
            fprintf('%20s', sprintf('%7.2f %9.2f', z));
            if all(abs(z) > 1.96), both{im}{end+1} = nm; end
        end
        fprintf('\n');
    end
end
fprintf(['\n|Z| > 1.96 rejects stationarity at the 5%% level. A rejection with the default\n' ...
    'lag and none with the longer one points to slow mixing; a rejection with both\n' ...
    'means the chain has not settled, or is too short for that parameter.\n']);
for im = 1:2
    if isempty(both{im})
        fprintf('%s: no parameter rejects with both lags.\n', mlabel{im});
    else
        fprintf('%s: rejects with both lags for %s.\n', mlabel{im}, strjoin(both{im}, '; '));
    end
end

function s = param_name(field, j, n, k, model, vnames)
% the name of column j of a block of res.draws
switch field
    case 'kappa'
        s = sprintf('kappa_%d', j);
    case {'phi', 'sig2', 'mu'}
        s = sprintf('%s, %s', field, vnames(j));
    case 'impact'
        if strcmp(model, 'CS')                  % columns are the free elements
            free = find(tril(true(n), -1));
            [r, c] = ind2sub([n n], free(j));
            s = sprintf('impact(%d,%d)', r, c);
        else
            [r, c] = ind2sub([n n], j);
            s = sprintf('B0(%d,%d)', r, c);
        end
    case 'A'
        [r, c] = ind2sub([k n], j);
        s = sprintf('A(%d,%d)', r, c);
end
end

