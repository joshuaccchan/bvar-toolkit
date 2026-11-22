%% your_data - choose the shrinkage hyperparameters of a BVAR on your data
%
% Set the file, the columns, their names, which of them are nonstationary and the
% lag lengths below, then run the script. It chooses the own-lag and other-lag
% shrinkage kappa2 and kappa3 of the asymmetric conjugate prior by maximizing the
% closed-form marginal likelihood, compares the result with the best symmetric
% prior (kappa2 = kappa3) and with the subjective values kappa2 = 0.04 and
% kappa3 = 0.0016, repeats the choice for each lag length up to pmax, plots the
% marginal likelihood around the optimum, forecasts H periods past the end of the
% sample under each prior, and writes everything to a csv and a mat file. The
% prior is set on the structural-form coefficients, as in Section 12.3 of
% Bayesian Macroeconometrics. The defaults use the 21-variable dataset of the
% tutorial and take a few seconds.
%
% Give cols as column numbers, which are read with readmatrix and count every
% column of the file, or as column names, which are read with readtable. The
% selected columns must have no missing values inside the sample; rows missing at
% either end are dropped. Setting datecol checks that the rows kept are evenly
% spaced. The series must be stationary, transformed as needed; that choice is
% separate from nonstationary below, which only centers the prior.
%
% The first n0 rows serve as initial conditions: at least max(pmax,4) of them.

repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(repo, 'setup.m'))
lastwarn('');       % the report records the last warning raised from here on

%% ---- settings ----
file    = fullfile(repo, 'replications', 'chan2019wp_acp', 'legacy', 'macrodata_Q_2018Q4.csv');
cols    = [1 2 18 22 34 35 57 59 76 81 95 97 120 123 133 138 145 148 152 160 245];
names   = ["GDP" "Consumption" "Disposable income" "Industrial production" ...
           "Capacity utilization" "Payroll employment" "Civilian employment" ...
           "Unemployment rate" "Hours" "Housing starts" "PCE prices" "GDP deflator" "CPI" ...
           "PPI" "Real earnings" "Productivity" "3-month T-bill" "10-year yield" ...
           "Baa spread" "Real M1" "S&P 500"];
datecol = [];       % the column holding the dates, by number or name; [] to skip the check
nonstationary = []; % positions in cols whose first own lag gets prior mean one
p       = 4;        % lags
pmax    = 8;        % longest lag length in the scan
n0      = 8;        % rows used as initial conditions
H       = 4;        % forecast horizons past the end of the sample; 0 for none
nsim_f  = 5000;     % posterior draws behind the forecasts
seed    = 1;
outdir  = tempdir;  % where the report goes; '' for none

%% ---- data ----
if isnumeric(cols)
    raw = readmatrix(file);                   % column numbers count every column of the file
    sel = raw(:, cols);
    dates = [];
    if ~isempty(datecol), dates = raw(:, datecol); end
else
    tbl = readtable(file, 'VariableNamingRule', 'preserve');
    sel = tbl{:, cols};
    dates = [];
    if ~isempty(datecol), dates = tbl{:, datecol}; end
end
assert(numel(names) == numel(cols), 'names and cols must have the same length');
assert(n0 >= max(pmax, 4), 'n0 must be at least max(pmax,4)');
keep = all(isfinite(sel), 2);
assert(any(keep), 'every row of the selected columns has a missing value');
lo = find(keep, 1);  hi = find(keep, 1, 'last');
assert(all(keep(lo:hi)), 'the selected columns have missing values inside the sample');
data = sel(lo:hi, :);
rows_kept = [lo hi];        % lo and hi are reused by the contour grid below
flat = std(data) == 0;
assert(~any(flat), 'these series are constant: %s', strjoin(names(flat), ', '));
if ~isempty(dates)
    d = dates(lo:hi);
    if isdatetime(d)
        step = days(diff(d));
        even = all(step > 0) && max(abs(step - median(step))) <= 3;   % months vary by up to 3 days
    else
        step = diff(double(d));
        even = all(step > 0) && max(abs(step - median(step))) <= 1e-6*max(1, abs(median(step)));
    end
    assert(even, 'the dates of the rows kept are not evenly spaced');
end
Y0 = data(1:n0, :);
Y  = data(n0+1:end, :);
[T, n] = size(Y);
fprintf('\n%d variables, %d observations after %d initial conditions\n', n, T, n0);

%% ---- the three priors at p lags ----
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end,:); Y], p);
sig2 = bvar.priors.resid_var_ar4(Y0, Y);
k_subj = [.04 .0016 1 100];
[ml_a, ka] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [.04 .04], 'stru', nonstationary);
[ml_s, ks] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [], 'stru', nonstationary, 'symmetric', true);
ml_0 = bvar.ml.acp(p, Y, Z, bvar.priors.acp_stru(n, p, k_subj, sig2, nonstationary));
    % the symmetric optimum is feasible under the asymmetric prior, so a search
    % started there should reach the same maximum
[ml_b, kb] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [ks(1) ks(1)], 'stru', nonstationary);
gap = ml_b - ml_a;
if gap > 0, ml_a = ml_b; ka = kb; end
fprintf('\nlog marginal likelihood at p = %d\n', p);
fprintf('%-42s %10s %10s %10s %10s\n', '', 'kappa2', 'kappa3', 'log ML', 'difference');
fprintf('%-42s %10.4f %10.4f %10.1f %10.1f\n', 'symmetric, chosen by marginal likelihood', ks(1), ks(2), ml_s, ml_s - ml_a);
fprintf('%-42s %10.4f %10.4f %10.1f %10.1f\n', 'subjective', .04, .0016, ml_0, ml_0 - ml_a);
fprintf('%-42s %10.4f %10.4f %10.1f %10.1f\n', 'asymmetric, chosen by marginal likelihood', ka(1), ka(2), ml_a, 0);
fprintf(['kappa2 and kappa3 are the own-lag and other-lag shrinkage of Section 12.3 of\n' ...
    'Bayesian Macroeconometrics; bvar.priors.acp_opt_kappa returns them as kappa(1:2).\n']);
fprintf('two starting values for the asymmetric search, (0.04, 0.04) and the symmetric\n');
fprintf('optimum (%.4f, %.4f), reach log marginal likelihoods %.4f apart\n', ks(1), ks(2), abs(gap));
fprintf('the estimates condition on these hyperparameters: their uncertainty is left out\n');

%% ---- the lag length ----
fprintf('\nlag length, with kappa2 and kappa3 chosen again for each\n');
fprintf('%4s %10s %10s %12s %14s\n', 'p', 'kappa2', 'kappa3', 'log ML', 'symmetric ML');
lml = zeros(pmax, 1);  lsym = zeros(pmax, 1);
kap = zeros(pmax, 2);  kaps = zeros(pmax, 2);
for pp = 1:pmax
    [~, Zp] = bvar.util.build_lags([Y0(end-pp+1:end,:); Y], pp);
    [lml(pp), kp] = bvar.priors.acp_opt_kappa(Y0, Y, Zp, pp, [.04 .04], 'stru', nonstationary);
    [lsym(pp), kq] = bvar.priors.acp_opt_kappa(Y0, Y, Zp, pp, [], 'stru', nonstationary, 'symmetric', true);
    kap(pp,:) = kp(1:2);  kaps(pp,:) = kq(1:2);
    fprintf('%4d %10.4f %10.5f %12.1f %14.1f\n', pp, kp(1), kp(2), lml(pp), lsym(pp));
end
[~, pbest] = max(lml);
fprintf('highest marginal likelihood at p = %d\n', pbest);
if pbest == pmax
    fprintf(['this is the longest lag length in the scan: raise pmax, and n0 with it,\n' ...
        'to see whether the marginal likelihood keeps rising\n']);
end

%% ---- forecasts from the end of the sample ----
% Given the coefficients and the error covariance matrix, the h-step forecast of a
% VAR is Gaussian, and bvar.forecast.predictive returns its mean and standard
% deviation for each posterior draw. Drawing one value from each of those normals
% gives the predictive distribution itself; its mean and standard deviation are
% printed, and the report also carries its 68% band under each of the three priors.
pname = ["symmetric"; "subjective"; "asymmetric"];
fcast = table();
if H > 0
    Kf = {ks, k_subj, ka};
    ylag = data(end:-1:end-p+1, :)';            % most recent lag first
    fm = zeros(n, H, 3);  fs = zeros(n, H, 3);  flo = zeros(n, H, 3);  fhi = zeros(n, H, 3);
    for j = 1:3
        prior = bvar.priors.acp_stru(n, p, Kf{j}, sig2, nonstationary);
        rng(seed, 'twister');                   % common random numbers across the priors
        [Alp, Beta, Sg] = bvar.samplers.acp_theta_sig(Y0, Y, p, prior, nsim_f);
        [Bt, St] = bvar.structural.reduced_form(Alp, Beta, Sg);
        [mu, sd] = bvar.forecast.predictive(Bt, St, ylag, H);
        yd = mu + sd.*randn(size(mu));
        q = quantile(yd, [.16 .84], 1);
        fm(:,:,j) = reshape(mean(mu, 1), n, H);
        fs(:,:,j) = reshape(sqrt(mean(sd.^2, 1) + var(mu, 0, 1)), n, H);
        flo(:,:,j) = reshape(q(1,:,:), n, H);
        fhi(:,:,j) = reshape(q(2,:,:), n, H);
    end
    fprintf('\nforecasts under the asymmetric prior, from %d posterior draws\n', nsim_f);
    fprintf('posterior mean, with the standard deviation of the predictive distribution\n');
    fprintf('%-24s', 'variable');
    for h = 1:H, fprintf('%18s', sprintf('h = %d', h)); end
    fprintf('\n');
    for i = 1:n
        fprintf('%-24.24s', names(i));
        for h = 1:H, fprintf('%10.2f (%5.2f)', fm(i,h,3), fs(i,h,3)); end
        fprintf('\n');
    end
    fprintf('the 68%% bands and the other two priors are in the report\n');
    fprintf('the draws condition on the chosen hyperparameters, so the spread leaves out\n');
    fprintf('the uncertainty about kappa2 and kappa3\n');
    vnames = names(:);
    [iv, ih, ip] = ndgrid(1:n, 1:H, 1:3);
    fcast = table(pname(ip(:)), vnames(iv(:)), ih(:), fm(:), fs(:), flo(:), fhi(:), ...
        'VariableNames', {'prior', 'variable', 'h', 'mean', 'sd', 'lo68', 'hi68'});
end

%% ---- the marginal likelihood around the optimum ----
    % a logarithmic grid from a tenth to ten times the optimum, widened to show
    % the symmetric and subjective priors
r2 = [min([ka(1)/10, ks(1), .04])/1.5, ka(1)*10];
r3 = [min([ka(2)/10, .0016])/1.5, max([ka(2)*10, ks(2)])*1.5];
g2 = logspace(log10(r2(1)), log10(r2(2)), 41);
g3 = logspace(log10(r3(1)), log10(r3(2)), 41);
[K2, K3] = meshgrid(g2, g3);
L = zeros(size(K2));
for i = 1:numel(K2)
    L(i) = bvar.ml.acp(p, Y, Z, bvar.priors.acp_stru(n, p, [K2(i), K3(i), 1, 100], sig2, nonstationary));
end
    % bands of the log marginal likelihood minus its maximum, at the edges the grid reaches
d = L - max(L(:));
lev = [-500 -200 -50 -10 -2];
lev = lev(lev > min(d(:)));
nb = numel(lev) + 1;
figure('Color', 'w'); hold on; box off
if nb > 1
    edges = [min(d(:)) - 1, lev, 1];                  % the bands, lowest first
    u = edges(1) + (1:4096)*(edges(end) - edges(1))/4096;   % a step colormap, one color per band
    cols = interp1([0; 1], [0.97 0.97 0.98; 0.16 0.33 0.58], linspace(0, 1, nb)');
    contourf(K2, K3, d, edges(1:end-1), 'LineStyle', 'none');
    colormap(gca, cols(discretize(u, edges), :));
    clim(edges([1 end]));
    cl = lev;
    if isscalar(cl), cl = [cl cl]; end                % one level, not a number of levels
    [C, hc] = contour(K2, K3, d, cl, 'LineColor', [0.35 0.35 0.35], 'LineWidth', 0.6);
    if any(lev < -2), clabel(C, hc, lev(lev < -2), 'FontSize', 9, 'LabelSpacing', 500); end
end
h = [plot(ka(1), ka(2), 'p', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'), ...
     plot(ks(1), ks(2), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w'), ...
     plot(.04, .0016, 's', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w')];
hl = {'asymmetric prior (maximum)', 'best symmetric prior', 'subjective prior'};
lo = max(g2(1), g3(1));  hi = min(g2(end), g3(end));
if lo < hi, h(end+1) = plot([lo hi], [lo hi], 'k--'); hl{end+1} = '\kappa_2 = \kappa_3'; end
set(gca, 'XScale', 'log', 'YScale', 'log', 'TickDir', 'out', 'Layer', 'top');
xlim(r2); ylim(r3);
xticklabels(compose('%g', xticks)); yticklabels(compose('%g', yticks));
legend(h, hl, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
xlabel('\kappa_2, own lags (smaller is tighter)'); ylabel('\kappa_3, other lags (smaller is tighter)');
title(sprintf('Log marginal likelihood minus its maximum, p = %d', p));

%% ---- the report ----
% One row per prior at p lags and per lag length in the scan, the forecasts in
% their own file, and the settings in the mat file.
if ~isempty(outdir)
    block = [repmat("prior comparison", 3, 1); repmat("lag length", 2*pmax, 1)];
    spec = [pname; repmat("asymmetric", pmax, 1); repmat("symmetric", pmax, 1)];
    pcol = [p; p; p; (1:pmax)'; (1:pmax)'];
    k2 = [ks(1); k_subj(1); ka(1); kap(:,1); kaps(:,1)];
    k3 = [ks(2); k_subj(2); ka(2); kap(:,2); kaps(:,2)];
    logml = [ml_s; ml_0; ml_a; lml; lsym];
    models = table(block, spec, pcol, k2, k3, logml, ...
        'VariableNames', {'block', 'prior', 'p', 'kappa2', 'kappa3', 'log_ML'});
    meta = struct('file', file, 'columns', string(cols(:)'), 'names', names, ...
        'nonstationary', nonstationary, 'p', p, 'pmax', pmax, 'n0', n0, 'n', n, 'T', T, ...
        'rows_kept', rows_kept, 'H', H, 'nsim_f', nsim_f, 'seed', seed, ...
        'kappa_asymmetric', ka, 'kappa_symmetric', ks, 'kappa_subjective', k_subj);
    out = struct('models', models);
    if H > 0, out.forecasts = fcast; end
    bvar.util.report('shrinkage_report', out, meta, outdir);
end
