%% build - regenerate the figures and numbers of tutorials/shrinkage/README.md
%
% The data and settings are those of the working-paper version of Chan (2022),
% archived in replications/chan2019wp_acp: 21 quarterly US variables, four lags,
% and the prior set on the structural-form coefficients. Part 1 evaluates the
% closed-form marginal likelihood under the three priors of the paper's Table 2
% and over the grid of its Figure 1, and checks the results against the capture in
% tests/golden. The figure it writes uses a logarithmic grid instead, which spans
% the three priors. Part 2 draws from the posterior under each prior and scans the
% lag length. Part 3 repeats the paper's recursive forecasting exercise, whose code
% is not in the package. The figures are written next to this file and everything
% printed goes to a log in tempdir.
%
% Usage, from anywhere:  run tutorials/shrinkage/build.m

tdir = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(tdir));
logf = fullfile(tempdir, 'bvar_shrinkage_build_log.txt');
if exist(logf, 'file'), delete(logf); end
diary(logf);
fprintf('tutorials/shrinkage/build.m, %s, MATLAB %s\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), version);
t0 = tic;
run(fullfile(repo, 'setup.m'));

    % the data as main_BVAR_ACP.m reads them
leg = fullfile(repo, 'replications', 'chan2019wp_acp', 'legacy');
raw = load(fullfile(leg, 'macrodata_Q_2018Q4.csv'));
var_id = [1,2,18,22,34,35,57,59,76,81,95,97,120,123,133,138,145,148,152,160,245];
data = raw(1:238, var_id);
vnames = ["GDP" "Consumption" "Disposable income" "Industrial production" ...
    "Capacity utilization" "Payroll employment" "Civilian employment" "Unemployment rate" ...
    "Hours" "Housing starts" "PCE prices" "GDP deflator" "CPI" "PPI" "Real earnings" ...
    "Productivity" "3-month T-bill" "10-year yield" "Baa spread" "Real M1" "S&P 500"];
qd = datetime(1959, 7, 1) + calquarters(0:size(data,1)-1);     % the first row is 1959Q3
[~, iu] = max(data(:,8));
assert(qd(iu) == datetime(1982, 10, 1), 'the unemployment rate should peak in 1982Q4');
p = 4;  n0 = 8;
Y0 = data(1:n0, :);  Y = data(n0+1:end, :);
[T, n] = size(Y);  k = n*p + 1;
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end,:); Y], p);
sig2 = bvar.priors.resid_var_ar4(Y0, Y);
k_subj = [.04, .0016, 1, 100];
labs = ["symmetric" "subjective" "asymmetric"];
fprintf('data: %s to %s, T = %d quarters after %d initial conditions, n = %d, p = %d\n', ...
    qlab(qd(n0+1)), qlab(qd(end)), T, n0, n, p);

%% ------------------------------------------------------------------
%  Part 1. The paper's Table 2 and Figure 1
%  ------------------------------------------------------------------
fprintf('\n=== Part 1: the marginal likelihood of the 21-variable VAR ===\n');
[ml_asym, k_asym] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [.04 .04], 'stru');
[ml_sym, k_sym] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [], 'stru', [], 'symmetric', true);
ml_subj = bvar.ml.acp(p, Y, Z, bvar.priors.acp_stru(n, p, k_subj, sig2));

    % the capture in tests/golden/chan2019wp_acp prints kappa to three decimals and
    % the log marginal likelihood to integers; Table 2 of the paper adds the subjective prior
assert(round(k_sym(1), 3) == 0.039 && round(k_asym(1), 3) == 0.406 && round(k_asym(2), 3) == 0.009 ...
    && round(ml_sym) == -9436 && round(ml_asym) == -9201, 'the optima differ from the golden capture');
assert(round(ml_subj) == -9372, 'the subjective prior does not reproduce Table 2 of the paper');
fprintf('optima equal the golden capture; Table 2 of the paper reproduced\n');

fprintf('\nTable 1\n');
fprintf('| Prior | kappa2 (own lags) | kappa3 (other lags) | log marginal likelihood | difference |\n');
fprintf('|---|---|---|---|---|\n');
K = {k_sym, k_subj, k_asym};  ML = [ml_sym, ml_subj, ml_asym];
for j = 1:3
    fprintf('| %s | %.4f | %.4f | %.1f | %.1f |\n', labs(j), K{j}(1), K{j}(2), ML(j), ML(j) - ml_asym);
end
fprintf('full precision: asymmetric (%.10g, %.10g) %.6f; symmetric %.10g %.6f; subjective %.6f\n', ...
    k_asym(1), k_asym(2), ml_asym, k_sym(1), ml_sym, ml_subj);
fprintf('kappa2 / kappa3: %.1f; kappa2 / symmetric kappa: %.1f; kappa3 / symmetric kappa: %.2f\n', ...
    k_asym(1)/k_asym(2), k_asym(1)/k_sym(1), k_asym(2)/k_sym(1));

    % the grid of the paper's Figure 1 (main_BVAR_ACP.m line 58); under a flat
    % prior on (kappa2, kappa3) the surface is proportional to their joint
    % posterior density, and the summaries below are grid approximations
[K2, K3] = meshgrid(0.25:.01:.65, .002:.0002:.02);
lml = zeros(size(K2));
for i = 1:numel(K2)
    lml(i) = bvar.ml.acp(p, Y, Z, bvar.priors.acp_stru(n, p, [K2(i), K3(i), 1, 100], sig2));
end
mlg = exp(lml - max(lml(:)));
[mx, im] = max(lml(:));
fprintf('grid of %d x %d over kappa2 in [%.2f, %.2f] and kappa3 in [%.3f, %.3f]: maximum %.3f at (%.2f, %.4f)\n', ...
    size(K2,1), size(K2,2), min(K2(:)), max(K2(:)), min(K3(:)), max(K3(:)), mx, K2(im), K3(im));
w = mlg / sum(mlg(:));
g2 = K2(1,:);  g3 = K3(:,1)';  m2 = sum(w, 1);  m3 = sum(w, 2)';
q = @(g, m, a) g(find(cumsum(m) >= a, 1));
fprintf('grid approximations on that support: kappa2 mean %.3f, 90%% interval [%.2f, %.2f]; kappa3 mean %.4f, 90%% interval [%.4f, %.4f]\n', ...
    g2*m2', q(g2, m2, .05), q(g2, m2, .95), g3*m3', q(g3, m3, .05), q(g3, m3, .95));
fprintf('mass at the edges of the grid: kappa2 %.2g, kappa3 %.2g\n', m2(1) + m2(end), m3(1) + m3(end));

    % the figure spans all three priors, which a linear grid around the optimum
    % leaves out, so it uses a logarithmic one; it shows the log marginal likelihood
    % minus its maximum in bands, so that the far priors sit on the surface
f2 = logspace(log10(0.02), log10(2), 101);
f3 = logspace(log10(0.0005), log10(0.2), 101);
[F2, F3] = meshgrid(f2, f3);
lf = zeros(size(F2));
for i = 1:numel(F2)
    lf(i) = bvar.ml.acp(p, Y, Z, bvar.priors.acp_stru(n, p, [F2(i), F3(i), 1, 100], sig2));
end
d = lf - ml_asym;                              % log marginal likelihood minus its maximum
lev = [-500 -200 -50 -10 -2];                  % the band edges
assert(min(d(:)) < lev(1) && max(d(:)) <= 1e-6, 'the grid must span every band, below the maximum');
nb = numel(lev) + 1;
edges = [min(d(:)) - 1, lev, 1];               % the bands, lowest first
    % a step colormap, one color per band: each entry takes the band of its upper end
u = edges(1) + (1:4096)*(edges(end) - edges(1))/4096;
cols = interp1([0; 1], [0.97 0.97 0.98; 0.16 0.33 0.58], linspace(0, 1, nb)');
fig = figure('Color', 'w', 'Position', [100 100 760 560]);
ax = axes(fig);
hold(ax, 'on');
contourf(ax, F2, F3, d, edges(1:end-1), 'LineStyle', 'none');
colormap(ax, cols(discretize(u, edges), :));
clim(ax, edges([1 end]));
[C, hc] = contour(ax, F2, F3, d, lev, 'LineColor', [0.35 0.35 0.35], 'LineWidth', 0.6);
clabel(C, hc, lev(lev < -2), 'FontSize', 10, 'Color', [0.15 0.15 0.15], 'LabelSpacing', 500);
dl = [max(f2(1), f3(1)), min(f2(end), f3(end))];
hm = [plot(ax, k_asym(1), k_asym(2), 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'w', ...
          'MarkerEdgeColor', 'k', 'LineWidth', 1), ...
      plot(ax, k_sym(1), k_sym(2), 'o', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w'), ...
      plot(ax, k_subj(1), k_subj(2), 's', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w'), ...
      plot(ax, dl, dl, 'k--', 'LineWidth', 1)];
xt = [0.05 0.1 0.2 0.5 1 2];  yt = [0.001 0.003 0.01 0.03 0.1];
set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12, 'TickDir', 'out', 'Layer', 'top', ...
    'XTick', xt, 'XTickLabel', compose('%g', xt), 'YTick', yt, 'YTickLabel', compose('%g', yt), ...
    'XMinorTick', 'off', 'YMinorTick', 'off');
xlim(ax, f2([1 end]));  ylim(ax, f3([1 end]));
xlabel(ax, '\kappa_2, own lags (smaller is tighter)', 'FontSize', 13);
ylabel(ax, '\kappa_3, other lags (smaller is tighter)', 'FontSize', 13);
box(ax, 'off');
legend(hm, {'asymmetric prior (maximum)', 'best symmetric prior', 'subjective prior', '\kappa_2 = \kappa_3'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
exportgraphics(fig, fullfile(tdir, 'fig_contour.png'), 'Resolution', 200);
close(fig)
fprintf('Figure 1: the symmetric prior is %.1f and the subjective prior %.1f below the maximum\n', ...
    ml_sym - ml_asym, ml_subj - ml_asym);

    % the same comparison in the 15-variable application of the published version
pkg22 = fullfile(repo, 'replications', 'chan2022_qe_acp');
pr22 = run_in(pkg22, 'preset');
d22 = xlsread(fullfile(pkg22, 'legacy', pr22.data_file)); %#ok<XLSRD>
Y0q = d22(1:pr22.n_init, pr22.d2.var_id);  Yq = d22(pr22.n_init+1:end, pr22.d2.var_id);
[~, Zq] = bvar.util.build_lags([Y0q(end-pr22.p+1:end,:); Yq], pr22.p);
lq_a = bvar.priors.acp_opt_kappa(Y0q, Yq, Zq, pr22.p, pr22.d2.kappa_init, 'redu', pr22.d2.idx_ns);
lq_s = bvar.priors.acp_opt_kappa(Y0q, Yq, Zq, pr22.p, [], 'redu', pr22.d2.idx_ns, 'symmetric', true);
fprintf('published version, 15 variables: log ML asymmetric %.1f, symmetric %.1f, difference %.1f\n', ...
    lq_a, lq_s, lq_a - lq_s);

%% ------------------------------------------------------------------
%  Part 2. What the choice changes: coefficients and lag length
%  ------------------------------------------------------------------
fprintf('\n=== Part 2: posterior means of the VAR coefficients, and the lag length ===\n');
nsim = 10000;
Abar = cell(1, 3);
for j = 1:3
    rng(20260919, 'twister');
    [Alp, Beta, Sg] = bvar.samplers.acp_theta_sig(Y0, Y, p, bvar.priors.acp_stru(n, p, K{j}, sig2), nsim);
    Bt = bvar.structural.reduced_form(Alp, Beta, Sg);
    Abar{j} = reshape(mean(Bt, 1), k, n);
end
own = false(n*p, n);
for l = 1:p, own((l-1)*n + (1:n), :) = logical(eye(n)); end
oth = ~own;
    % coefficient of y_{r,t-l} in equation i, times s_r/s_i; all prior means are zero
s = sqrt(sig2(:));
scl = repmat(s, p, 1) ./ s';
fprintf('%d posterior draws under each prior; mean absolute posterior mean of the lag\n', nsim);
fprintf('coefficients in standard-deviation units\n');
fprintf('| Prior | own lags | other lags | other lags, lag 1 | other lags, lags 2-%d |\n', p);
fprintf('|---|---|---|---|---|\n');
first = false(n*p, n);  first(1:n, :) = true;
Ld = cell(1, 3);
for j = 1:3
    Ld{j} = Abar{j}(2:end, :).*scl;
    fprintf('| %s | %.4f | %.4f | %.4f | %.4f |\n', labs(j), mean(abs(Ld{j}(own))), mean(abs(Ld{j}(oth))), ...
        mean(abs(Ld{j}(oth & first))), mean(abs(Ld{j}(oth & ~first))));
end
fprintf('own lags: asymmetric / symmetric %.2f, subjective / symmetric %.2f\n', ...
    mean(abs(Ld{3}(own)))/mean(abs(Ld{1}(own))), mean(abs(Ld{2}(own)))/mean(abs(Ld{1}(own))));
fprintf('other lags: asymmetric / symmetric %.2f, subjective / symmetric %.2f\n', ...
    mean(abs(Ld{3}(oth)))/mean(abs(Ld{1}(oth))), mean(abs(Ld{2}(oth)))/mean(abs(Ld{1}(oth))));

pmax = n0;
lag = zeros(pmax, 5);
for pp = 1:pmax
    [~, Zp] = bvar.util.build_lags([Y0(end-pp+1:end,:); Y], pp);
    [la, ka] = bvar.priors.acp_opt_kappa(Y0, Y, Zp, pp, [.04 .04], 'stru');
    [ls, ks] = bvar.priors.acp_opt_kappa(Y0, Y, Zp, pp, [], 'stru', [], 'symmetric', true);
    lag(pp,:) = [ka(1), ka(2), la, ks(1), ls];
end
assert(lag(p,3) == ml_asym, 'the lag scan at p = %d does not reproduce Part 1', p);
fprintf('\nTable 2\n');
fprintf('| p | kappa2 | kappa3 | log ML, asymmetric | kappa, symmetric | log ML, symmetric |\n');
fprintf('|---|---|---|---|---|---|\n');
for pp = 1:pmax
    fprintf('| %d | %.3f | %.4f | %.1f | %.3f | %.1f |\n', pp, lag(pp,:));
end
[~, pa] = max(lag(:,3));  [~, ps] = max(lag(:,5));
fprintf('lag length with the highest log marginal likelihood: asymmetric p = %d, symmetric p = %d\n', pa, ps);
fprintf('asymmetric minus symmetric, by lag length: %s\n', sprintf('%.1f ', lag(:,3) - lag(:,5)));

%% ------------------------------------------------------------------
%  Part 3. Recursive forecasts, 1985Q1 to the end of the sample
%  ------------------------------------------------------------------
fprintf('\n=== Part 3: recursive forecasts ===\n');
fprintf(['pseudo-out-of-sample: the exercise truncates this one vintage at each origin, so the\n' ...
    'early samples hold the revised values, and each forecast conditions on the hyperparameters\n' ...
    'that maximize the marginal likelihood at that origin\n']);
hs = [1 4];  H = max(hs);
origins = find(qd == datetime(1984, 10, 1)):size(data,1)-1;    % 1984Q4 to 2018Q3
nor = numel(origins);  nsim_f = 10000;
point = nan(nor, n, 2, 3);  lpl = nan(nor, n, 2, 3);  actual = nan(nor, n, 2);
jlpl = nan(nor, 3);  kpath = nan(nor, 3);
fprintf('%d forecast origins, %s to %s; %d posterior draws per prior and origin; p = %d throughout\n', ...
    nor, qlab(qd(origins(1))), qlab(qd(origins(end))), nsim_f, p);
for io = 1:nor
    t = origins(io);
    Yf = data(n0+1:t, :);  Yf1 = data(n0+1:t+1, :);
    [~, Zf] = bvar.util.build_lags([Y0(end-p+1:end,:); Yf], p);
    [~, Zf1] = bvar.util.build_lags([Y0(end-p+1:end,:); Yf1], p);
    s2f = bvar.priors.resid_var_ar4(Y0, Yf);
    [~, ka] = bvar.priors.acp_opt_kappa(Y0, Yf, Zf, p, [.04 .04], 'stru');
    [~, ks] = bvar.priors.acp_opt_kappa(Y0, Yf, Zf, p, [], 'stru', [], 'symmetric', true);
    kpath(io,:) = [ka(1), ka(2), ks(1)];
    Kf = {ks, k_subj, ka};
    ylag = data(t:-1:t-p+1, :)';
    for ih = 1:2
        if t + hs(ih) <= size(data,1), actual(io,:,ih) = data(t+hs(ih), :); end
    end
    for j = 1:3
        prior = bvar.priors.acp_stru(n, p, Kf{j}, s2f);
            % all 21 variables one quarter ahead, exactly: with the prior known at t
            % held fixed, p(y_t+1 | y_1:t) = p(y_1:t+1) / p(y_1:t)
        jlpl(io,j) = bvar.ml.acp(p, Yf1, Zf1, prior) - bvar.ml.acp(p, Yf, Zf, prior);
        rng(20260919 + t, 'twister');       % common random numbers across the priors
        [Alp, Beta, Sg] = bvar.samplers.acp_theta_sig(Y0, Yf, p, prior, nsim_f);
        [Bt, St] = bvar.structural.reduced_form(Alp, Beta, Sg);
        [mu, sd] = bvar.forecast.predictive(Bt, St, ylag, H);
        for ih = 1:2
            if isnan(actual(io,1,ih)), continue; end
            m = mu(:,:,hs(ih));  sv = sd(:,:,hs(ih));
            point(io,:,ih,j) = mean(m, 1);
            lpl(io,:,ih,j) = logmeanexp(-0.5*log(2*pi) - log(sv) - 0.5*((actual(io,:,ih) - m)./sv).^2);
        end
    end
end
fprintf('optimal kappa over the origins: kappa2 %.3f to %.3f, kappa3 %.4f to %.4f, symmetric kappa %.3f to %.3f\n', ...
    min(kpath(:,1)), max(kpath(:,1)), min(kpath(:,2)), max(kpath(:,2)), min(kpath(:,3)), max(kpath(:,3)));
fprintf('number of forecasts: h = 1: %d, h = 4: %d\n', nnz(~isnan(actual(:,1,1))), nnz(~isnan(actual(:,1,2))));

    % gains of the asymmetric prior as in the paper: 100(1 - RMSFE ratio) and
    % 100 x (ALPL difference); DM positive when the asymmetric prior is more accurate
gR = zeros(n, 2, 2);  gA = zeros(n, 2, 2);  zR = zeros(n, 2, 2);  zA = zeros(n, 2, 2);
for v = 1:n
    for ih = 1:2
        ok = ~isnan(actual(:,v,ih));
        e = squeeze(point(ok,v,ih,:)) - actual(ok,v,ih);
        ls = squeeze(lpl(ok,v,ih,:));
        r = sqrt(mean(e.^2, 1));  a = mean(ls, 1);
        for b = 1:2                          % benchmark: 1 symmetric, 2 subjective
            gR(v,ih,b) = 100*(1 - r(3)/r(b));
            gA(v,ih,b) = 100*(a(3) - a(b));
            zR(v,ih,b) = dm_stat(e(:,b).^2 - e(:,3).^2, hs(ih));
            zA(v,ih,b) = dm_stat(ls(:,3) - ls(:,b), hs(ih));
        end
    end
end
bn = ["symmetric" "subjective"];
fprintf('\nTable 3 (gains of the asymmetric prior: RMSFE in percent, ALPL as 100 x the\n');
fprintf('difference in average log predictive likelihood; significant = two-sided DM test at 5%%)\n');
fprintf('| Benchmark | h | median RMSFE gain | RMSFE gains / significant / significant losses | median 100 x ALPL difference | ALPL gains / significant / significant losses |\n');
fprintf('|---|---|---|---|---|---|\n');
for b = 1:2
    for ih = 1:2
        fprintf('| %s | %d | %.2f | %d / %d / %d | %.2f | %d / %d / %d |\n', bn(b), hs(ih), ...
            median(gR(:,ih,b)), nnz(gR(:,ih,b) > 0), nnz(zR(:,ih,b) > 1.96), nnz(zR(:,ih,b) < -1.96), ...
            median(gA(:,ih,b)), nnz(gA(:,ih,b) > 0), nnz(zA(:,ih,b) > 1.96), nnz(zA(:,ih,b) < -1.96));
    end
end
fprintf('\nTable 4 (gains of the asymmetric prior by variable: RMSFE in percent, ALPL as 100 x the\n');
fprintf('difference; stars: two-sided DM test at 10/5/1%%)\n');
fprintf('| Variable | RMSFE h=1, vs symmetric | RMSFE h=4, vs symmetric | ALPL h=1, vs symmetric | ALPL h=4, vs symmetric | RMSFE h=1, vs subjective | RMSFE h=4, vs subjective | ALPL h=1, vs subjective | ALPL h=4, vs subjective |\n');
fprintf('|---|---|---|---|---|---|---|---|---|\n');
for v = 1:n
    fprintf('| %s |', vnames(v));
    for b = 1:2
        for mtr = 1:2
            for ih = 1:2
                if mtr == 1, g = gR(v,ih,b); z = zR(v,ih,b); else, g = gA(v,ih,b); z = zA(v,ih,b); end
                fprintf(' %.1f%s |', g, stars(z));
            end
        end
    end
    fprintf('\n');
end
fprintf('\njoint log predictive likelihood of all %d variables one quarter ahead, exact (a difference\n', n);
fprintf('of log marginal likelihoods), summed over %d origins: symmetric %.1f, subjective %.1f, asymmetric %.1f\n', ...
    nor, sum(jlpl, 1));
fprintf('DM, asymmetric against: symmetric %.2f, subjective %.2f\n', ...
    dm_stat(jlpl(:,3) - jlpl(:,1), 1), dm_stat(jlpl(:,3) - jlpl(:,2), 1));

fig = figure('Color', 'w', 'Position', [100 100 1100 640]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ttl = {'RMSFE gain against the symmetric prior', 'RMSFE gain against the subjective prior'; ...
       'ALPL gain against the symmetric prior', 'ALPL gain against the subjective prior'};
for mtr = 1:2
    for b = 1:2
        nexttile; hold on; box off
        if mtr == 1, G = gR(:,:,b); else, G = gA(:,:,b); end
        bh = bar(G, 'grouped', 'BarWidth', 0.9);
        bh(1).FaceColor = [0 0.447 0.741];  bh(2).FaceColor = [0.85 0.325 0.098];
        yline(0, 'k-');
        set(gca, 'XTick', 1:n, 'XTickLabel', vnames, 'XTickLabelRotation', 60, 'FontSize', 8);
        xlim([0.4 n+0.6]);
        title(ttl{mtr,b}, 'FontSize', 10);
        if mtr == 1 && b == 1, legend({'one quarter ahead', 'four quarters ahead'}, 'Location', 'northwest', 'Box', 'off'); end
    end
end
exportgraphics(fig, fullfile(tdir, 'fig_forecasts.png'), 'Resolution', 150);
close(fig)

fprintf('\nbuild finished in %.1f minutes\n', toc(t0)/60);
diary off

function out = run_in(folder, fname)
% call a function that exists under the same name in several packages, from its folder
od = cd(folder);  back = onCleanup(@() cd(od));
out = feval(fname);
end

function s = stars(z)
% significance stars for a two-sided test at the 10, 5 and 1% levels
s = repmat('*', 1, (abs(z) > 1.645) + (abs(z) > 1.96) + (abs(z) > 2.576));
end

function v = logmeanexp(x)
% log(mean(exp(x))) down each column, without overflow
m = max(x, [], 1);
v = m + log(mean(exp(x - m), 1));
end

function z = dm_stat(d, h)
% Diebold-Mariano statistic for the loss differential d, with a Newey-West
% long-run variance using h-1 lags
T = numel(d);  u = d - mean(d);
lrv = (u'*u)/T;
for l = 1:h-1, lrv = lrv + 2*(1 - l/h)*(u(1+l:end)'*u(1:end-l))/T; end
z = mean(d)/sqrt(lrv/T);
end

function s = qlab(d)
s = sprintf('%dQ%d', year(d), quarter(d));
end
