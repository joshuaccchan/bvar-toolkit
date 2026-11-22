%% guide - the code of guides/var_forecast/README.md: forecasts of the next 20 quarters from a
% BVAR estimated on a sample from the true model of guides/var_simulation, one simulated path per
% posterior draw. Sections 2 to 4 are the code to adapt. Saves fig_fan.png and fig_widths.png
% next to this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
names = {'GDP growth', 'PCE inflation', 'unemployment'};

%% 1. the true model of guides/var_simulation, and the sample of 240 quarters simulated there
D = readmatrix(fullfile(fileparts(here), 'data', 'macro3_Q.csv'));
p = 4;
[Yd, Zd] = bvar.util.build_lags(D(:, 2:4), p);
Atrue = Zd\Yd;
Ud = Yd - Zd*Atrue;
Sigtrue = Ud'*Ud/(size(Yd, 1) - size(Zd, 2));
rng(2, 'twister');
data = bvar.util.simulate_var(Atrue, Sigtrue, 240);   % GDP growth, inflation, unemployment

%% 2. the posterior draws: the code to adapt
[Y, Z] = bvar.util.build_lags(data, p);           % the first four quarters are initial values
[T, n] = size(Y); k = size(Z, 2);
s2 = bvar.priors.resid_var_ar4(data(1:4, :), Y);  % AR(4) residual variances
nu0 = n + 2; S0 = diag(s2);                       % Sig ~ IW(nu0, S0)
c3 = 10^2;                                        % the intercepts
VA = @(c1) [c3; kron(1./(1:p)'.^2, c1./s2)];      % the diagonal of V_A
lc = fminbnd(@(lc) -bvar.ml.ncp(Y, Z, VA(exp(lc)), nu0, S0), log(1e-3), log(10));
c1 = exp(lc);                                     % the shrinkage chosen by the marginal likelihood
[~, Ahat, KA, Shat] = bvar.ml.ncp(Y, Z, VA(c1), nu0, S0);
nuhat = nu0 + T;
R = 10000;
LK = chol(KA, 'lower');                           % KA = LK*LK'
Adraw = zeros(k, n, R); Sdraw = zeros(n, n, R);
rng(3, 'twister');
for r = 1:R
    Sdraw(:, :, r) = iwishrnd(Shat, nuhat);       % Sig ~ IW(nuhat, Shat)
    Adraw(:, :, r) = Ahat + (LK'\randn(k, n))*chol(Sdraw(:, :, r), 'lower')';   % A given Sig
end

%% 3. one path per draw: the code to adapt
H = 20;
ylast = data(end-p+1:end, :);                     % the last four quarters, the oldest first
paths = zeros(H, n, R);
for r = 1:R
    paths(:, :, r) = bvar.util.simulate_var(Adraw(:, :, r), Sdraw(:, :, r), H, ylast);
end

%% 4. the forecasts: the code to adapt
q = quantile(paths, [.05 .16 .5 .84 .95], 3);     % H x n x 5: 90% and 68% bounds, median
fmean = mean(paths, 3);                           % the point forecast under squared error loss
pneg = mean(paths(1, 1, :) < 0);                  % P(GDP growth < 0) next quarter

%% 5. the same forecasts from the true parameters, and from the posterior means
Sigbar = Shat/(nuhat - n - 1);                    % the posterior mean of Sig
ptrue = zeros(H, n, R); pplug = zeros(H, n, R);
rng(4, 'twister');
for r = 1:R
    ptrue(:, :, r) = bvar.util.simulate_var(Atrue, Sigtrue, H, ylast);
    pplug(:, :, r) = bvar.util.simulate_var(Ahat, Sigbar, H, ylast);
end
qtrue = quantile(ptrue, [.05 .16 .5 .84 .95], 3);
qplug = quantile(pplug, [.05 .16 .5 .84 .95], 3);

%% the results
fprintf('\nVAR(%d) of %d variables on %d simulated quarters, T = %d\n', p, n, size(data, 1), T);
fprintf('last four quarters of the sample:\n'); disp(ylast);
fprintf(['c1 chosen by the marginal likelihood: %.4f (sqrt %.3f); log ML %.2f, ' ...
    'at c1 = 0.04 %.2f\n'], ...
    c1, sqrt(c1), bvar.ml.ncp(Y, Z, VA(c1), nu0, S0), bvar.ml.ncp(Y, Z, VA(.04), nu0, S0));
Ols = Z\Y;
fprintf('first own lags: true %s, least squares %s, posterior mean %s\n', ...
    mat2str(diag(Atrue(2:n+1, :))', 3), mat2str(diag(Ols(2:n+1, :))', 3), ...
    mat2str(diag(Ahat(2:n+1, :))', 3));
mods = zeros(R, 1);
for r = 1:R
    mods(r) = max(abs(eig([Adraw(2:end, :, r)'; eye(n*(p-1)), zeros(n*(p-1), n)])));
end
qst = quantile(paths(:, :, mods < 1), [.05 .16 .5 .84 .95], 3);
fprintf(['draws with largest modulus above one: %d of %d; without them, the largest change ' ...
    'in a quantile %.3f\n'], ...
    sum(mods >= 1), R, max(abs(qst - q), [], 'all'));
[mtrue, sdtrue] = stationary(Atrue, Sigtrue); [mhat, sdhat] = stationary(Ahat, Sigbar);
fprintf('means: true %s, at the posterior means %s, sample %s\n', mat2str(mtrue, 3), ...
    mat2str(mhat, 3), mat2str(mean(data), 3));
fprintf('sds: true %s, at the posterior means %s, sample %s\n', mat2str(sdtrue, 3), ...
    mat2str(sdhat, 3), mat2str(std(data), 3));
Wb = q(:, :, 5) - q(:, :, 1);                     % the widths of the 90% intervals
Wp = qplug(:, :, 5) - qplug(:, :, 1);
Wt = qtrue(:, :, 5) - qtrue(:, :, 1);
for h = [1 4 8 20]
    fprintf('\nh = %2d                 %26s %26s %26s\n', h, names{:});
    fprintf('  BVAR    median, 90%%:  %s\n', fmt3(squeeze(q(h, :, [3 1 5]))));
    fprintf('  true    median, 90%%:  %s\n', fmt3(squeeze(qtrue(h, :, [3 1 5]))));
    fprintf('  plug-in median, 90%%:  %s\n', fmt3(squeeze(qplug(h, :, [3 1 5]))));
    fprintf('  90%% widths, BVAR / plug-in / true: %s\n', sprintf('%7.2f %6.2f %6.2f   ', ...
        [Wb(h, :); Wp(h, :); Wt(h, :)]));
end
fprintf('\nlargest gap between the mean and the median of the paths: %.3f\n', ...
    max(abs(fmean - q(:, :, 3)), [], 'all'));
fprintf('largest gap between the BVAR and true medians over the %d horizons: %s\n', H, ...
    mat2str(max(abs(q(:, :, 3) - qtrue(:, :, 3))), 3));
fprintf('BVAR median of GDP growth below the true one at every horizon: %d\n', ...
    all(q(:, 1, 3) < qtrue(:, 1, 3)));
fprintf('90%% widths, BVAR over plug-in: h = 1 %s, h = 20 %s\n', mat2str(Wb(1, :)./Wp(1, :), 3), ...
    mat2str(Wb(H, :)./Wp(H, :), 3));
fprintf('90%% widths, BVAR over true: h = 1 %s, h = 20 %s\n', mat2str(Wb(1, :)./Wt(1, :), 3), ...
    mat2str(Wb(H, :)./Wt(H, :), 3));
fprintf('90%% widths, plug-in over true: h = 1 %s, h = 20 %s\n', mat2str(Wp(1, :)./Wt(1, :), 3), ...
    mat2str(Wp(H, :)./Wt(H, :), 3));
fprintf('probability of negative GDP growth next quarter: BVAR %.3f, true %.3f, plug-in %.3f\n', ...
    pneg, mean(ptrue(1, 1, :) < 0), mean(pplug(1, 1, :) < 0));
x = [1, reshape(ylast(end:-1:1, :)', 1, [])];      % x_{T+1} = (1, y_T', ..., y_{T-3}')
m1 = squeeze(pagemtimes(x, Adraw));               % n x R: the one-step means, one per draw
fprintf('one-step means: true %s, posterior mean %s, posterior sd %s\n', mat2str(x*Atrue, 3), ...
    mat2str(mean(m1, 2)', 3), mat2str(std(m1, 0, 2)', 3));

%% 6. the checks
% the one-step predictive density is multivariate t with nuhat - n + 1 degrees of freedom
df = nuhat - n + 1;
loc = x*Ahat; scl = sqrt(diag(Shat)'*(1 + x*(KA\x'))/df);
pr = [.05 .5 .95]';
qt = loc + scl.*tinv(pr, df);
q1 = squeeze(q(1, :, [1 3 5]))';
psd = squeeze(std(paths, 0, 3));                  % H x n, the predictive sds
gapt = max(abs(q1 - qt)./psd(1, :), [], 'all');
fprintf('\none-step quantiles against the t: largest gap %.3f predictive sds\n', gapt);
% the quantiles at every horizon against the mixture of normals from the same draws
[mu, sd] = bvar.forecast.predictive(reshape(Adraw, k*n, R)', permute(Sdraw, [3 1 2]), ...
    ylast(end:-1:1, :)', H);
gapq = zeros(H, n, 3);
for h = 1:H
    for i = 1:n
        gapq(h, i, :) = (bvar.forecast.mixquantile(mu(:, i, h), sd(:, i, h), pr) - ...
            squeeze(q(h, i, [1 3 5])))/psd(h, i);
    end
end
fprintf(['quantiles at every horizon against bvar.forecast.mixquantile: largest gap %.3f ' ...
    'predictive sds\n'], ...
    max(abs(gapq(:))));
% the draws: means against the posterior means of A and Sig
Asd = sqrt(diag(KA\eye(k))*diag(Sigbar)');       % k x n posterior sds of A
fprintf('mean of the draws of A against Ahat: largest gap %.3f posterior sds\n', ...
    max(abs(mean(Adraw, 3) - Ahat)./Asd, [], 'all'));
fprintf('mean of the draws of Sig against Shat/(nuhat - n - 1): largest relative gap %.4f\n', ...
    max(abs(mean(Sdraw, 3) - Sigbar)./abs(Sigbar), [], 'all'));
% the true model's paths against its exact predictive means and sds
[mut, sdt] = bvar.forecast.predictive(Atrue(:)', reshape(Sigtrue, [1 n n]), ylast(end:-1:1, :)', H);
mut = squeeze(mut)'; sdt = squeeze(sdt)';
fprintf('true paths against the exact predictive: means %.3f sds, sds %.4f (largest gaps)\n', ...
    max(abs(mean(ptrue, 3) - mut)./sdt, [], 'all'), max(abs(std(ptrue, 0, 3)./sdt - 1), [], 'all'));
assert(gapt < .1 && max(abs(gapq(:))) < .1, 'the simulated quantiles disagree with the exact ones');

%% 7. the probability of a recession: the adaptation
g = squeeze(paths(1:8, 1, :));                    % GDP growth, 8 quarters by R paths
prec = mean(any(g(1:7, :) < 0 & g(2:8, :) < 0));  % two negative quarters in a row
gt = squeeze(ptrue(1:8, 1, :)); gp = squeeze(pplug(1:8, 1, :));
fprintf(['\nprobability of a recession in the next eight quarters: BVAR %.3f, true %.3f, ' ...
    'plug-in %.3f\n'], ...
    prec, mean(any(gt(1:7, :) < 0 & gt(2:8, :) < 0)), mean(any(gp(1:7, :) < 0 & gp(2:8, :) < 0)));
fprintf('probability of negative growth in each of the eight quarters, BVAR: %s\n', ...
    mat2str(mean(g < 0, 2)', 2));

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 8]);
tiledlayout(1, n, 'TileSpacing', 'compact');
th = 201:240; tf = 240 + (0:H);
for i = 1:n
    nexttile; hold on
    h90 = bvar.util.shaded_band(tf, [data(end, i); q(:, i, 1)], [data(end, i); q(:, i, 5)], .85);
    h68 = bvar.util.shaded_band(tf, [data(end, i); q(:, i, 2)], [data(end, i); q(:, i, 4)], .7);
    hs = plot(th, data(th, i), 'k');
    hm = plot(tf, [data(end, i); q(:, i, 3)], 'k', 'LineWidth', 1);
    ht = plot(tf, [data(end, i); qtrue(:, i, 3)], 'k--');
    ht90 = plot(tf, [data(end, i); qtrue(:, i, 1)], 'k:');
    plot(tf, [data(end, i); qtrue(:, i, 5)], 'k:');
    box off; xlim([th(1) tf(end)]);
    title(names{i}, 'FontWeight', 'normal');
end
lg = legend([hs hm h68 h90 ht ht90], {'Sample', 'BVAR median', 'BVAR 68% interval', ...
    'BVAR 90% interval', 'True model median', 'True model 90% interval'}, 'NumColumns', 3, ...
    'Box', 'off');
lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(here, 'fig_fan.png'), 'Resolution', 150); close(fig);

fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 7]);
tiledlayout(1, n, 'TileSpacing', 'compact');
for i = 1:n
    nexttile; hold on
    hb = plot(1:H, Wb(:, i), 'k', 'LineWidth', 1);
    hp = plot(1:H, Wp(:, i), 'k-.');
    ht = plot(1:H, Wt(:, i), 'k--');
    box off; xlim([1 H]); ylim([0, 1.05*max([Wb(:, i); Wp(:, i); Wt(:, i)])]);
    xlabel('quarters ahead'); title(names{i}, 'FontWeight', 'normal');
end
lg = legend([hb hp ht], {'BVAR', 'Plug-in', 'True model'}, 'Orientation', 'horizontal', ...
    'Box', 'off');
lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(here, 'fig_widths.png'), 'Resolution', 150); close(fig);

function [m, s] = stationary(B, S)
% the means and standard deviations of the stationary distribution of the VAR with coefficients
% B, intercept first, and error covariance matrix S
[k, n] = size(B); p = (k - 1)/n;
F = [B(2:end, :)'; eye(n*(p-1)), zeros(n*(p-1), n)];
mu = (eye(n*p) - F)\[B(1, :)'; zeros(n*(p-1), 1)];
Q = blkdiag(S, zeros(n*(p-1)));
V = reshape((eye((n*p)^2) - kron(F, F))\Q(:), n*p, n*p);
m = mu(1:n)'; s = sqrt(diag(V(1:n, 1:n)))';
end

function s = fmt3(M)
% the rows of M, median, lower and upper bound, as 'median (lower, upper)'
s = '';
for i = 1:size(M, 1)
    s = [s, sprintf('%26s', sprintf('%.2f (%.2f, %.2f)', M(i, 1), M(i, 2), M(i, 3)))]; %#ok<AGROW>
end
end
