%% guide - the code of guides/var_lag_shrinkage/README.md: the lag length and the shrinkage of the
% BVAR of guides/var_forecast, chosen by the marginal likelihood at each forecast origin. Section
% 2 and the functions at the end are the code to adapt. Saves fig_choice.png next to this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));

%% 1. the true model of guides/var_simulation, and the sample of 240 quarters simulated there
D = readmatrix(fullfile(fileparts(here), 'data', 'macro3_Q.csv'));
[Yd, Zd] = bvar.util.build_lags(D(:, 2:4), 4);
Atrue = Zd\Yd;
Ud = Yd - Zd*Atrue;
Sigtrue = Ud'*Ud/(size(Yd, 1) - size(Zd, 2));
rng(2, 'twister');
data = bvar.util.simulate_var(Atrue, Sigtrue, 240);   % GDP growth, inflation, unemployment
[Tall, n] = size(data);

%% 2. the choices at each origin: the code to adapt
pmax = 8; c1fix = [.2^2, 10];                     % c1 = 0.2^2, and a nearly flat prior
origins = 120:Tall-1; no = numel(origins);
LS = NaN(no, pmax, 3, 2);                         % joint log scores: p, the c1, the window
SE = NaN(no, n, pmax, 3, 2);                      % squared forecast errors
pstar = zeros(no, 2); c1star = zeros(no, 2);      % the choices of the marginal likelihood
for w = 1:2
    for j = 1:no
        t = origins(j);
        if w == 1, W = data(1:t, :); else, W = data(t-79:t, :); end   % all data, or the last 80
        [LS(j, :, :, w), SE(j, :, :, :, w), pstar(j, w), c1star(j, w)] = ...
            choose_and_forecast(W, data(t+1, :), pmax, c1fix);
    end
end

%% 3. the scores
als = squeeze(mean(LS, 1));                       % average joint log scores, pmax x 3 x 2
rmsfe = squeeze(sqrt(mean(SE, 1)));               % n x pmax x 3 x 2
lsml = zeros(no, 2);                              % the scores of the choices of the ML
for w = 1:2
    lsml(:, w) = LS(sub2ind(size(LS), (1:no)', pstar(:, w), ones(no, 1), w*ones(no, 1)));
end
lst = zeros(no, 1);                               % the scores of the true model
for j = 1:no
    t = origins(j);
    x = [1, reshape(data(t:-1:t-3, :)', 1, [])];
    lst(j) = mvt_logpdf(data(t+1, :), x*Atrue, Sigtrue, Inf);
end

%% the results
wn = {'all data up to the origin', 'the last 80 quarters'};
for w = 1:2
    fprintf('\n%s: average joint log scores by lag length\n', wn{w});
    fprintf('  p                 %s\n', sprintf('%8d', 1:pmax));
    fprintf('  c1 by the ML      %s\n', sprintf('%8.3f', als(:, 1, w)));
    fprintf('  c1 = 0.2^2        %s\n', sprintf('%8.3f', als(:, 2, w)));
    fprintf('  c1 = 10           %s\n', sprintf('%8.3f', als(:, 3, w)));
    [best, ib] = max(als(:, :, w), [], 'all', 'linear'); [pb, sb] = ind2sub([pmax 3], ib);
    fprintf(['  the choice of the marginal likelihood %.3f; the best alternative %.3f ' ...
        '(p = %d, c1 %d)\n'], mean(lsml(:, w)), best, pb, sb);
    fprintf('  chosen lag lengths: %s (counts of p = 1 to %d)\n', ...
        mat2str(histcounts(pstar(:, w), .5:pmax+.5)), pmax);
    fprintf('  chosen c1: median %.3f, from %.3f to %.3f\n', median(c1star(:, w)), ...
        min(c1star(:, w)), max(c1star(:, w)));
    for p = [4 8]
        fprintf('  RMSFE, p = %d: c1 by the ML %s, c1 = 0.2^2 %s, c1 = 10 %s\n', p, ...
            mat2str(rmsfe(:, p, 1, w)', 3), mat2str(rmsfe(:, p, 2, w)', 3), ...
            mat2str(rmsfe(:, p, 3, w)', 3));
    end
end
fprintf('\nthe true model: average joint log score %.3f\n', mean(lst));
fprintf('last 80 quarters, eight lags: RMSFE with c1 = 10 over that with c1 chosen %s\n', ...
    mat2str(rmsfe(:, 8, 3, 2)'./rmsfe(:, 8, 1, 2)', 3));

%% 4. the checks
% the log marginal likelihood is the sum of the one-step log predictive densities of the
% observations, each from the posterior given the observations before it
W = data(1:origins(1), :); p = 4;
s2 = bvar.priors.resid_var_ar4(W(pmax-3:pmax, :), W(pmax+1:end, :));
[Y, Z] = bvar.util.build_lags(W(pmax-p+1:end, :), p);
VA = [10^2; kron(1./(1:p)'.^2, c1star(1, 1)./s2)];
T = size(Y, 1); S0 = diag(s2); nu0 = n + 2;
lsum = 0;
for t = 1:T
    if t == 1                                     % the prior predictive of the first observation
        Ahat = zeros(size(Z, 2), n); KA = diag(1./VA); Shat = S0;
    else
        [~, Ahat, KA, Shat] = bvar.ml.ncp(Y(1:t-1, :), Z(1:t-1, :), VA, nu0, S0);
    end
    df = nu0 + t - 1 - n + 1;
    x = Z(t, :);
    lsum = lsum + mvt_logpdf(Y(t, :), x*Ahat, (1 + x*(KA\x'))*Shat/df, df);
end
gapml = lsum - bvar.ml.ncp(Y, Z, VA, nu0, S0);
fprintf(['\nfirst origin, p = 4: log ML %.6f, sum of %d one-step log predictive densities ' ...
    '%.6f, gap %.1e\n'], bvar.ml.ncp(Y, Z, VA, nu0, S0), T, lsum, gapml);
assert(abs(gapml) < 1e-8, ...
    'the marginal likelihood is not the product of the one-step predictive densities');

%% 5. forecasts averaged over c1: the adaptation
lmode = .2; lsd = .4;                             % lambda = sqrt(c1): gamma, mode 0.2, sd 0.4
th = (-lmode + sqrt(lmode^2 + 4*lsd^2))/2; a = lmode/th + 1;   % its scale and shape
grid = exp(linspace(log(1e-3), log(10), 81));     % c1, evenly spaced in log c1
lpri = log(gampdf(sqrt(grid), a, th)) + log(sqrt(grid)/2);   % the prior density of log c1
lsav = zeros(no, 2);
for w = 1:2
    for j = 1:no
        t = origins(j);
        if w == 1, W = data(1:t, :); else, W = data(t-79:t, :); end
        lsav(j, w) = average_over_c1(W, data(t+1, :), 4, pmax, grid, lpri);
    end
    fprintf(['%s, p = 4: average joint log score with c1 chosen %.4f, averaged over c1 %.4f, ' ...
        'difference %.4f\n'], wn{w}, als(4, 1, w), mean(lsav(:, w)), ...
        mean(lsav(:, w)) - als(4, 1, w));
end

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 8]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
yl = [min([als(:); mean(lst)]) - .03, max([als(:); mean(lst)]) + .03];
for w = 1:2
    nexttile; hold on
    hm = plot(1:pmax, als(:, 1, w), 'k-o', 'LineWidth', 1, 'MarkerSize', 4);
    hs = plot(1:pmax, als(:, 2, w), 'k--');
    hf = plot(1:pmax, als(:, 3, w), 'k:', 'LineWidth', 1);
    ht = yline(mean(lst), 'Color', [.6 .6 .6]);
    box off; xlim([1 pmax]); ylim(yl); xlabel('lag length');
    title(['Estimated on ' wn{w}], 'FontWeight', 'normal');
end
lg = legend([hm hs hf ht], {'c_1 chosen by the marginal likelihood', 'c_1 = 0.2^2', ...
    'c_1 = 10', 'True model'}, 'NumColumns', 2, 'Box', 'off');
lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(here, 'fig_choice.png'), 'Resolution', 150); close(fig);

function [ls, se, pstar, c1star] = choose_and_forecast(W, y1, pmax, c1fix)
% for lag lengths 1 to pmax, the log marginal likelihood maximized over c1, and the joint log
% scores and the squared errors of the one-step forecasts of the outcome y1, with c1 chosen and
% with each value in c1fix, all estimated on the window W; pstar and c1star are the choices of the
% marginal likelihood
n = size(W, 2); ns = 1 + numel(c1fix);
ls = zeros(pmax, ns); se = zeros(n, pmax, ns);
lml = zeros(pmax, 1); c1 = zeros(pmax, 1);
s2 = bvar.priors.resid_var_ar4(W(pmax-3:pmax, :), W(pmax+1:end, :));   % AR(4) variances
for p = 1:pmax
    [Y, Z] = bvar.util.build_lags(W(pmax-p+1:end, :), p);    % the same Y for every p
    VA = @(c) [10^2; kron(1./(1:p)'.^2, c./s2)];
    [lc, f] = fminbnd(@(lc) -bvar.ml.ncp(Y, Z, VA(exp(lc)), n + 2, diag(s2)), log(1e-3), log(10));
    lml(p) = -f; c1(p) = exp(lc);
    x = [1, reshape(W(end:-1:end-p+1, :)', 1, [])];          % x_{t+1}
    cs = [c1(p), c1fix];
    for s = 1:ns
        [~, Ahat, KA, Shat] = bvar.ml.ncp(Y, Z, VA(cs(s)), n + 2, diag(s2));
        df = size(Y, 1) + 3;                                  % nuhat - n + 1
        m = x*Ahat; Psi = (1 + x*(KA\x'))*Shat/df;            % the Student-t forecast
        ls(p, s) = mvt_logpdf(y1, m, Psi, df);
        se(:, p, s) = (y1 - m).^2;
    end
end
[~, pstar] = max(lml);                            % the lag length with the largest ML
c1star = c1(pstar);
end

function ls = average_over_c1(W, y1, p, pmax, grid, lpri)
% the joint log score of the one-step forecast of y1 from a VAR(p) estimated on W, with the
% forecasts for the values of c1 in grid averaged by their posterior probabilities under the
% prior density lpri of log c1
n = size(W, 2);
s2 = bvar.priors.resid_var_ar4(W(pmax-3:pmax, :), W(pmax+1:end, :));
[Y, Z] = bvar.util.build_lags(W(pmax-p+1:end, :), p);
x = [1, reshape(W(end:-1:end-p+1, :)', 1, [])];
lw = zeros(size(grid)); lp = lw;
for g = 1:numel(grid)
    VA = [10^2; kron(1./(1:p)'.^2, grid(g)./s2)];
    [lw(g), Ahat, KA, Shat] = bvar.ml.ncp(Y, Z, VA, n + 2, diag(s2));
    df = size(Y, 1) + 3;
    lp(g) = mvt_logpdf(y1, x*Ahat, (1 + x*(KA\x'))*Shat/df, df);
end
lw = lw + lpri; lw = lw - bvar.util.logsumexp(lw(:));   % the posterior probabilities of the grid
ls = bvar.util.logsumexp(lw(:) + lp(:));          % the log of the averaged predictive density
end

function lp = mvt_logpdf(y, m, Psi, df)
% the log density at y of the multivariate Student-t with location m, scale matrix Psi and df
% degrees of freedom; df = Inf gives the Gaussian
n = numel(y);
C = chol(Psi, 'lower');
z = C\(y - m)';
if isinf(df)
    lp = -n/2*log(2*pi) - sum(log(diag(C))) - z'*z/2;
else
    lp = gammaln((df + n)/2) - gammaln(df/2) - n/2*log(df*pi) - sum(log(diag(C))) ...
        - (df + n)/2*log1p(z'*z/df);
end
end
