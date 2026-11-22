%% guide - the code of guides/var_forecast_evaluation/README.md: the BVAR of guides/var_forecast
% and an AR(4), re-estimated at each quarter and scored out of sample. Section 2 and the
% function ncp_forecast at the end are the code to adapt. Saves fig_scores.png next to this
% file.

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
[Tall, n] = size(data);

%% 2. the recursive forecasts: the code to adapt
H = 4; R = 1000;
origins = 120:Tall-1;                             % forecasts made at the end of quarters 120 to 239
no = numel(origins);
Yout = NaN(no, n, H);                             % the outcomes, NaN after the end of the sample
F = NaN(no, n, H, 3); LP = F; PIT = F;            % forecasts and scores: BVAR, AR(4), true model
rng(3, 'twister');
for j = 1:no
    t = origins(j);
    Yt = data(1:t, :);                            % the data up to the origin, and nothing after it
    hmax = min(H, Tall - t);
    yout = NaN(H, n); yout(1:hmax, :) = data(t+1:t+hmax, :);
    Yout(j, :, :) = yout';
    [F(j, :, :, 1), LP(j, :, :, 1), PIT(j, :, :, 1)] = ncp_forecast(Yt, yout, p, H, R);
    for i = 1:n                                   % an AR(4) for each variable
        [F(j, i, :, 2), LP(j, i, :, 2), PIT(j, i, :, 2)] = ...
            ncp_forecast(Yt(:, i), yout(:, i), p, H, R);
    end
    [F(j, :, :, 3), LP(j, :, :, 3), PIT(j, :, :, 3)] = gauss_forecast(Atrue, Sigtrue, Yt, yout);
end

%% 3. the scores
E = Yout - F;                                     % forecast errors, no x n x H x 3
rmsfe = squeeze(sqrt(mean(E.^2, 1, 'omitnan')));  % n x H x 3
alpl = squeeze(mean(LP, 1, 'omitnan'));           % average log predictive likelihood
cov90 = squeeze(sum(PIT > .05 & PIT < .95, 1)./sum(~isnan(PIT), 1));   % the 90% coverage

%% the results
mods = {'BVAR', 'AR(4)', 'true'};
fprintf('\n%d forecast origins, the ends of quarters %d to %d\n', no, origins(1), origins(end));
for h = [1 4]
    fprintf('\nh = %d, %d forecasts %24s %14s %14s\n', h, sum(~isnan(Yout(:, 1, h))), names{:});
    for m = 1:3
        fprintf('  RMSFE, %-8s %25.3f %14.3f %14.3f\n', mods{m}, rmsfe(:, h, m));
    end
    for m = 1:3
        fprintf('  average log score, %-8s %13.3f %14.3f %14.3f\n', mods{m}, alpl(:, h, m));
    end
    for m = 1:3
        fprintf('  90%% coverage, %-8s %18.3f %14.3f %14.3f\n', mods{m}, cov90(:, h, m));
    end
    fprintf('  RMSFE over that of the AR(4): BVAR %s, true %s\n', ...
        mat2str(rmsfe(:, h, 1)'./rmsfe(:, h, 2)', 3), mat2str(rmsfe(:, h, 3)'./rmsfe(:, h, 2)', 3));
    fprintf('  average log score minus that of the AR(4): BVAR %s, true %s\n', ...
        mat2str(alpl(:, h, 1)' - alpl(:, h, 2)', 3), mat2str(alpl(:, h, 3)' - alpl(:, h, 2)', 3));
end
fprintf(['\na correct 90%% interval covers between %.3f and %.3f of %d independent outcomes, ' ...
    'with probability 0.95\n'], binoinv(.025, no, .9)/no, binoinv(.975, no, .9)/no, no);
c = cumsum(LP(:, :, 1, 1) - LP(:, :, 1, 2));      % cumulative log score gains of the BVAR, h = 1
fprintf(['cumulative log score gains of the BVAR over the AR(4), h = 1: to quarter 180 %s, ' ...
    'to quarter 240 %s\n'], mat2str(c(origins == 179, :), 3), mat2str(c(end, :), 3));

%% 4. the checks
% one quarter ahead, the predictive density is Student-t in closed form: against the mixture
gt = zeros(no, n, 2); gp = gt;                    % gaps in log densities and PITs: BVAR, AR(4)
for j = 1:no
    Yt = data(1:origins(j), :);
    [lpt, pitt] = ncp_onestep(Yt, Yout(j, :, 1), p);
    gt(j, :, 1) = lpt - LP(j, :, 1, 1); gp(j, :, 1) = pitt - PIT(j, :, 1, 1);
    for i = 1:n
        [lpt, pitt] = ncp_onestep(Yt(:, i), Yout(j, i, 1), p);
        gt(j, i, 2) = lpt - LP(j, i, 1, 2); gp(j, i, 2) = pitt - PIT(j, i, 1, 2);
    end
end
gavg = max(abs(mean(gt, 1)), [], 'all');          % the gap in the average log scores
[~, ix] = max(abs(gt(:))); [jw, iv, mw] = ind2sub(size(gt), ix);   % the largest gap
Yw = data(1:origins(jw), :); yw = squeeze(Yout(jw, :, :))'; kw = iv;
if mw == 2, Yw = Yw(:, iv); yw = yw(:, iv); kw = 1; end
rng(5, 'twister');
[~, lpw] = ncp_forecast(Yw, yw, p, H, 20000);     % the same forecast from 20,000 draws
lpt = ncp_onestep(Yw, yw(1, :), p);
fprintf(['\none quarter ahead, the closed form against the mixture of %d draws: average log ' ...
    'scores %.4f, single log densities %.4f, PITs %.4f (largest gaps)\n'], R, gavg, ...
    max(abs(gt(:))), max(abs(gp(:))));
fprintf('  the largest gap, an outcome with PIT %.4f, is %.4f with 20,000 draws\n', ...
    PIT(jw, iv, 1, mw), abs(lpt(kw) - lpw(kw, 1)));
% the simulation error of the mixture four quarters ahead, at the first origin
rng(6, 'twister');
[~, lp20, pit20] = ncp_forecast(data(1:origins(1), :), squeeze(Yout(1, :, :))', p, H, 20000);
fprintf(['first origin, four quarters ahead, %d draws against 20,000: log densities %.4f, ' ...
    'PITs %.4f (largest gaps)\n'], R, max(abs(lp20(:, H)' - LP(1, :, H, 1))), ...
    max(abs(pit20(:, H)' - PIT(1, :, H, 1))));
assert(gavg < .005 && abs(lpt(kw) - lpw(kw, 1)) < .02, ...
    'the mixture disagrees with the closed form');

%% 5. the PIT histogram: the adaptation
counts = zeros(n, 10);
for i = 1:n
    counts(i, :) = histcounts(PIT(:, i, 1, 1), 0:.1:1);   % BVAR, one quarter ahead, ten bins
end
fprintf('\nPITs of the BVAR one quarter ahead, counts in ten bins:\n'); disp(counts);
fprintf(['counts from %d to %d; one bin holds between %d and %d of %d uniform PITs with ' ...
    'probability 0.95\n'], min(counts(:)), max(counts(:)), binoinv(.025, no, .1), ...
    binoinv(.975, no, .1), no);

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 7]);
tiledlayout(1, n, 'TileSpacing', 'compact');
for i = 1:n
    nexttile; hold on
    hb = plot(origins + 1, cumsum(LP(:, i, 1, 1) - LP(:, i, 1, 2)), 'k', 'LineWidth', 1);
    ht = plot(origins + 1, cumsum(LP(:, i, 1, 3) - LP(:, i, 1, 2)), 'k--');
    yline(0, 'k:'); box off; xlim([origins(1) + 1, Tall]);
    title(names{i}, 'FontWeight', 'normal');
end
lg = legend([hb ht], {'BVAR minus AR(4)', 'True model minus AR(4)'}, 'Orientation', ...
    'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(here, 'fig_scores.png'), 'Resolution', 150); close(fig);

function [m, lp, pit] = ncp_forecast(Yt, yout, p, H, R)
% the forecasts of a VAR(p), or of an AR(p) when Yt has one column, under the natural conjugate
% prior with c1 chosen by the marginal likelihood, all estimated on Yt alone: the predictive means
% m, and the log predictive densities lp and the PITs of the outcomes yout (H x n, NaN where there
% is none), each n x H, from R draws of the parameters
[Y, Z] = bvar.util.build_lags(Yt, p);
[T, n] = size(Y); k = size(Z, 2);
s2 = bvar.priors.resid_var_ar4(Yt(1:4, :), Y);     % AR(4) residual variances, from Yt only
VA = @(c1) [10^2; kron(1./(1:p)'.^2, c1./s2)];
lc = fminbnd(@(lc) -bvar.ml.ncp(Y, Z, VA(exp(lc)), n + 2, diag(s2)), log(1e-3), log(10));
[~, Ahat, KA, Shat] = bvar.ml.ncp(Y, Z, VA(exp(lc)), n + 2, diag(s2));
LK = chol(KA, 'lower');
Ad = zeros(R, k*n); Sd = zeros(R, n, n);
for r = 1:R
    S = iwishrnd(Shat, n + 2 + T);
    A = Ahat + (LK'\randn(k, n))*chol(S, 'lower')';
    Ad(r, :) = A(:)'; Sd(r, :, :) = S;
end
[mu, sd] = bvar.forecast.predictive(Ad, Sd, Yt(end:-1:end-p+1, :)', H);   % R x n x H
z = (reshape(yout', 1, n, H) - mu)./sd;
m = reshape(mean(mu, 1), n, H);                    % the mean of the mixture over the draws
lp = reshape(bvar.util.logsumexp(-log(sd) - z.^2/2, 1) - log(R) - log(2*pi)/2, n, H);
pit = reshape(mean(normcdf(z), 1), n, H);
end

function [lp, pit] = ncp_onestep(Yt, y1, p)
% the log predictive density and the PIT of the outcome y1 (1 x n) one quarter ahead, from the
% Student-t predictive distribution of the same estimate as in ncp_forecast
[Y, Z] = bvar.util.build_lags(Yt, p);
[T, n] = size(Y);
s2 = bvar.priors.resid_var_ar4(Yt(1:4, :), Y);
VA = @(c1) [10^2; kron(1./(1:p)'.^2, c1./s2)];
lc = fminbnd(@(lc) -bvar.ml.ncp(Y, Z, VA(exp(lc)), n + 2, diag(s2)), log(1e-3), log(10));
[~, Ahat, KA, Shat] = bvar.ml.ncp(Y, Z, VA(exp(lc)), n + 2, diag(s2));
x = [1, reshape(Yt(end:-1:end-p+1, :)', 1, [])];
df = T + 3;                                        % nuhat - n + 1, with nuhat = n + 2 + T
scl = sqrt(diag(Shat)'*(1 + x*(KA\x'))/df);
z = (y1 - x*Ahat)./scl;
lp = log(tpdf(z, df)) - log(scl);
pit = tcdf(z, df);
end

function [m, lp, pit] = gauss_forecast(A, Sig, Yt, yout)
% the same for the VAR with known coefficients A and error covariance matrix Sig
[k, n] = size(A); p = (k - 1)/n; H = size(yout, 1);
[mu, sd] = bvar.forecast.predictive(A(:)', reshape(Sig, [1 n n]), Yt(end:-1:end-p+1, :)', H);
m = reshape(mu, n, H); sd = reshape(sd, n, H);
z = (yout' - m)./sd;
lp = -log(sd) - z.^2/2 - log(2*pi)/2;
pit = normcdf(z);
end
