%% guide - the code of guides/soe_conditional_forecast/README.md: forecasts of the Australian
% variables in the VAR of guides/soe_var, conditional on a path for the US variables. Section 4
% and the function cond_forecast at the end are the code to adapt. Saves fig_scenario.png next to
% this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
rng(1, 'twister');
names = {'US GDP', 'US CPI', 'federal funds rate', 'excess bond premium', ...
    'AU GDP', 'AU CPI', 'AU interest rate', 'Australian dollar'};

%% 1. the true model of guides/soe_var: the VAR with block exogeneity fitted to the data
D = readmatrix(fullfile(fileparts(here), 'data', 'soe_Q.csv'));
data = D(:, 2:9);                         % 1984Q1-2019Q4: the US block, then Australia
nus = 4; nau = 4; n = nus + nau;
us = 1:nus; au = nus+1:n;
lev = [3 4 7];                            % the variables in levels: the two rates and the premium
c3 = 10^2; p = 3;                         % the intercepts, and the lag length of guides/soe_var
[Y, Z] = bvar.util.build_lags(data, p);   % x_t = (1, y_{t-1}', ..., y_{t-p}')'
[T, k] = size(Z);
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];   % the intercept and the US lags in x_t
Xus = Z(:, ius); Xau = [Z, Y(:, us)];
Bus = Xus\Y(:, us);                       % least squares, block by block
Bau = Xau\Y(:, au);
Uus = Y(:, us) - Xus*Bus;  Sus = Uus'*Uus/(T - size(Xus, 2));
Uau = Y(:, au) - Xau*Bau;  Om = Uau'*Uau/(T - size(Xau, 2));
[Atrue, Sigtrue] = bvar.util.soe_reduced_form(Bus, Bau, Sus, Om);

%% 2. the sample of guides/soe_var
q0 = max(p, 4);                           % presample quarters: p lags, four for the AR(4) fits
Yfull = bvar.util.simulate_var(Atrue, Sigtrue, T + q0);     % the presample, then y_1, ..., y_T

%% 3. draws from the posterior of the two blocks, as in guides/soe_var
R = 10000;
res = bvar.models.var_soe(Yfull(1:q0, :), Yfull(q0+1:end, :), p, nus, 'lev', lev, 'c3', c3, ...
    'nsim', R);

%% 4. the forecasts given a US path: the code to adapt
H = 8;                                    % quarters ahead
ig = [1 2 4];                             % the Australian variables in growth rates
ylag = Yfull(end:-1:end-p+1, :)';         % the last p quarters, the most recent first
[mu_u, sd_u] = bvar.forecast.predictive(res.A, res.Sig, ylag, H);   % unconditional, R x n x H
yf_base = squeeze(mean(mu_u(:, us, :), 1))';   % the baseline US path: the mean forecast, H x nus
yf_scen = yf_base;
yf_scen(1:4, 4) = yf_scen(1:4, 4) + 1;    % the scenario: the premium one point higher for a year
m_base = zeros(nau, H, R); s_base = m_base; m_scen = m_base;
for r = 1:R
    A = reshape(res.A(r, :), k, n); Sig = squeeze(res.Sig(r, :, :));   % draw r
    [m_base(:, :, r), s_base(:, :, r)] = cond_forecast(A, Sig, ylag, yf_base, us, au);
    m_scen(:, :, r) = cond_forecast(A, Sig, ylag, yf_scen, us, au);
end
eff = level(m_scen - m_base, ig);         % the effect of the scenario, levels for growth rates
qeff = quantile(eff, [.05 .5 .95], 3);    % the median and the 90 percent credible band
eff_true = level(cond_forecast(Atrue, Sigtrue, ylag, yf_scen, us, au) ...
    - cond_forecast(Atrue, Sigtrue, ylag, yf_base, us, au), ig);
% the 90 percent predictive bands of the Australian variables, without and with the US path
wid = zeros(nau, H, 2);
for i = 1:nau
    for h = 1:H
        wid(i, h, 1) = diff(bvar.forecast.mixquantile(mu_u(:, au(i), h), sd_u(:, au(i), h), ...
            [.05 .95]));
        wid(i, h, 2) = diff(bvar.forecast.mixquantile(m_base(i, h, :), s_base(i, h, :), ...
            [.05 .95]));
    end
end

%% 5. the premium of 2008 and 2009: the adaptation
yf_gfc = yf_base;
yf_gfc(:, 4) = data(97:104, 4);           % the premium in 2008Q1-2009Q4
m_gfc = zeros(nau, H, R);
for r = 1:R
    A = reshape(res.A(r, :), k, n); Sig = squeeze(res.Sig(r, :, :));
    m_gfc(:, :, r) = cond_forecast(A, Sig, ylag, yf_gfc, us, au);
end
q_gfc = quantile(level(m_gfc - m_base, ig), [.05 .5 .95], 3);
gfc_true = level(cond_forecast(Atrue, Sigtrue, ylag, yf_gfc, us, au) ...
    - cond_forecast(Atrue, Sigtrue, ylag, yf_base, us, au), ig);

%% the results
fprintf(['\nthe sample: c1 = %.3f; the scenario raises the premium by one point in quarters 1 ' ...
    'to 4 and keeps the other US variables at their baseline\n'], res.c1);
fprintf('the baseline premium path: %s\n', mat2str(round(yf_base(:, 4)', 2)));
fprintf(['effect of the scenario after 1, 4 and 8 quarters: true | posterior median (90%% ' ...
    'credible band); GDP, CPI and the dollar in log points\n']);
for i = 1:nau
    fprintf('  %-20s', names{au(i)});
    for h = [1 4 8]
        fprintf(' %6.2f | %6.2f (%6.2f, %6.2f)', eff_true(i, h), qeff(i, h, 2), qeff(i, h, 1), ...
            qeff(i, h, 3));
    end
    fprintf('\n');
end
[~, hpk] = max(abs(eff_true), [], 2);
for i = 1:nau
    fprintf('  %-20s largest true effect %.2f after %d quarters\n', names{au(i)}, ...
        eff_true(i, hpk(i)), hpk(i));
end
fprintf(['a true depreciation of %.2f log points after 4 quarters: the Australian dollar buys ' ...
    '%.2f percent fewer US dollars\n'], eff_true(4, 4), 100*(1 - exp(-eff_true(4, 4)/100)));
out = eff_true < qeff(:, :, 1) | eff_true > qeff(:, :, 3);   % true effects outside the band
fprintf('horizons, of %d, at which the band misses the true effect, by variable: %s\n', H, ...
    mat2str(sum(out, 2)'));
fprintf(['\nwidth of the 90%% predictive band, as the variables enter the VAR, after 1, 4 and ' ...
    '8 quarters: without | with the US path (percent narrower)\n']);
for i = 1:nau
    fprintf('  %-20s', names{au(i)});
    for h = [1 4 8]
        fprintf(' %6.2f | %6.2f (%4.1f)', wid(i, h, 1), wid(i, h, 2), ...
            100*(1 - wid(i, h, 2)/wid(i, h, 1)));
    end
    fprintf('\n');
end
fprintf('\nthe premium of 2008Q1-2009Q4: %s\n', mat2str(round(yf_gfc(:, 4)', 2)));
fprintf('its effect after 4 and 8 quarters: true | posterior median (90%% credible band)\n');
for i = 1:nau
    fprintf('  %-20s', names{au(i)});
    for h = [4 8]
        fprintf(' %6.2f | %6.2f (%6.2f, %6.2f)', gfc_true(i, h), q_gfc(i, h, 2), ...
            q_gfc(i, h, 1), q_gfc(i, h, 3));
    end
    fprintf('\n');
end
out_gfc = gfc_true < q_gfc(:, :, 1) | gfc_true > q_gfc(:, :, 3);
fprintf('horizons, of %d, at which the band misses the true effect, by variable: %s\n', H, ...
    mat2str(sum(out_gfc, 2)'));
for i = find(any(out_gfc, 2))'
    fprintf('  %-20s after %s quarters\n', names{au(i)}, mat2str(find(out_gfc(i, :))));
end

%% the checks
% the conditional distribution of the joint forecast of the VAR, from its stacked mean and
% covariance over the H quarters, for the true model and the first five draws
gap_joint = 0;
for r = 0:5
    if r == 0
        A = Atrue; Sig = Sigtrue;
    else
        A = reshape(res.A(r, :), k, n); Sig = squeeze(res.Sig(r, :, :));
    end
    [mc, sc] = cond_forecast(A, Sig, ylag, yf_scen, us, au);
    [mj, sj] = joint_conditional(A, Sig, ylag, yf_scen, us, au);
    gap_joint = max([gap_joint, max(abs(mc - mj), [], 'all'), max(abs(sc - sj), [], 'all')]);
end
% given each draw's own forecast of the US variables, the conditional forecast equals the
% unconditional one, for the first 100 draws
gap_own = 0;
for r = 1:100
    A = reshape(res.A(r, :), k, n); Sig = squeeze(res.Sig(r, :, :));
    m_own = cond_forecast(A, Sig, ylag, squeeze(mu_u(r, us, :))', us, au);
    gap_own = max(gap_own, max(abs(m_own - squeeze(mu_u(r, au, :))), [], 'all'));
end
fprintf('\nagainst the conditional distribution of the joint forecast: largest gap %.1e\n', ...
    gap_joint);
fprintf(['given each draw''s own US forecast, against the unconditional forecast: largest gap ' ...
    '%.1e\n'], gap_own);
assert(gap_joint < 1e-8 && gap_own < 1e-8, 'the conditional forecast fails a check');

%% the figure: the effect of the scenario on the Australian variables
fig = figure('Visible', 'off', 'Position', [100 100 760 520]);
tiledlayout(2, 2, 'TileSpacing', 'compact');
hz = 0:H;
z = zeros(nau, 1);
for i = 1:nau
    nexttile; hold on
    hb = bvar.util.shaded_band(hz, [0, qeff(i, :, 1)], [0, qeff(i, :, 3)]);
    plot(hz([1 end]), [0 0], 'k-', 'LineWidth', .5);
    hm = plot(hz, [z(i), qeff(i, :, 2)], 'k', 'LineWidth', 1);
    ht = plot(hz, [z(i), eff_true(i, :)], 'k--', 'LineWidth', 1);
    box off; xlim(hz([1 end])); title(names{au(i)}, 'FontWeight', 'normal');
    if i >= 3, xlabel('quarters ahead'); end
end
lg = legend([ht hm hb], {'true effect', 'posterior median', '90 percent credible band'}, ...
    'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(here, 'fig_scenario.png'), 'Resolution', 150);
close(fig);

%% local functions
function [m, sd] = cond_forecast(A, Sig, ylag, yf, us, au)
% the mean m and standard deviation sd of the variables au, 1 to H quarters ahead, given the
% path yf (H x numel(us)) of the variables us, in a reduced-form VAR with block exogeneity: A is
% k x n with the intercept first, Sig is n x n and ylag is n x p, the most recent column first.
% The errors of au are G0 times those of us plus a remainder with covariance Om, so given the
% path the forecast of au is the regression of the Australian block run forward.
[n, p] = size(ylag);  H = size(yf, 1);  nd = numel(au);
G0 = Sig(au, us)/Sig(us, us);            % the projection of the Australian errors on the US errors
Om = Sig(au, au) - G0*Sig(us, au);       % the covariance of the remainder
Phi = reshape(A(2:end, :)', n, n, p);    % Phi(:, :, l): the coefficients on lag l
m = zeros(nd, H);  sd = zeros(nd, H);
Psi = zeros(nd, nd, H);  Psi(:, :, 1) = eye(nd);
yl = ylag;  V = zeros(nd);
for h = 1:H
    yh = A(1, :)';
    for l = 1:p, yh = yh + Phi(:, :, l)*yl(:, l); end    % the forecast one quarter ahead
    y = yh;
    y(us) = yf(h, :)';                                   % the US values on the path
    y(au) = yh(au) + G0*(y(us) - yh(us));                % the Australian mean given them
    yl = [y, yl(:, 1:p-1)];
    for l = 1:min(h-1, p), Psi(:, :, h) = Psi(:, :, h) + Phi(au, au, l)*Psi(:, :, h-l); end
    V = V + Psi(:, :, h)*Om*Psi(:, :, h)';
    m(:, h) = y(au);
    sd(:, h) = sqrt(diag(V));
end
end

function [m, sd] = joint_conditional(A, Sig, ylag, yf, us, au)
% the same moments from the joint Gaussian forecast of all n variables over H quarters, by the
% formulas of the conditional normal distribution: the check of cond_forecast
[n, p] = size(ylag);  H = size(yf, 1);
Phi = reshape(A(2:end, :)', n, n, p);
mu = zeros(n, H);  yl = ylag;
for h = 1:H
    yh = A(1, :)';
    for l = 1:p, yh = yh + Phi(:, :, l)*yl(:, l); end
    mu(:, h) = yh;  yl = [yh, yl(:, 1:p-1)];
end
Psi = zeros(n, n, H);  Psi(:, :, 1) = eye(n);          % the moving-average coefficients
for h = 2:H
    for l = 1:min(h-1, p), Psi(:, :, h) = Psi(:, :, h) + Phi(:, :, l)*Psi(:, :, h-l); end
end
M = zeros(n*H);
for h = 1:H
    for g = 1:h, M((h-1)*n+(1:n), (g-1)*n+(1:n)) = Psi(:, :, h-g+1); end
end
C = M*kron(eye(H), Sig)*M';                              % the covariance of the stacked forecast
iu = reshape(us' + (0:H-1)*n, [], 1);  ia = reshape(au' + (0:H-1)*n, [], 1);
x = reshape(yf', [], 1) - mu(iu);
mj = mu(ia) + C(ia, iu)*(C(iu, iu)\x);
Cj = C(ia, ia) - C(ia, iu)*(C(iu, iu)\C(iu, ia));
m = reshape(mj, numel(au), H);
sd = reshape(sqrt(diag(Cj)), numel(au), H);
end

function x = level(x, ig)
% the responses of the rows ig, which are 400 times log differences, as responses of the level
% in log points times 100; x may have pages, one per draw
x(ig, :, :) = cumsum(x(ig, :, :), 2)/4;
end
