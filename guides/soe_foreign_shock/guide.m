%% guide - the code of guides/soe_foreign_shock/README.md: the responses of a small open economy,
% Australia with a US block, to a shock to the US excess bond premium, in the VAR with block
% exogeneity of guides/soe_var. The shock is identified recursively in the US block, with the
% premium after GDP growth and inflation and before the federal funds rate, as in Gilchrist and
% Zakrajsek (2012). The guide computes the responses and the shares of the forecast error
% variances for every draw from the posterior of the two blocks and compares them with those of
% the true model. The checks compare the responses with a simulation of the VAR, with the other
% order of the premium and the funds rate, and confirm that no Australian shock moves a US
% variable. Section 4 and the function responses at the end are the code to adapt. Saves
% fig_responses.png next to this file.

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
[Atrue, Sigtrue] = reduced_form(Bus, Bau, Sus, Om, ius, k);

%% 2. the sample of guides/soe_var
q0 = max(p, 4);                           % presample quarters: p lags, four for the AR(4) fits
Yfull = bvar.util.simulate_var(Atrue, Sigtrue, T + q0);     % the presample, then y_1, ..., y_T

%% 3. draws from the posterior of the two blocks, as in guides/soe_var
R = 10000;
[A_d, Sig_d, c1] = block_draws(Yfull, p, us, au, lev, c3, R);

%% 4. the responses to a shock to the excess bond premium: the code to adapt
H = 21;                                   % horizons 0 to 20
ord = [1 2 4 3 5:n];                      % the recursive order: the premium before the funds rate
js = 3;                                   % the premium shock in that order
ig = [1 2 5 6 8];                         % the variables in growth rates
resp_true = responses(Atrue, Sigtrue, ord, H);          % all shocks of the true model
irf_true = level(squeeze(resp_true(:, js, :)), ig);
share = bvar.structural.fevd(resp_true);                 % variance shares at every horizon
fev_true = share(:, js, 8);                              % the premium shock, 8 quarters ahead
irf_d = zeros(n, H, R); fev_d = zeros(n, R); us_au = 0;
for r = 1:R
    resp = responses(A_d(:, :, r), Sig_d(:, :, r), ord, H);
    irf_d(:, :, r) = level(squeeze(resp(:, js, :)), ig);
    share = bvar.structural.fevd(resp);
    fev_d(:, r) = share(:, js, 8);                      % the premium shock, 8 quarters ahead
    us_au = max(us_au, max(abs(resp(us, au, :)), [], 'all'));   % US variables, Australian shocks
end
qirf = quantile(irf_d, [.05 .5 .95], 3);  % the median and the 90 percent band
qfev = quantile(fev_d, [.05 .5 .95], 2);

%% 5. a shock of one percentage point: the adaptation
irf_pp = irf_d./irf_d(4, 1, :);           % each draw divided by its response of the premium
qpp = quantile(irf_pp, [.05 .5 .95], 3);
irf_pp_true = irf_true/irf_true(4, 1);

%% the results
show = [1 5 9 21];
fprintf(['\nthe sample: c1 = %.3f; a shock of one standard deviation raises the premium by ' ...
    '%.2f\n'], c1, irf_true(4, 1));
fprintf(['responses after 0, 4, 8 and 20 quarters: true | posterior median (90%% band); GDP, ' ...
    'CPI and the dollar in percent of the level\n']);
for i = 1:n
    fprintf('  %-20s', names{i});
    for h = show
        fprintf(' %6.2f | %6.2f (%6.2f, %6.2f)', irf_true(i, h), qirf(i, h, 2), qirf(i, h, 1), ...
            qirf(i, h, 3));
    end
    fprintf('\n');
end
[~, hpk] = max(irf_true(8, :));
fprintf('the true response of the dollar peaks at %.2f after %d quarters\n', irf_true(8, hpk), ...
    hpk - 1);
out = irf_true < qirf(:, :, 1) | irf_true > qirf(:, :, 3);   % true responses outside the band
fprintf('horizons, of %d, at which the band misses the true response, by variable: %s\n', H, ...
    mat2str(sum(out, 2)'));
for i = find(any(out, 2))'
    fprintf('  %-20s after %s quarters\n', names{i}, mat2str(find(out(i, :)) - 1));
end
fprintf('\nshare of the 8-quarter forecast error variance due to the shock, percent:\n');
for i = 1:n
    fprintf('  %-20s true %5.1f | posterior median %5.1f (%5.1f, %5.1f)\n', names{i}, ...
        100*fev_true(i), 100*qfev(i, 2), 100*qfev(i, 1), 100*qfev(i, 3));
end
fprintf('\na shock of one point to the premium, after 4 quarters: true | median (band)\n');
for i = [1 8 5]
    fprintf('  %-20s %6.2f | %6.2f (%6.2f, %6.2f)\n', names{i}, irf_pp_true(i, 5), ...
        qpp(i, 5, 2), qpp(i, 5, 1), qpp(i, 5, 3));
end

%% the checks
% the responses of the true model against a simulation of the VAR after the shock
L = chol(Sigtrue(ord, ord), 'lower'); L(ord, :) = L;
ysim = zeros(H + p, n); ysim(p + 1, :) = L(:, js)';
for h = 2:H
    ysim(p + h, :) = reshape(ysim(p+h-1:-1:h, :)', 1, [])*Atrue(2:end, :);
end
gap_sim = max(abs(squeeze(resp_true(:, js, :)) - ysim(p+1:end, :)'), [], 'all');
% the other order of the premium and the funds rate
resp_alt = responses(Atrue, Sigtrue, 1:n, H);            % the premium after the funds rate
irf_alt = level(squeeze(resp_alt(:, 4, :)), ig);
gap_ord = max(abs(irf_alt - irf_true), [], 2);
fprintf('\nresponses against a simulation of the VAR: largest gap %.1e\n', gap_sim);
fprintf('the premium ordered after the funds rate: largest change in a response, by variable:\n');
fprintf('  %s\n', mat2str(round(gap_ord', 2)));
fprintf(['responses of the US variables to the Australian shocks, largest over the draws: ' ...
    '%.1e\n'], us_au);
assert(gap_sim < 1e-10 && us_au == 0, 'the responses fail a check');

%% the figure: the responses to the shock, the posterior against the true model
fig = figure('Visible', 'off', 'Position', [100 100 900 430]);
tiledlayout(2, 4, 'TileSpacing', 'compact');
hz = 0:H-1;
for i = [4 1 2 3 8 5 6 7]
    nexttile; hold on
    hb = bvar.util.shaded_band(hz, qirf(i, :, 1), qirf(i, :, 3));
    plot(hz([1 end]), [0 0], 'k-', 'LineWidth', .5);
    hm = plot(hz, qirf(i, :, 2), 'k', 'LineWidth', 1);
    ht = plot(hz, irf_true(i, :), 'k--', 'LineWidth', 1);
    box off; xlim(hz([1 end])); title(names{i}, 'FontWeight', 'normal');
    if i >= 5, xlabel('quarters after the shock'); end
end
lg = legend([ht hm hb], {'true response', 'posterior median', '90 percent band'}, ...
    'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(here, 'fig_responses.png'), 'Resolution', 150);
close(fig);

%% local functions
function resp = responses(A, Sig, ord, H)
% the responses to shocks of one standard deviation identified recursively in the order ord:
% resp(i, j, h + 1) is the response of variable i, in the order of the data, to shock j of that
% order after h quarters
L = chol(Sig(ord, ord), 'lower');        % the impact matrix in the recursive order
L(ord, :) = L;                           % its rows back in the order of the data
resp = bvar.structural.irf_redu(A(2:end, :), L, H, size(Sig, 1));
end

function x = level(x, ig)
% the responses of the variables in rows ig, which are 400 times log differences, as responses
% of the level in log points times 100
x(ig, :) = cumsum(x(ig, :), 2)/4;
end

function [A_d, Sig_d, c1] = block_draws(Yfull, p, us, au, lev, c3, R)
% R draws of the reduced form from the posterior of the two blocks, as in guides/soe_var, from a
% sample whose first max(p, 4) rows are the presample; c1 maximizes the marginal likelihood
nus = numel(us); nau = numel(au); n = nus + nau;
q0 = max(p, 4);
[Ys, Zs] = bvar.util.build_lags(Yfull(q0-p+1:end, :), p);
[T, k] = size(Zs);
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];
s2 = bvar.priors.resid_var_ar4(Yfull(1:q0, :), Ys);
lc = fminbnd(@(lc) -soe_lml(Ys, Zs, us, au, lev, exp(lc), c3, s2), log(1e-3), log(10));
c1 = exp(lc);
[VA, B0] = soe_prior(p, n, lev, c1, c3, s2);
[Bh_us, K_us, Sh_us] = ncp_post(Ys(:, us), Zs(:, ius), B0(ius, us), VA(ius), diag(s2(us)));
[Bh_au, K_au, Sh_au] = ncp_post(Ys(:, au), [Zs, Ys(:, us)], [B0(:, au); zeros(nus, nau)], ...
    [VA; 1./s2(us)], diag(s2(au)));
Cus = chol(K_us, 'lower'); Cau = chol(K_au, 'lower');
A_d = zeros(k, n, R); Sig_d = zeros(n, n, R);
for r = 1:R
    Sus = iwishrnd(Sh_us, nus + 2 + T);                 % the US block
    Bus = Bh_us + (Cus'\randn(numel(ius), nus))*chol(Sus, 'lower')';
    Om = iwishrnd(Sh_au, n + 2 + T);                    % Australia, independent of the US block
    Bau = Bh_au + (Cau'\randn(k + nus, nau))*chol(Om, 'lower')';
    [A_d(:, :, r), Sig_d(:, :, r)] = reduced_form(Bus, Bau, Sus, Om, ius, k);
end
end

function [VA, B0] = soe_prior(p, n, lev, c1, c3, s2)
% the prior of the coefficients on x_t = (1, y_{t-1}', ..., y_{t-p}')': their variances VA, and
% their means B0, one on the own first lag of the variables lev, in levels, and zero elsewhere
VA = [c3; kron(1./(1:p)'.^2, c1./s2)];
B0 = zeros(1 + n*p, n);
B0(sub2ind(size(B0), 1 + lev, lev)) = 1;
end

function lml = soe_lml(Y, Z, us, au, lev, c1, c3, s2)
% the log marginal likelihood of the VAR with block exogeneity: the sum over the two blocks
n = size(Y, 2); p = (size(Z, 2) - 1)/n; nus = numel(us); nau = numel(au);
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];
[VA, B0] = soe_prior(p, n, lev, c1, c3, s2);
lml = ncp_lml(Y(:, us), Z(:, ius), B0(ius, us), VA(ius), nus + 2, diag(s2(us))) ...
    + ncp_lml(Y(:, au), [Z, Y(:, us)], [B0(:, au); zeros(nus, nau)], [VA; 1./s2(us)], ...
    n + 2, diag(s2(au)));
end

function [B, K, Shat] = ncp_post(Y, X, B0, VA, S0)
% the posterior of the multivariate regression Y = X*B + E under vec(B) | Sig ~ N(vec(B0), Sig
% kron diag(VA)) and Sig ~ IW(nu0, S0): vec(B) | Sig, Y ~ N(vec(B), Sig kron inv(K)) and
% Sig | Y ~ IW(nu0 + T, Shat)
K = diag(1./VA) + X'*X;
B = K\(B0./VA + X'*Y);
E = Y - X*B; D = B - B0;
Shat = S0 + E'*E + D'*(D./VA);           % S0 + Y'*Y + B0'*(B0./VA) - B'*K*B, computed stably
Shat = (Shat + Shat')/2;
end

function lml = ncp_lml(Y, X, B0, VA, nu0, S0)
% the log marginal likelihood of the multivariate regression Y = X*B + E under the prior of
% ncp_post, in closed form, as in guides/natural_conjugate_ml
[T, m] = size(Y);
[~, K, Shat] = ncp_post(Y, X, B0, VA, S0);
lml = -m*T/2*log(pi) - m/2*(sum(log(VA)) + bvar.util.ldet(K)) ...
    + nu0/2*bvar.util.ldet(S0) - (nu0 + T)/2*bvar.util.ldet(Shat) ...
    + bvar.util.mgammaln(m, (nu0 + T)/2) - bvar.util.mgammaln(m, nu0/2);
end

function [A, Sig] = reduced_form(Bus, Bau, Sus, Om, ius, k)
% the reduced form of the two blocks: the US block yf_t' = x_t(ius)'*Bus + uf_t', and Australia
% yd_t' = (x_t', yf_t')*Bau + e_t', with e_t independent of uf_t
nus = size(Bus, 2);
G0 = Bau(k+1:end, :)';                   % the current US coefficients
A = zeros(k, nus + size(Bau, 2));
A(ius, 1:nus) = Bus;                     % zeros for the Australian lags
A(:, nus+1:end) = Bau(1:k, :) + A(:, 1:nus)*G0';
Sig = [Sus, Sus*G0'; G0*Sus, G0*Sus*G0' + Om];
end
