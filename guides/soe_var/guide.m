%% guide - the code of guides/soe_var/README.md: a VAR for a small open economy, Australia with a
% US block, in which the Australian variables do not enter the US equations, estimated block by
% block. Sections 1 and 3 are the code to adapt. Saves fig_us_coef.png next to this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
rng(1, 'twister');
names = {'US GDP growth', 'US CPI inflation', 'federal funds rate', 'excess bond premium', ...
    'AU GDP growth', 'AU CPI inflation', 'AU interest rate', 'AUD depreciation'};

%% 1. the VAR with block exogeneity fitted to the data: the code to adapt
D = readmatrix(fullfile(fileparts(here), 'data', 'soe_Q.csv'));
data = D(:, 2:9);                         % 1984Q1-2019Q4: the US block, then Australia
nus = 4; nau = 4; n = nus + nau;
us = 1:nus; au = nus+1:n;
lev = [3 4 7];                            % the variables in levels: the two rates and the premium
c3 = 10^2;                                % the prior variance of the intercepts
pmax = 6; lml = zeros(pmax, 1); c1q = zeros(pmax, 1);
for q = 1:pmax                            % the lag length and the shrinkage, on the same quarters
    rq = bvar.models.var_soe(data(1:pmax, :), data(pmax+1:end, :), q, nus, 'lev', lev, ...
        'c3', c3, 'nsim', 0);
    c1q(q) = rq.c1; lml(q) = rq.lml;      % the chosen c1 and its log marginal likelihood
end
[~, p] = max(lml);
[Y, Z] = bvar.util.build_lags(data, p);   % x_t = (1, y_{t-1}', ..., y_{t-p}')'
[T, k] = size(Z);
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];   % the intercept and the US lags in x_t
Xus = Z(:, ius);                          % the US block: its own lags
Xau = [Z, Y(:, us)];                      % Australia: all lags and the current US variables
Bus = Xus\Y(:, us);                       % least squares, block by block
Bau = Xau\Y(:, au);
Uus = Y(:, us) - Xus*Bus;  Sus = Uus'*Uus/(T - size(Xus, 2));
Uau = Y(:, au) - Xau*Bau;  Om = Uau'*Uau/(T - size(Xau, 2));
[Atrue, Sigtrue] = bvar.util.soe_reduced_form(Bus, Bau, Sus, Om);   % the reduced form

%% 2. a sample from the true model
F = [Atrue(2:end, :)'; eye(n*(p-1)), zeros(n*(p-1), n)];   % the companion matrix
rho = max(abs(eig(F)));
q0 = max(p, 4);                           % presample quarters: p lags, four for the AR(4) fits
Yfull = bvar.util.simulate_var(Atrue, Sigtrue, T + q0);     % the presample, then y_1, ..., y_T
[Ys, Zs] = bvar.util.build_lags(Yfull(q0-p+1:end, :), p);   % y_1, ..., y_T and their regressors

%% 3. the posterior of each block: the code to adapt
res = bvar.models.var_soe(Yfull(1:q0, :), Ys, p, nus, 'lev', lev, 'c3', c3, 'nsim', 0);
c1 = res.c1;                              % the shrinkage chosen by the marginal likelihood
s2 = res.s2;                              % AR(4) variances of the sample
VA = res.VA; B0 = res.B0;                 % the prior variances and means on x_t
VAus = VA(ius); B0us = B0(ius, us);       % the US block
VAau = [VA; 1./s2(us)];                   % Australia: the current US variables, not shrunk
B0au = [B0(:, au); zeros(nus, nau)];
nu0us = nus + 2; S0us = diag(s2(us));     % Sig_us ~ IW(nu0us, S0us)
nu0au = n + 2; S0au = diag(s2(au));       % Omega ~ IW(nu0au, S0au)
[~, Bhat_us, K_us, Shat_us] = bvar.ml.ncp(Ys(:, us), Zs(:, ius), VAus, nu0us, S0us, B0us);
[~, Bhat_au, K_au, Shat_au] = bvar.ml.ncp(Ys(:, au), [Zs, Ys(:, us)], VAau, nu0au, S0au, B0au);
R = 10000;                                % independent draws, block by block
Bus_d = zeros(numel(ius), nus, R); Sus_d = zeros(nus, nus, R);
Bau_d = zeros(k + nus, nau, R);    Om_d = zeros(nau, nau, R);
A_d = zeros(k, n, R);
Cus = chol(K_us, 'lower'); Cau = chol(K_au, 'lower');
st = rng;                                 % the state of the generator before the draws
for r = 1:R
    Sus_d(:, :, r) = iwishrnd(Shat_us, nu0us + T);      % the US block
    Bus_d(:, :, r) = Bhat_us + (Cus'\randn(numel(ius), nus))*chol(Sus_d(:, :, r), 'lower')';
    Om_d(:, :, r) = iwishrnd(Shat_au, nu0au + T);       % Australia, independent of the US block
    Bau_d(:, :, r) = Bhat_au + (Cau'\randn(k + nus, nau))*chol(Om_d(:, :, r), 'lower')';
    A_d(:, :, r) = bvar.util.soe_reduced_form(Bus_d(:, :, r), Bau_d(:, :, r), Sus_d(:, :, r), ...
        Om_d(:, :, r));
end

% the unrestricted VAR, under the natural conjugate prior with the same VA and B0
nu0un = n + 2;
[~, Bhat_un, K_un, Shat_un] = bvar.ml.ncp(Ys, Zs, VA, nu0un, diag(s2), B0);
Cun = chol(K_un, 'lower');
Bun_d = zeros(k, n, R);
for r = 1:R
    Bun_d(:, :, r) = Bhat_un + (Cun'\randn(k, n))*chol(iwishrnd(Shat_un, nu0un + T), 'lower')';
end

%% the results
fprintf(['\nthe data: %d quarters; log marginal likelihood by lag length, 1 to %d, on the ' ...
    'last %d:\n'], size(data, 1), pmax, size(data, 1) - pmax);
fprintf('  %s\n', mat2str(round(lml', 2)));
fprintf('  with c1 = %s\n', mat2str(round(c1q', 3)));
fprintf(['the chosen lag length %d, %.2f above %d lags and %.2f above %d; %d quarters after ' ...
    'the first %d\n'], p, lml(p) - lml(p-1), p - 1, lml(p) - lml(p+1), p + 1, T, p);
fprintf('%d coefficients in each US equation, %d in each Australian one\n', numel(ius), k + nus);
fprintf('the true model: largest modulus of the companion matrix %.3f\n', rho);
fprintf('  own first lags of the variables in levels: %s\n', ...
    mat2str(round(Atrue(sub2ind(size(Atrue), 1 + lev, lev)), 2)));
fprintf('true error sd: %s\n', mat2str(round(sqrt(diag(Sigtrue))', 2)));
fprintf('current US coefficients in the Australian equations (rows AU GDP, CPI, rate, AUD):\n');
disp(round(Bau(end-nus+1:end, :)', 2));
r04 = bvar.models.var_soe(Yfull(1:q0, :), Ys, p, nus, 'lev', lev, 'c3', c3, 'c1', .04, 'nsim', 0);
fprintf('the sample: c1 = %.3f by the marginal likelihood, %.2f log points above c1 = 0.04\n', ...
    c1, res.lml - r04.lml);

iaulag = setdiff(2:k, ius);              % the Australian lags in x_t
cun = Bun_d(iaulag, us, :);              % their coefficients in the US equations
mun = mean(cun, 3); qun = quantile(cun, [.05 .95], 3);
nout = sum(qun(:, :, 1) > 0 | qun(:, :, 2) < 0, 'all');
fprintf('\nthe US equations of the unrestricted VAR: %d coefficients on Australian lags\n', ...
    numel(mun));
fprintf('  posterior means from %.2f to %.2f; %d of the 90%% intervals exclude 0\n', ...
    min(mun(:)), max(mun(:)), nout);
% the mean absolute error of the posterior means over the coefficients of the US equations that
% both VARs estimate: the intercepts and the coefficients on the US lags
mae = @(B) mean(abs(B(:) - reshape(Atrue(ius, us), [], 1)));
Ahat = mean(A_d, 3); Ahat_un = mean(Bun_d, 3);
fprintf(['US equations, mean absolute error of the posterior means over the %d coefficients ' ...
    'of both VARs: restricted %.3f, unrestricted %.3f\n'], numel(ius)*nus, mae(Ahat(ius, us)), ...
    mae(Ahat_un(ius, us)));
% the shortcut that Zha (1999, pp. 298-299) calls incorrect: unrestricted draws with the
% Australian lags of the US equations set to zero, the other coefficients left as they are
gap = (Ahat_un(ius, us) - Ahat(ius, us))./std(A_d(ius, us, :), 0, 3);
fprintf(['US equations, the other coefficients: zeroed unrestricted against restricted, ' ...
    'largest gap %.2f posterior sd\n'], max(abs(gap(:))));

%% the checks
% least squares by block against the first-order conditions of maximum likelihood of the
% restricted VAR, whose covariance matrices divide by T: GLS on the stacked equations at the
% covariance matrix of the block fit returns the coefficients of the block fit, and those
% coefficients return its covariance matrix
[Ablock, Sblock] = bvar.util.soe_reduced_form(Bus, Bau, Uus'*Uus/T, Uau'*Uau/T);
Agls = restricted_gls(Y, Z, ius, nus, Sblock);
Ugls = Y - Z*Agls;
fprintf(['\nleast squares by block against maximum likelihood: GLS coefficients %.1e, ' ...
    'covariance %.1e\n'], max(abs(Ablock(:) - Agls(:))), ...
    max(abs(Sblock - Ugls'*Ugls/T), [], 'all'));
% the draws against the posterior moments: E(B) = Bhat, E(Sig) = Shat/(nu0 + T - m - 1)
gm = @(Bd, Bh) max(abs(mean(Bd, 3) - Bh)./std(Bd, 0, 3), [], 'all');
gs = @(Sd, Sh, nu, m) max(abs(mean(Sd, 3) - Sh/(nu - m - 1)), [], 'all') ...
    /max(abs(Sh/(nu - m - 1)), [], 'all');
fprintf(['draws against the posterior: means %.3f and %.3f posterior sd; covariances %.1e ' ...
    'and %.1e relative\n'], gm(Bus_d, Bhat_us), gm(Bau_d, Bhat_au), ...
    gs(Sus_d, Shat_us, nu0us + T, nus), gs(Om_d, Shat_au, nu0au + T, nau));
% bvar.models.var_soe takes the same draws in one call, from the same state of the generator
rng(st);
rd = bvar.models.var_soe(Yfull(1:q0, :), Ys, p, nus, 'lev', lev, 'c3', c3, 'nsim', R);
gap_core = max(abs(rd.A - reshape(permute(A_d, [3 1 2]), R, [])), [], 'all');
fprintf('bvar.models.var_soe against these draws: largest gap %.1e\n', gap_core);
assert(gap_core == 0, 'bvar.models.var_soe takes other draws');
% the prior of the error covariance matrix: with Sig_us ~ IW(nus + 2, S0us), Omega ~ IW(n + 2,
% S0au) and the current US coefficients G0' | Omega ~ N(0, Omega kron diag(1./s2(us))), the
% matrix Sig of the reduced form is IW(n + 2, diag(s2)). The two log densities agree at the true
% values, with the Jacobian |Sig_us|^nau of the map from (Sig_us, G0, Omega) to Sig
lp_blocks = lniw(Sus, nu0us, S0us) + bvar.ml.lniwpdf(Bau(k+1:end, :), Om, zeros(nus, nau), ...
    sparse(diag(s2(us))), nu0au, S0au);
lp_full = lniw(Sigtrue, n + 2, diag(s2)) + nau*bvar.util.ldet(Sus);
fprintf('the prior of Sig from the two blocks against IW(n + 2, diag(s2)): gap %.1e\n', ...
    lp_blocks - lp_full);

%% the figure: coefficients of the US equations on the Australian lags, in the unrestricted VAR
fig = figure('Visible', 'off', 'Position', [100 100 760 330]);
x = 1:numel(mun);
hold on
plot([0, numel(mun) + 1], [0 0], 'k-', 'LineWidth', .5);
h1 = plot([x; x], [reshape(qun(:, :, 1), 1, []); reshape(qun(:, :, 2), 1, [])], '-', ...
    'Color', [.6 .6 .6], 'LineWidth', 1.5);
h2 = plot(x, mun(:)', 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
nlag = nau*p;
for j = 1:nus-1
    xline(j*nlag + .5, ':', 'Color', [.5 .5 .5]);
end
xlim([0, numel(mun) + 1]); box off
set(gca, 'XTick', (nlag + 1)/2 + (0:nus-1)*nlag, 'XTickLabel', names(us), 'FontSize', 8);
ylabel('coefficient');
legend([h2, h1(1)], {'posterior mean', '90% credible interval'}, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'Box', 'off');
exportgraphics(fig, fullfile(here, 'fig_us_coef.png'), 'Resolution', 150);
close(fig);

%% local functions
function l = lniw(Sig, nu, S)
% the log density of the inverse-Wishart distribution IW(nu, S) at Sig
m = size(Sig, 1);
l = nu/2*bvar.util.ldet(S) - nu*m/2*log(2) - bvar.util.mgammaln(m, nu/2) ...
    - (nu + m + 1)/2*bvar.util.ldet(Sig) - .5*trace(Sig\S);
end

function A = restricted_gls(Y, Z, ius, nus, Sig)
% generalized least squares on the stacked equations of the reduced-form VAR with zeros on the
% Australian lags of the US equations, given the error covariance matrix Sig
[T, n] = size(Y); k = size(Z, 2);
X = blkdiag(kron(eye(nus), Z(:, ius)), kron(eye(n - nus), Z));
W = kron(inv(Sig), speye(T));
b = (X'*W*X)\(X'*W*Y(:));
A = zeros(k, n);
A(ius, 1:nus) = reshape(b(1:numel(ius)*nus), [], nus);
A(:, nus+1:end) = reshape(b(numel(ius)*nus+1:end), k, []);
end
