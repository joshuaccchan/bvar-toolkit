%% guide - the code of guides/natural_conjugate_ml/README.md: the shrinkage of the VAR of
% guides/var_stacking chosen by the closed-form marginal likelihood under the natural conjugate
% prior. Section 2 is the code to adapt. Saves fig_marglik.png next to this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
rng(1, 'twister');

%% 1. the data of the earlier guides
n = 3; p = 2; T = 200;
a = [.5; .2; .1];                    % intercepts
A1 = [.5 .1 0; .2 .4 .1; 0 .2 .5];   % lag 1
A2 = [.2 0 0; 0 .1 0; .1 0 .2];      % lag 2
Sig = [1 .3 .2; .3 .8 .1; .2 .1 .6];
A = [a'; A1'; A2'];   % k x n, A = (a, A_1, A_2)'
C = chol(Sig, 'lower'); burn = 100;
Yfull = zeros(burn + p + T, n);
for t = p+1:burn+p+T
    Yfull(t,:) = a' + Yfull(t-1,:)*A1' + Yfull(t-2,:)*A2' + (C*randn(n,1))';
end
Yfull = Yfull(burn-1:end, :);   % y_{-3}, ..., y_0, y_1, ..., y_T
Y0 = Yfull(1:4, :);             % the presample of the AR(4) fits
[Y, Z] = bvar.util.build_lags(Yfull(3:end, :), p);
k = size(Z, 2);
Atrue = A; Sigtrue = Sig;
[~, ~, s2] = bvar.priors.minn(p, .2^2, .1^2, 10^2, Y0, Y, 4);   % the AR(4) residual variances

%% 2. the choice of c1: the code to adapt
nu0 = n + 2; S0 = diag(s2);                    % Sig ~ IW(nu0, S0)
c3 = 10^2;                                     % the intercepts
VA = @(c1) [c3; kron(1./(1:p)'.^2, c1./s2)];   % the diagonal of V_A
lml = @(c1) bvar.ml.ncp(Y, Z, VA(c1), nu0, S0);
lc = fminbnd(@(lc) -lml(exp(lc)), log(1e-3), log(10));   % maximize over log(c1)
c1hat = exp(lc);

%% 3. the marginal likelihood over a grid, and the posterior means
c1grid = logspace(-3, 1, 161);
lmlgrid = arrayfun(lml, c1grid);
pathA = zeros(k, n, numel(c1grid));
for j = 1:numel(c1grid)
    [~, pathA(:, :, j)] = lml(c1grid(j));
end
Ols = (Z'*Z)\(Z'*Y);
[lml_hat, Ahat] = lml(c1hat);
[lml_std, Astd] = lml(.04);
mae = @(B) mean(abs(B(2:end, :) - Atrue(2:end, :)), 'all');   % over the 18 lag coefficients
[~, jmax] = max(lmlgrid);
fprintf('\nthe chosen c1: %.4f (grid maximum at %.4f), log ML %.2f\n', c1hat, c1grid(jmax), ...
    lml_hat);
fprintf('prior standard deviation of the first own lags, sqrt(c1): %.3f\n', sqrt(c1hat));
fprintf('log ML at c1 = 0.04: %.2f; at c1 = 0.001: %.2f; at c1 = 10: %.2f\n', lml_std, ...
    lml(1e-3), lml(10));
onepeak = all(diff(lmlgrid(1:jmax)) > 0) && all(diff(lmlgrid(jmax:end)) < 0);
fprintf('on the grid, the log ML rises up to its maximum and falls after it: %d\n', onepeak);
fprintf('at c1 = 10, lag coefficients against least squares, largest gap: %.1e\n', ...
    max(abs(pathA(2:end, :, end) - Ols(2:end, :)), [], 'all'));
assert(onepeak && abs(log10(c1hat/c1grid(jmax))) <= 4/160, ...
    'the maximum of fminbnd is not the single peak of the grid');
fprintf('first own lags: least squares %s, chosen c1 %s, c1 = 0.04 %s\n', ...
    mat2str(diag(Ols(2:n+1, :))', 3), mat2str(diag(Ahat(2:n+1, :))', 3), ...
    mat2str(diag(Astd(2:n+1, :))', 3));
fprintf('estimation error, mean absolute: least squares %.3f, chosen c1 %.3f, c1 = 0.04 %.3f\n', ...
    mae(Ols), mae(Ahat), mae(Astd));

%% 4. the same choice on 200 additional datasets
rng(2, 'twister');
R = 200; mc = zeros(R, 4);   % chosen c1, then the error of least squares, chosen c1, c1 = 0.04
for r = 1:R
    Yr = simulate_var(a, A1, A2, C, burn, T);
    [Yr1, Zr] = bvar.util.build_lags(Yr(3:end, :), p);
    [~, ~, s2r] = bvar.priors.minn(p, .2^2, .1^2, 10^2, Yr(1:4, :), Yr1, 4);
    VAr = @(c1) [c3; kron(1./(1:p)'.^2, c1./s2r)];
    lcr = fminbnd(@(lc) -bvar.ml.ncp(Yr1, Zr, VAr(exp(lc)), nu0, diag(s2r)), log(1e-3), log(10));
    [~, Bhat] = bvar.ml.ncp(Yr1, Zr, VAr(exp(lcr)), nu0, diag(s2r));
    [~, Bstd] = bvar.ml.ncp(Yr1, Zr, VAr(.04), nu0, diag(s2r));
    mc(r, :) = [exp(lcr), mae((Zr'*Zr)\(Zr'*Yr1)), mae(Bhat), mae(Bstd)];
end
fprintf('\n%d additional datasets: chosen c1 median %.3f, quartiles %.3f and %.3f\n', R, ...
    median(mc(:, 1)), quantile(mc(:, 1), .25), quantile(mc(:, 1), .75));
fprintf('  average error: least squares %.3f, chosen c1 %.3f, c1 = 0.04 %.3f\n', mean(mc(:, 2:4)));
fprintf('  chosen c1 more accurate than least squares in %d, than c1 = 0.04 in %d\n', ...
    sum(mc(:, 3) < mc(:, 2)), sum(mc(:, 3) < mc(:, 4)));

%% 5. the checks: the identity of the marginal likelihood
gap = zeros(size(c1grid));
for j = 1:numel(c1grid)
    [l, Aj, KA, Shat] = bvar.ml.ncp(Y, Z, VA(c1grid(j)), nu0, S0);
    pri = struct('A0', zeros(k, n), 'VA0', VA(c1grid(j)), 'nu0', nu0, 'S0', S0);
    est = struct('Ahat', Aj, 'Shat', Shat, 'KA', KA);
    gap(j) = l - bvar.ml.kron_bvar(Y, Z, pri, est);
end
[l, ~, KA, Shat] = bvar.ml.ncp(Y, Z, VA(c1hat), nu0, S0);
CS = chol(Sigtrue, 'lower');
E = (Y - Z*Atrue)/CS';
llike = -T*n/2*log(2*pi) - T*sum(log(diag(CS))) - .5*sum(E(:).^2);   % log p(Y | A, Sig)
lpri = bvar.ml.lniwpdf(Atrue, Sigtrue, zeros(k, n), sparse(diag(1./VA(c1hat))), nu0, S0);
lpost = bvar.ml.lniwpdf(Atrue, Sigtrue, Ahat, KA, nu0 + T, Shat);
fprintf('\nclosed form against bvar.ml.kron_bvar, over the grid: largest gap %.1e\n', ...
    max(abs(gap)));
fprintf('at the true values, log p(Y | A, Sig) + log p(A, Sig) - log p(A, Sig | Y) = %.6f\n', ...
    llike + lpri - lpost);
fprintf('  against the closed form %.6f: gap %.1e\n', l, llike + lpri - lpost - l);
assert(max(abs(gap)) < 1e-8 && abs(llike + lpri - lpost - l) < 1e-8, ...
    'the closed form disagrees with the identity of the marginal likelihood');

%% 6. the lag length: the adaptation
pmax = 4; lmlp = zeros(pmax, 1); c1p = zeros(pmax, 1);
for q = 1:pmax
    [Yq, Zq] = bvar.util.build_lags(Yfull(5-q:end, :), q);   % the same T rows of Y for every q
    assert(isequal(Yq, Y), 'the lag lengths explain different observations');
    VAq = @(c1) [c3; kron(1./(1:q)'.^2, c1./s2)];
    [lcq, f] = fminbnd(@(lc) -bvar.ml.ncp(Yq, Zq, VAq(exp(lc)), nu0, S0), log(1e-3), log(10));
    c1p(q) = exp(lcq); lmlp(q) = -f;
end
fprintf('\nlag length: chosen c1 %s, log ML %s\n', mat2str(c1p', 3), mat2str(lmlp', 5));

%% 7. the asymmetric conjugate prior of tutorials/shrinkage: own and other lags apart
[~, kappa] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [.04 .01], 'redu');
[~, kappa_b] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [.2 .2], 'redu');   % another start
fprintf('\nasymmetric conjugate prior: own lags %.4f, other lags %.4f\n', kappa(1:2));
fprintf('  from another starting value: %.4f, %.4f\n', kappa_b(1:2));

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 7]);
subplot(1, 2, 1);
semilogx(c1grid, lmlgrid, 'k'); hold on; box off
xline(c1hat, 'k--'); xline(.04, 'k:'); xlim(c1grid([1 end]));
xlabel('c_1'); title('Log marginal likelihood', 'FontWeight', 'normal');
subplot(1, 2, 2);
own = false(k, n); own(2:n+1, :) = logical(eye(n)); own(n+2:2*n+1, :) = logical(eye(n));
lagrows = repmat((1:k)' > 1, 1, n);
paths = reshape(pathA, k*n, []);
semilogx(c1grid, paths(lagrows(:) & ~own(:), :), 'Color', [.6 .6 .6]); hold on; box off
semilogx(c1grid, paths(own(:), :), 'k');
xline(c1hat, 'k--'); xline(.04, 'k:'); xlim(c1grid([1 end]));
xlabel('c_1'); title('Posterior means', 'FontWeight', 'normal');
exportgraphics(fig, fullfile(here, 'fig_marglik.png'), 'Resolution', 150); close(fig);

function Yfull = simulate_var(a, A1, A2, C, burn, T)
% the data of section 1 from the current state of the generator, with the same presample rows
n = numel(a);
Yfull = zeros(burn + 2 + T, n);
for t = 3:burn+2+T
    Yfull(t,:) = a' + Yfull(t-1,:)*A1' + Yfull(t-2,:)*A2' + (C*randn(n,1))';
end
Yfull = Yfull(burn-1:end, :);
end
