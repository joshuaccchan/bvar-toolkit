%% guide - the code of guides/coefficient_draw/README.md: the posterior mean and a draw of the
% coefficients of the VAR of guides/var_stacking under the Minnesota prior. The sections called
% "the code to adapt" contain everything a sampler uses. Saves fig_shrinkage.png next to this
% file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
rng(1, 'twister');
rel = @(a, b) max(abs(a(:) - b(:)))/max(abs(b(:)));

%% 1. the data: the VAR(2) of guides/var_stacking, with four presample rows
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

%% 2. the prior: the code to adapt
c1 = .2^2; c2 = .1^2; c3 = 10^2;   % own lags, other lags, intercepts
[beta0, V, s2] = bvar.priors.minn(p, c1, c2, c3, Y0, Y, 4);
Sighat = diag(s2);                 % Sig replaced by the AR(4) residual variances

%% 2. the checks
Vmat = reshape(V, k, n);           % column i holds the prior variances of equation i
own = false(k, n); other = false(k, n); lag = zeros(k, 1);
for r = 2:k
    lag(r) = ceil((r - 1)/n);
    j = mod(r - 2, n) + 1;         % the variable whose lag row r holds
    own(r, j) = true; other(r, setdiff(1:n, j)) = true;
end
fprintf('\nthe AR(4) residual variances s2: %s\n', mat2str(s2', 4));
fprintf('prior standard deviations, equation 1: intercept, lags 1 and 2 of variables 1-3\n');
disp(sqrt(Vmat(:, 1))');

%% 3. the posterior mean and a draw: the code to adapt
iV = spdiags(1./V, 0, k*n, k*n);   % the prior precision, diagonal
[Adraw, Ahat] = bvar.samplers.var_coef_joint(Y, Z, Sighat, ones(T, 1), iV);

%% 3. the checks
X = bvar.util.surform2(Z, n);
y = reshape(Y', T*n, 1);
W = kron(speye(T), inv(Sighat));
K = full(diag(1./V) + X'*W*X);
bhat = K\(X'*W*y);
gapM = rel(Ahat(:), bhat);
fprintf('\nthe posterior mean against dense algebra on the stacked form: %.1e (relative)\n', gapM);
R = 10000;
B = zeros(k*n, R);
for r = 1:R
    B(:, r) = reshape(bvar.samplers.var_coef_joint(Y, Z, Sighat, ones(T, 1), iV), [], 1);
end
Kinv = inv(K);
sd = sqrt(diag(Kinv));
dmean = max(abs(mean(B, 2) - bhat)./sd);
dsd = max(abs(std(B, 0, 2)./sd - 1));
dcorr = max(abs(corrcoef(B') - Kinv./(sd*sd')), [], 'all');
fprintf('%d draws: mean against the posterior mean, largest gap %.2f posterior sd\n', R, dmean);
fprintf('%d draws: sd against the posterior sd, largest relative gap %.3f\n', R, dsd);
fprintf('%d draws: correlations against the posterior correlations, largest gap %.3f\n', R, dcorr);
Aeq = zeros(k, n);                 % one regression per equation, since Sighat is diagonal
for i = 1:n
    Ki = diag(1./Vmat(:, i)) + Z'*Z/s2(i);
    Aeq(:, i) = Ki\(Z'*Y(:, i)/s2(i));
end
gapE = rel(Aeq, Ahat);
fprintf('equation by equation against the joint posterior mean: %.1e (relative)\n', gapE);
Ols = (Z'*Z)\(Z'*Y);
[~, Aloose] = bvar.samplers.var_coef_joint(Y, Z, Sighat, ones(T, 1), 1e-10*speye(k*n));
gapL = rel(Aloose, Ols);
fprintf('a nearly flat prior against least squares: %.1e (relative)\n', gapL);
Wtrue = kron(speye(T), inv(Sig));
bhat_true = (diag(1./V) + X'*Wtrue*X)\(X'*Wtrue*y);
[~, Ahat_true] = bvar.samplers.var_coef_joint(Y, Z, Sig, ones(T, 1), iV);
gapT = rel(Ahat_true(:), bhat_true);
fprintf('with the true, full Sig, against dense algebra: %.1e (relative)\n', gapT);
assert(gapM < 1e-10 && gapE < 1e-10 && gapL < 1e-8 && gapT < 1e-10, ...
    'the posterior mean disagrees with dense algebra');
assert(dmean < .1 && dsd < .05 && dcorr < .1, 'the draws do not match the posterior');

%% 4. the posterior mean against least squares
fprintf('\nsum of absolute values, posterior mean over least squares:\n');
groups = {own & lag == 1, 'own lag 1'; own & lag == 2, 'own lag 2'; ...
    other & lag == 1, 'other lag 1'; other & lag == 2, 'other lag 2'; ...
    repmat((1:k)' == 1, 1, n), 'intercepts'};
for g = 1:size(groups, 1)
    m = groups{g, 1};
    fprintf('  %-12s %.2f\n', groups{g, 2}, sum(abs(Ahat(m)))/sum(abs(Ols(m))));
end
fprintf('first own lags: least squares %s, posterior mean %s\n', ...
    mat2str(diag(Ols(2:n+1, :))', 3), mat2str(diag(Ahat(2:n+1, :))', 3));
[~, Vfree] = bvar.priors.minn(p, c1, 100, c3, Y0, Y, 4);   % c2 = 100: other lags almost unshrunk
[~, Afree] = bvar.samplers.var_coef_joint(Y, Z, Sighat, ones(T, 1), ...
    spdiags(1./Vfree, 0, k*n, k*n));
Rz = corrcoef(Z(:, 2:end));
fprintf('first own lags with c2 = 100: %s; correlation of the first lags of y1 and y2: %.2f\n', ...
    mat2str(diag(Afree(2:n+1, :))', 3), Rz(1, 2));
fprintf('estimation error, mean absolute: least squares %.3f, posterior mean %.3f\n', ...
    mean(abs(Ols(2:end,:) - A(2:end,:)), 'all'), mean(abs(Ahat(2:end,:) - A(2:end,:)), 'all'));

%% 5. a prior mean of one on the first own lags: the adaptation
A0 = zeros(k, n); A0(2:n+1, :) = eye(n);   % one on each first own lag
[G, Ghat] = bvar.samplers.var_coef_joint(Y - Z*A0, Z, Sighat, ones(T, 1), iV);
Adraw0 = A0 + G;
Ahat0 = A0 + Ghat;
bhat0 = K\(diag(1./V)*A0(:) + X'*W*y);
gap0 = rel(Ahat0(:), bhat0);
fprintf('\nprior mean of one on the first own lags, against dense algebra: %.1e (relative)\n', ...
    gap0);
assert(gap0 < 1e-10, 'the draw with a nonzero prior mean disagrees with dense algebra');
fprintf('posterior means of the first own lags: %s (prior mean zero), %s (prior mean one)\n', ...
    mat2str(diag(Ahat(2:n+1, :))', 3), mat2str(diag(Ahat0(2:n+1, :))', 3));

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 10 9]);
hold on
lims = [-.2 .6];
plot(lims, lims, ':', 'Color', [.5 .5 .5], 'HandleVisibility', 'off');
mk = {'o', 's', '^', 'v'};
for g = 1:4
    m = groups{g, 1};
    plot(Ols(m), Ahat(m), mk{g}, 'Color', 'k', 'MarkerSize', 6, 'DisplayName', groups{g, 2});
end
xlim(lims); ylim(lims); axis square; box off
xlabel('least squares'); ylabel('posterior mean');
legend('Location', 'northwest'); legend boxoff
exportgraphics(fig, fullfile(here, 'fig_shrinkage.png'), 'Resolution', 150); close(fig);
