%% guide - the code of guides/gibbs_sampler/README.md: a Gibbs sampler for the VAR of
% guides/var_stacking under the independent normal and inverse-Wishart prior. The section called
% "the code to adapt" contains the sampler. Saves fig_chain.png next to this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
rng(1, 'twister');

%% 1. the data and the prior of guides/coefficient_draw
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
[~, V, s2] = bvar.priors.minn(p, .2^2, .1^2, 10^2, Y0, Y, 4);

%% 2. the Gibbs sampler: the code to adapt
nu0 = n + 2; S0 = diag(s2);        % Sig ~ IW(nu0, S0), with prior mean S0
iV = spdiags(1./V, 0, k*n, k*n);   % beta ~ N(0, V), the Minnesota variances
nsim = 20000; burnin = 1000;
Sig = S0;                          % the starting value
store_A = zeros(k, n, nsim); store_Sig = zeros(n, n, nsim);
for isim = 1:burnin + nsim
    A = bvar.samplers.var_coef_joint(Y, Z, Sig, ones(T, 1), iV);   % beta given Sig
    U = Y - Z*A;
    Sig = iwishrnd(S0 + U'*U, nu0 + T);                            % Sig given beta
    if isim > burnin
        store_A(:, :, isim - burnin) = A;
        store_Sig(:, :, isim - burnin) = Sig;
    end
end

%% 3. the checks of the chain
low = tril(true(n));
draws = [reshape(store_A, k*n, nsim)', reshape(store_Sig(repmat(low, 1, 1, nsim)), [], nsim)'];
truth = [Atrue(:); Sigtrue(low)];
L = 100;                           % Bartlett weights up to lag 100
IF = bvar.diag.inefficiency_factor(draws, L);
se = bvar.diag.mcse(draws, L);
[Zg, pg] = bvar.diag.geweke(draws);
q = quantile(draws, [.05 .95]);
fprintf('\n%d draws after %d burn-in, %d parameters\n', nsim, burnin, size(draws, 2));
fprintf('inefficiency factors: coefficients %.2f to %.2f, Sig %.2f to %.2f\n', ...
    min(IF(1:k*n)), max(IF(1:k*n)), min(IF(k*n+1:end)), max(IF(k*n+1:end)));
[zmax, jmax] = max(abs(Zg));
fprintf('Geweke Z: largest |Z| %.2f, for parameter %d of 27; %d of %d beyond 1.96\n', zmax, jmax, ...
    sum(abs(Zg) > 1.96), numel(Zg));
mz = zeros(1, 10);                 % the same chain from ten other seeds
for s = 1:10
    rng(s + 1, 'twister');
    mz(s) = max(abs(bvar.diag.geweke(gibbs_draws(Y, Z, iV, S0, nu0, nsim, burnin))));
end
fprintf('the largest |Z| over ten other seeds: %.2f to %.2f\n', min(mz), max(mz));
fprintf('true values inside the 90%% intervals: %d of %d\n', ...
    sum(truth' > q(1,:) & truth' < q(2,:)), numel(truth));
names = {'A(2,1), own lag 1 of eq 1', 'A(2,2), lag 1 of y1 in eq 2', 'Sig(1,1)', 'Sig(2,1)'};
cols = [2, k + 2, k*n + 1, k*n + 2];
fprintf('%-30s %6s %8s %18s %6s %8s\n', '', 'true', 'mean', '90% interval', 'IF', 'MCSE');
for j = 1:numel(cols)
    c = cols(j);
    fprintf('%-30s %6.2f %8.3f   (%6.3f, %6.3f) %6.2f %8.4f\n', names{j}, truth(c), ...
        mean(draws(:, c)), q(1, c), q(2, c), IF(c), se(c));
end

%% 4. the natural conjugate prior: the checks
VA = [100; kron(1./(1:p)'.^2, .04./s2)];   % prior variances of the rows of A, times Sig
iVA = 1./VA;
Sig = S0;
ncp_A = zeros(k, n, nsim); ncp_Sig = zeros(n, n, nsim);
for isim = 1:burnin + nsim
    iVb = kron(sparse(inv(Sig)), spdiags(iVA, 0, k, k));   % prior precision, inv(Sig) kron inv(VA)
    A = bvar.samplers.var_coef_joint(Y, Z, Sig, ones(T, 1), iVb);
    U = Y - Z*A;
    Sig = iwishrnd(S0 + U'*U + A'*(A.*iVA), nu0 + T + k);
    if isim > burnin
        ncp_A(:, :, isim - burnin) = A;
        ncp_Sig(:, :, isim - burnin) = Sig;
    end
end
KA = diag(iVA) + Z'*Z;             % the closed-form posterior
Ahat = KA\(Z'*Y);
Shat = S0 + Y'*Y - Ahat'*KA*Ahat;
ESig = Shat/(nu0 + T - n - 1);
sdA = sqrt(kron(diag(ESig), diag(inv(KA))));
ncp = [reshape(ncp_A, k*n, nsim)', reshape(ncp_Sig(repmat(low, 1, 1, nsim)), [], nsim)'];
exact_mean = [Ahat(:); ESig(low)];
se_ncp = bvar.diag.mcse(ncp, L);
fprintf('\nnatural conjugate prior, Gibbs against the closed form:\n');
zgap = max(abs(mean(ncp)' - exact_mean)./se_ncp');
sdgap = max(abs(std(ncp(:, 1:k*n))'./sdA - 1));
fprintf('  posterior means: largest gap %.2f Monte Carlo standard errors\n', zgap);
fprintf('  posterior sd of the coefficients: largest relative gap %.3f\n', sdgap);
assert(zgap < 4.5 && sdgap < .05, 'the Gibbs sampler disagrees with the closed-form posterior');

%% 5. a forecast of y_{T+1} from the draws: the adaptation
xT1 = [1, Y(end, :), Y(end-1, :)];   % x_{T+1}' = (1, y_T', y_{T-1}')
yT1 = zeros(nsim, n);
for r = 1:nsim
    yT1(r, :) = xT1*store_A(:, :, r) + (chol(store_Sig(:, :, r), 'lower')*randn(n, 1))';
end
qf = quantile(yT1, [.05 .95]);
condmean = zeros(nsim, n); dSig = zeros(nsim, n);
for r = 1:nsim
    condmean(r, :) = xT1*store_A(:, :, r);
    dSig(r, :) = diag(store_Sig(:, :, r))';
end
fprintf('\nthe forecast of y_{T+1}: mean %s, 90%% intervals %s to %s\n', ...
    mat2str(mean(yT1), 3), mat2str(qf(1,:), 3), mat2str(qf(2,:), 3));
fprintf('  mean against the average of the conditional means x''A: %.1e (largest)\n', ...
    max(abs(mean(yT1) - mean(condmean))));
fprintf('  variance against mean(diag(Sig)) + var(x''A): largest relative gap %.3f\n', ...
    max(abs(var(yT1)./(mean(dSig) + var(condmean)) - 1)));

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 8]);
labels = {'A_{21}: own lag 1, eq. 1', 'A_{22}: lag 1 of y_1, eq. 2', '\Sigma_{11}', '\Sigma_{21}'};
for j = 1:numel(cols)
    c = cols(j);
    subplot(2, 4, j);
    plot(draws(:, c), 'k'); box off; xlim([1 nsim]);
    title(labels{j}, 'FontWeight', 'normal');
    subplot(2, 4, 4 + j);
    histogram(draws(:, c), 50, 'Normalization', 'pdf', 'FaceColor', [.75 .75 .75], ...
        'EdgeColor', 'none');
    hold on; xline(truth(c), 'k--'); box off
end
exportgraphics(fig, fullfile(here, 'fig_chain.png'), 'Resolution', 150); close(fig);

function draws = gibbs_draws(Y, Z, iV, S0, nu0, nsim, burnin)
% the loop of Step 3, returning the draws of vec(A) and of the lower triangle of Sig, one row each
[T, n] = size(Y); k = size(Z, 2); low = tril(true(n));
Sig = S0;
draws = zeros(nsim, k*n + n*(n+1)/2);
for isim = 1:burnin + nsim
    A = bvar.samplers.var_coef_joint(Y, Z, Sig, ones(T, 1), iV);
    U = Y - Z*A;
    Sig = iwishrnd(S0 + U'*U, nu0 + T);
    if isim > burnin
        draws(isim - burnin, :) = [A(:); Sig(low)]';
    end
end
end
