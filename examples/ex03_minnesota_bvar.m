%% ex03 - A small BVAR with a Minnesota / natural-conjugate prior
%
% THE MODEL. A VAR(p) with an intercept, written a row at a time:
%
%       y_t' = z_t' A + eps_t',      eps_t ~ N(0, Sig),
%       z_t  = [1, y_{t-1}', ..., y_{t-p}']'
%
% so A is the k x n coefficient matrix (k = 1 + n*p) and Sig the n x n error
% covariance. With n = 20 variables and p = 4 lags there are k*n = 1620
% coefficients for ~200 observations, which is why a shrinkage prior is not
% optional. This example uses n = 3 so you can read the numbers.
%
% THE PRIOR. Two constructors in core/+bvt/+priors/, and it is worth being
% clear about how they differ:
%
%   bvt.priors.minn  - the classic Minnesota prior. Each coefficient gets its
%       OWN prior variance: c1/l^2 on own lag l, c2*sig2_i/(l^2*sig2_j) on the
%       lag of variable j in equation i, c3 on the intercept. Because the
%       variance depends on the equation i, this prior does NOT factor as a
%       Kronecker product, so Sig must be fixed (or drawn separately) and the
%       posterior is not available in closed form.
%
%   bvt.priors.niw   - the natural-conjugate (normal-inverse-Wishart) prior:
%       A | Sig ~ MN(A0, diag(VA0), Sig),  Sig ~ IW(nu0, S0). The prior
%       variance c1/(l^2*sig2_j) drops the equation index i, which is exactly
%       the restriction that buys the Kronecker structure - and with it an
%       ANALYTIC posterior from which we draw directly, no MCMC at all.
%
% Both scale their hyperparameters by sig2, the residual variances of
% univariate AR(4) fits (bvt.priors.resid_var_ar4). That is what makes a
% single scalar c1 mean the same thing for an interest rate and for GDP growth.
%
% DATA. Read-only from replications/chan2020_jbes_kronecker/legacy/data_Q.csv,
% the quarterly US macro panel of Chan (2020, JBES), 1959Q1-2013Q4. We take
% the first three columns of that file and hold out the LAST observation to
% score a one-step-ahead forecast against.
%
% WHAT TO LOOK AT: the prior-variance comparison in section 3 (how much the
% Kronecker restriction costs you in flexibility), the posterior coefficient
% table, and the one-step forecast intervals against the realized values.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

rng(20260903, 'twister')
repo = fileparts(fileparts(mfilename('fullpath')));
fprintf('\n=== ex03: a small BVAR with a Minnesota-type prior ===\n');

%% ------------------------------------------------------------------
%  1. Data (read-only from the legacy replication folder)
%  ------------------------------------------------------------------
csv  = fullfile(repo, 'replications', 'chan2020_jbes_kronecker', 'legacy', 'data_Q.csv');
raw  = load(csv);                    % never written to; the legacy folders are frozen
data = raw(:, 1:3);                  % first three of the paper's 20 series
p    = 4;                            % lag length (the replication's setting)

    % hold out the final quarter, so we have something to score against
data_est = data(1:end-1, :);
y_real   = data(end, :);             % the realized values we will forecast

Y0     = data_est(1:4, :);           % initial conditions (4 rows: the AR(4)
shortY = data_est(5:end, :);         %   prior fits need exactly four)
[T, n] = size(shortY);
k      = 1 + n*p;

fprintf('\nsample: T = %d quarters used for estimation, n = %d variables, p = %d lags\n', T, n, p);
fprintf('        k = 1 + n*p = %d coefficients per equation, k*n = %d in total\n', k, k*n);
fprintf('        series means %s, sds %s\n', mat2str(round(mean(shortY),2)), ...
    mat2str(round(std(shortY),2)));

%% ------------------------------------------------------------------
%  2. Design matrix: bvt.util.build_lags
%     Convention used everywhere in this toolkit: intercept FIRST, then the
%     complete lag-1 block, then the lag-2 block, ...
%  ------------------------------------------------------------------
[Y, X] = bvt.util.build_lags([Y0(end-p+1:end,:); shortY], p);
fprintf('\nbvt.util.build_lags -> Y is %dx%d, X is %dx%d; X(1,:) = %s\n', ...
    size(Y,1), size(Y,2), size(X,1), size(X,2), mat2str(round(X(1,1:4),2)));

%% ------------------------------------------------------------------
%  3. Priors
%  ------------------------------------------------------------------
sig2 = bvt.priors.resid_var_ar4(Y0, shortY);       % AR(4) residual variances
fprintf('\nAR(4) residual variances sig2 (the prior''s units of measurement):\n');
fprintf('  %s\n', mat2str(round(sig2', 3)));

c1 = 0.2^2;    % overall shrinkage on the lag coefficients
c2 = 100;      % (near-)flat prior on the intercepts
c3 = 100;      % same, for the bvt.priors.minn call below

    % (a) natural-conjugate NIW prior - the one we actually use to sample
[A0, VA0, nu0, S0] = bvt.priors.niw(p, [c1 c2], Y0, shortY, 'largebvar_nc');

    % (b) the classic Minnesota prior, for comparison only
[beta_Minn, V_Minn] = bvt.priors.minn(p, c1, c1, c3, Y0, shortY, 4);

fprintf('\nprior standard deviations on the three LAG-1 coefficients, by equation:\n');
fprintf('  %-30s %10s %10s %10s\n', 'prior / equation', 'on y1(t-1)', 'on y2(t-1)', 'on y3(t-1)');
for i = 1:n
    fprintf('  %-30s %10.4f %10.4f %10.4f\n', ...
        sprintf('natural conjugate VA0, eq %d', i), sqrt(VA0(2)), sqrt(VA0(3)), sqrt(VA0(4)));
end
for i = 1:n
    b = (i-1)*k;
    fprintf('  %-30s %10.4f %10.4f %10.4f\n', ...
        sprintf('Minnesota V_Minn,      eq %d', i), ...
        sqrt(V_Minn(b+2)), sqrt(V_Minn(b+3)), sqrt(V_Minn(b+4)));
end
fprintf('  -> VA0 repeats: the natural-conjugate prior is the SAME in every\n');
fprintf('     equation, which is exactly the Kronecker restriction. V_Minn\n');
fprintf('     rescales row i by sig2_i/sig2_j and treats the own lag specially.\n');
fprintf('     That extra freedom is what costs you the analytic posterior.\n');
fprintf('     (Read the two rows in their own units: VA0 is the variance of A\n');
fprintf('     GIVEN Sig, so the prior variance of the (i,j) element is really\n');
fprintf('     VA0(i)*Sig(j,j) - the equation-to-equation scaling that V_Minn\n');
fprintf('     writes down explicitly enters the NIW prior through Sig instead.\n');
fprintf('     What the NIW prior cannot do is treat own and cross lags\n');
fprintf('     differently, which is why c1 is used for both here.)\n');
fprintf('\ninverse-Wishart prior on Sig: nu0 = %d, S0 = diag(%s)\n', ...
    nu0, mat2str(round(diag(S0)', 3)));

%% ------------------------------------------------------------------
%  4. Posterior. With the NIW prior it is available in closed form:
%       A | Sig, y ~ MN(Ahat, KA^-1, Sig),   Sig | y ~ IW(nu0 + T, Shat)
%     so we draw Sig first and then A given Sig - INDEPENDENT draws, no
%     Markov chain, no burn-in, no convergence to worry about.
%     (These are the same six lines as replications/chan2020_jbes_kronecker/
%      run_all.m, subfunction post_bvar.)
%  ------------------------------------------------------------------
XX     = X'*X;
KA     = sparse(1:k,1:k,1./VA0) + XX;
Ahat   = KA\(sparse(1:k,1:k,VA0)\A0 + X'*Y);
Shat   = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Y'*Y - Ahat'*KA*Ahat;
Shat   = (Shat + Shat')/2;                       % symmetrize against rounding

nsim = 2000;
store_A   = zeros(nsim, k*n);
store_Sig = zeros(n, n);
store_fc  = zeros(nsim, n);
z_next    = [1 shortY(end,:) shortY(end-1,:) shortY(end-2,:) shortY(end-3,:)];
CKA       = chol(KA, 'lower');

for isim = 1:nsim
    Sig  = iwishrnd(Shat, nu0 + T);              % Sig | y
    CSig = chol(Sig, 'lower');
    A    = Ahat + (CKA'\randn(k,n))*CSig';       % A | Sig, y  (matrix normal)

    store_A(isim,:)  = A(:)';
    store_Sig        = store_Sig + Sig;

        % one-step-ahead predictive draw: y_{T+1} = z_{T+1}'A + eps
    store_fc(isim,:) = z_next*A + (CSig*randn(n,1))';
end
A_mean   = reshape(mean(store_A)', k, n);
A_sd     = reshape(std(store_A)',  k, n);
Sig_mean = store_Sig/nsim;

%% ------------------------------------------------------------------
%  5. Read the posterior
%  ------------------------------------------------------------------
rownames = cell(k,1);
rownames{1} = 'intercept';
for l = 1:p
    for j = 1:n
        rownames{1+(l-1)*n+j} = sprintf('y%d(t-%d)', j, l);
    end
end

fprintf('\nposterior mean (posterior sd) of A, first five rows:\n');
fprintf('  %-12s', '');
for j = 1:n, fprintf('%18s', sprintf('eq %d', j)); end
fprintf('\n');
for i = 1:5
    fprintf('  %-12s', rownames{i});
    for j = 1:n
        fprintf('%10.3f (%.3f)', A_mean(i,j), A_sd(i,j));
    end
    fprintf('\n');
end

    % is the posterior-mean VAR stable? build the companion matrix
F     = [A_mean(2:end,:)'; eye(n*(p-1)) zeros(n*(p-1), n)];
maxev = max(abs(eig(F)));
if maxev < 1, verdict = 'stable'; else, verdict = 'explosive'; end
fprintf('\nlargest companion-matrix eigenvalue at A_mean: %.3f (%s)\n', maxev, verdict);

fprintf('\nposterior mean of Sig (error covariance):\n');
disp(round(Sig_mean, 3));
fprintf('implied error correlations:\n');
disp(round(corrcov(Sig_mean), 3));

%% ------------------------------------------------------------------
%  6. The one-step-ahead forecast, scored against the held-out quarter
%  ------------------------------------------------------------------
fc_mean = mean(store_fc);
fc_ci   = quantile(store_fc, [.05 .95]);
fprintf('\none-step-ahead forecast of the held-out quarter:\n');
fprintf('  %-8s %10s %10s %22s %10s\n', 'variable', 'realized', 'pred mean', '90% interval', 'covered');
for j = 1:n
    covered = y_real(j) >= fc_ci(1,j) && y_real(j) <= fc_ci(2,j);
    fprintf('  %-8s %10.3f %10.3f   [%8.3f, %8.3f] %10s\n', sprintf('y%d', j), ...
        y_real(j), fc_mean(j), fc_ci(1,j), fc_ci(2,j), string(covered));
end
fprintf('  (the intervals are wide because they include BOTH coefficient\n');
fprintf('   uncertainty and the one-period error eps_{T+1}.)\n');

%% ------------------------------------------------------------------
%  7. Figure (non-blocking)
%  ------------------------------------------------------------------
figure('Name','ex03 Minnesota BVAR');
for j = 1:n
    subplot(n,1,j)
    tail = 24;                       % last six years of the estimation sample
    plot(1:tail, shortY(end-tail+1:end, j), 'k-o', 'MarkerSize', 3); hold on
    plot(tail+1, y_real(j), 'ks', 'MarkerFaceColor', 'k');
    plot(tail+1, fc_mean(j), 'rd', 'MarkerFaceColor', 'r');
    plot([tail+1 tail+1], fc_ci(:,j), 'r-', 'LineWidth', 1.5);
    hold off; box off
    title(sprintf('variable %d: data, realized (black square), forecast + 90%% interval (red)', j))
end
drawnow

fprintf('\nex03 done. Next: ex04_bvar_sv_blocks (the same VAR, now with SV).\n');
