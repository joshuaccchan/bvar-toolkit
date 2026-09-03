%% ex04 - A REDUCED-FORM BVAR with stochastic volatility, equation by equation
%
% This example teaches the sampler of
%
%       replications/chan2023_joe_mlvarsv/legacy/VAR_ARSV_redu.m
%
% (Chan, 2023, JoE 235(2): 1419-1446, "Comparing Stochastic Volatility
% Specifications for Large Bayesian VARs"), on simulated data small enough to
% check every number against the truth.
%
% THE MODEL. Unlike ex03's structural triangular form, the VAR here is written
% in REDUCED form - each equation's regressors are lags only, no contemporaneous
% variables:
%
%       Y = X*A + E,        A is k x n,  k = 1 + n*p,  intercept first,
%       B0*eps_t = u_t,     u_{it} ~ N(0, exp(h_{it})),
%       h_{it} = mu_i + phi_i*(h_{i,t-1} - mu_i) + v_{it},  v_{it} ~ N(0, sig2_i),
%       h_{i1} ~ N(mu_i, sig2_i/(1 - phi_i^2)).
%
% E is T x n with rows eps_t'. B0 is LOWER UNITRIANGULAR: ones on the diagonal,
% n*(n-1)/2 free elements below it. So the reduced-form errors are correlated -
% Var(eps_t) = B0^{-1} diag(exp(h_t)) B0^{-1}' - and it is the ORTHOGONALIZED
% errors B0*eps_t that carry the n independent AR(1) stochastic volatilities.
%
% WHY THIS IS HARDER THAN ex03, AND WHY IT IS WORTH IT. In the structural form
% of ex03/MAHP, conditioning on the log-volatilities makes the n equations
% independent, so each equation is a self-contained weighted regression
% (bvt.samplers.eq_gauss). Here A is a reduced-form object: change equation ii's
% coefficients - column ii of A, see the layout note below - and you move the
% orthogonalized errors of SEVERAL equations at once, because the
% orthogonalization mixes them. The naive conditional for the whole of A is a
% k*n x k*n system - 21 x 21 here, but 3660 x 3660 at the paper's n = 15, p = 4,
% and a dense Cholesky of that every sweep is what you are trying to avoid.
%
% THE ALGORITHM. Draw A EQUATION BY EQUATION (section 5, block 1), so the
% expensive step becomes n Cholesky factorizations of size k x k. The trick that
% makes each of those systems small is triangularity - see the long comment in
% block 1.
%
% BE PRECISE ABOUT THE LAYOUT. A is k x n with beta = vec(A) - the convention of
% the legacy code and of the rest of this toolkit - so equation ii's k
% coefficients are COLUMN ii of A:  Y(:,ii) = X*A(:,ii) + error. The loop below
% therefore draws A(:,ii), which is ONE EQUATION per iteration. Store the same
% coefficients the other way round, as the n x k matrix A', and the identical
% draw would be row by row. Column versus row here is a layout convention, not a
% difference in the algorithm; a column of A is an equation, not a lag and not
% one variable across equations. (Do not switch the layout: keeping A as k x n
% is what makes this script diff line by line against the legacy file.) The four
% blocks:
%
%   1. A equation by equation       (the centrepiece; written out inline here)
%   2. the free elements of B0      bvt.samplers.alp_tri_cs
%   3. the n log-volatility paths   bvt.sv.ksc_ar1_mean
%   4. (mu, phi, sig2) per equation bvt.sv.sv_params
%
% CONTRAST WITH bvt.samplers.eq_svar_oi. That function is the same idea for the
% order-invariant SVAR-SV of Chan, Koop and Yu (2024, JBES): same standardized
% stacking, same per-equation Cholesky. The one difference is decisive - its B0
% is a GENERAL rotation, not triangular, so every equation's residual depends on
% every column of A and its Wi = kron(B0(:,ii),X) keeps all n*T rows. Here
% B0(1:ii-1,ii) = 0, so rows 1..ii-1 drop out exactly. Read its header next to
% block 1 below: the two blocks differ in exactly one index range.
%
% TWO HONEST CAVEATS. This script ILLUSTRATES VAR_ARSV_redu; it does not
% reproduce it bitwise, for two reasons.
%   (a) bvt.sv.sv_params canonicalizes the OISV copy of sample_SVpara, NOT the
%       ml_varsv copy that VAR_ARSV_redu calls. The two differ - phi truncation
%       bound .999 vs .998, and the handling of the h columns - and are on the
%       never-merge list in tests/variant_map.md. Same conditional, different
%       floating-point path and different rng consumption.
%   (b) The shrinkage hyperparameters kappa are held FIXED at the paper's preset
%       values. VAR_ARSV_redu draws kappa1, kappa2 and kappa4 from
%       generalized-inverse-Gaussian conditionals with gigrnd every sweep, which
%       is why its prior variances are rebuilt inside the loop. Fixing them costs
%       one Gibbs block and buys a much shorter script. (ex04's predecessor, and
%       replications/chan2021_ijf_mahp/run_all.m, show that block switched on.)
%
% A NAMING TRAP, worth two lines because it has caught people. In THIS paper's
% code `alp` is vec(A), the VAR coefficients, and `beta` collects the free
% elements of B0 - the opposite of the MAHP convention used in ex03 and in
% bvt.priors.vtheta, where Vbeta holds the VAR coefficients and Valp the impact
% matrix. The variable names below follow the ml_varsv legacy file, so the code
% diffs line by line against it. Read the comments, not the letters.
%
% WHAT TO LOOK AT: the row-count line printed in block 1 (the truncation, in
% numbers), the coefficient RMSE against ordinary least squares, the recovered
% B0 elements, and the correlation between the estimated and true log-volatility
% paths.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

rng(20260903, 'twister')
fprintf('\n=== ex04: a reduced-form BVAR-SV, drawn equation by equation ===\n');

%% ------------------------------------------------------------------
%  1. Simulate a 3-variable reduced-form VAR(2) with AR(1) stochastic
%     volatility on the orthogonalized errors.
%
%     Data are on a "percent" scale (variances of order one), which keeps the
%     log(u^2 + 1e-4) offset of the SV step harmless - see ex02, section 5.
%  ------------------------------------------------------------------
n = 3; p = 2; T = 200; n0 = 4; nburn_sim = 100;
k = 1 + n*p;                                % regressors per equation
k_alp  = k*n;                               % elements of A
k_beta = n*(n-1)/2;                         % free elements of B0

    % true reduced-form coefficient matrix A (k x n; COLUMN ii = equation ii)
B1 = [ .5  .1  .0;                          % B1(i,j) multiplies y_j(t-1) in eq i
       .0  .6  .1;
       .1  .0  .5];
B2 = .1*eye(n);                             % lag-2 block
c  = [.2; .1; .0];                          % intercepts
A_true = [c'; B1'; B2'];                    % note the transposes: A is k x n

    % true impact matrix: lower unitriangular, free elements (2,1) (3,1) (3,2)
beta_true = [.4; -.3; .5];
B0_true = eye(n);
B0_id   = nonzeros(tril(reshape(1:n^2,n,n),-1)');   % legacy line 32, verbatim
B0_true(B0_id) = beta_true;                 % -> [1 0 0; .4 1 0; -.3 .5 1]

    % true SV parameters, one AR(1) per orthogonalized error
mu_true   = log([1.0; 0.8; 1.2]);           % mean log variance
phi_true  = [.95; .98; .90];                % persistence
sig2_true = [.10; .05; .15];                % innovation variance of h

Ttot = T + n0 + nburn_sim;
h_full = zeros(Ttot,n);
h_full(1,:) = mu_true' + (sqrt(sig2_true./(1-phi_true.^2)).*randn(n,1))';
for t = 2:Ttot
    h_full(t,:) = mu_true' + phi_true'.*(h_full(t-1,:) - mu_true') ...
        + (sqrt(sig2_true).*randn(n,1))';
end
U_full = exp(h_full/2).*randn(Ttot,n);      % the ORTHOGONALIZED errors u_t
E_full = U_full/B0_true';                   % eps_t = B0^{-1} u_t, one period per row

Yfull = zeros(Ttot,n);
for t = p+1:Ttot
    Yfull(t,:) = [1, Yfull(t-1,:), Yfull(t-2,:)]*A_true + E_full(t,:);
end
Yfull  = Yfull(nburn_sim+1:end, :);
h_true = h_full(nburn_sim+1:end, :);
Y0 = Yfull(1:n0, :);                        % initial conditions (4 rows: the
Y  = Yfull(n0+1:end, :);                    %   AR(4) prior fits need exactly 4)
h_true = h_true(n0+1:end, :);

    % companion matrix of the truth, so we know the data are stationary
F = [A_true(2:end,:)'; eye(n*(p-1)) zeros(n*(p-1), n)];
fprintf('\nsimulated n = %d, p = %d, T = %d\n', n, p, T);
fprintf('  series sds %s; largest companion eigenvalue of A_true %.3f\n', ...
    mat2str(round(std(Y),2)), max(abs(eig(F))));
fprintf('  true log-volatility paths range over [%.2f, %.2f]\n', ...
    min(h_true(:)), max(h_true(:)));

%% ------------------------------------------------------------------
%  2. Design matrix. Toolkit convention: intercept first, then the complete
%     lag-1 block, then the lag-2 block, ... (identical to main_varsv.m's
%     inline construction, lines 45-50).
%  ------------------------------------------------------------------
[~, X] = bvt.util.build_lags([Y0(end-p+1:end,:); Y], p);
fprintf('\nbvt.util.build_lags -> X is %dx%d, k = 1 + n*p = %d\n', size(X,1), size(X,2), k);

%% ------------------------------------------------------------------
%  3. Priors, at the paper's preset values (main_varsv.m lines 56-62, 86-91).
%
%     Two constructors, both extracted verbatim from this paper's legacy
%     utility folder:
%       bvt.priors.minn      - the Minnesota prior on vec(A). Called with
%                              n0pre = 4, it IS ml_varsv's prior_Minn.m. Its
%                              third output sig2_hat is the vector of univariate
%                              AR(4) residual variances (the same numbers
%                              bvt.priors.resid_var_ar4 returns), which is what
%                              makes one scalar kappa mean the same thing for an
%                              interest rate and for GDP growth.
%       bvt.priors.impact_B0 - the prior on the free elements of B0, scaled
%                              kappa4*sig2_i/sig2_j. (This is ml_varsv's
%                              prior_B0.m. bvt.priors.vtheta builds the SAME
%                              object under the MAHP parameterization, as its
%                              Valp output - see the naming trap in the header.)
%  ------------------------------------------------------------------
kappa1 = .2^2;          % own lags
kappa2 = (.2^2)^2;      % other lags - an order of magnitude tighter
kappa3 = 100;           % intercepts: (near-)flat
kappa4 = .2^2;          % free elements of B0

[alp0,  Valp,  sig2_hat] = bvt.priors.minn(p, kappa1, kappa2, kappa3, Y0, Y, 4);
[beta0, Vbeta] = bvt.priors.impact_B0(Y0, Y, kappa4);

    % priors of the SV state equation (main_varsv.m case 3)
Hyper.nuh  = 3*ones(n,1);   Hyper.Sh   = .1*(Hyper.nuh-1);   % sig2 ~ IG(3,.2)
Hyper.mu0  = zeros(n,1);    Hyper.Vmu  = 100*ones(n,1);      % mu ~ N(0,100)
Hyper.phi0 = .98*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);    % phi ~ N(.98,.05^2)

fprintf('\nAR(4) residual variances sig2_hat = %s\n', mat2str(round(sig2_hat',3)));
fprintf('  prior sd on an own lag-1 coefficient  %.3f\n', sqrt(Valp(2)));
fprintf('  prior sd on a cross lag-1 coefficient %.3f  (kappa2/kappa1 = %.2f)\n', ...
    sqrt(Valp(3)), kappa2/kappa1);
fprintf('  prior sd on the free elements of B0   %s\n', mat2str(round(sqrt(Vbeta)',3)));

%% ------------------------------------------------------------------
%  4. Initialize the Markov chain - VAR_ARSV_redu lines 17-33, step for step.
%     A starts at a ridge estimate, h at the crude one-component log-chi^2
%     approximation of the squared residuals (bvt.sv.init_approx1N, the paper's
%     getARh_approx1N), and mu at the sample mean of log(residual^2).
%  ------------------------------------------------------------------
mu   = log(sig2_hat);
sig2 = 1./gamrnd(Hyper.nuh, 1./Hyper.Sh);
phi  = Hyper.phi0;
A    = zeros(k,n);
h    = zeros(T,n);
XX   = X'*X;
for ii = 1:n
    iValpi   = sparse(1:k, 1:k, 1./Valp((ii-1)*k+1:ii*k));
    A(:,ii)  = (XX + iValpi)\(X'*Y(:,ii));
    s2i      = (Y(:,ii) - X*A(:,ii)).^2;
    mu(ii)   = mean(log(s2i));
    h(:,ii)  = bvt.sv.init_approx1N(s2i, mu(ii), phi(ii), sig2(ii));
end
beta = zeros(k_beta,1);
B0   = eye(n);

nsim = 1000; burnin = 500;
store_A     = zeros(nsim, k_alp);
store_beta  = zeros(nsim, k_beta);
store_h     = zeros(nsim, T, n);
store_hpara = zeros(nsim, 3*n);
count_phi   = zeros(n,1);

%% ------------------------------------------------------------------
%  5. The Gibbs sampler: four blocks.
%  ------------------------------------------------------------------
fprintf('\nrunning %d sweeps (%d kept)...\n', nsim+burnin, nsim);
tic
for isim = 1:nsim + burnin

    %% ---- BLOCK 1: the reduced-form coefficients A, EQUATION BY EQUATION ----
    %                (equation ii is COLUMN ii of A, since A is k x n with
    %                 beta = vec(A); see the layout note in the header)
    %
    % This is the block worth studying. Fix the log-volatilities and B0, and
    % stack the ORTHOGONALIZED residuals of all n equations into one long
    % regression. Writing e_t = eps_t for the reduced-form error, the model says
    % B0*e_t ~ N(0, diag(exp(h_t))), so
    %
    %       vec((Y - X*A)*B0') ./ vec(exp(h/2))  ~  N(0, I_{nT}).
    %
    % Now single out column ii of A. Zero it, so that (Y - X*A) contains the
    % contribution of every OTHER column and none of this one, and read off
    % what column ii adds back: it enters equation j through B0(j,ii)*X*A(:,ii).
    % That gives regressor matrix kron(B0(:,ii), X) against the standardized
    % left-hand side above - a Gaussian linear model in the k unknowns A(:,ii).
    %
    % WHY ONLY ROWS ii:n. B0 is LOWER TRIANGULAR, so B0(j,ii) = 0 for every
    % j < ii: column ii of A cannot affect the orthogonalized errors of
    % equations 1..ii-1. Those rows of the stacked system carry an all-zero
    % regressor block, so they contribute nothing to Wi'*Wi and nothing to
    % Wi'*yi. Dropping them is EXACT, not an approximation - the posterior is
    % identical, the arithmetic is smaller. Equation ii then costs
    % (n-ii+1)*T rows instead of n*T, and the last equation costs just T.
    %
    % That truncation is what makes reduced-form estimation feasible at scale.
    % Contrast bvt.samplers.eq_svar_oi, the same block for the order-invariant
    % SVAR-SV of Chan, Koop and Yu (2024): there B0 is a general rotation with
    % no structural zeros, so it must keep all n*T rows for every equation -
    % `yi = vec((Y-X*A)*B0')./Lambda` and `Wi = kron(B0(:,ii),X)./Lambda`, with
    % the full column B0(:,ii). The two blocks differ in exactly one index range.
    %
    % Note also that A(:,ii) = 0 must be re-imposed EVERY equation, and Y-X*A
    % recomputed, because the previous equations have already been updated: this
    % is a Gibbs sweep through the equations (the columns of A), not a parallel
    % update.
    sqrt_exph = exp(h/2);
    for ii = 1:n
        A(:,ii)  = 0;
        Lambda   = bvt.util.vec(sqrt_exph(:,ii:n));           % (n-ii+1)*T x 1
        yi       = bvt.util.vec((Y - X*A)*B0(ii:n,:)')./Lambda;
        Wi       = kron(B0(ii:n,ii), X)./Lambda;              % (n-ii+1)*T x k
        iValpi   = sparse(1:k, 1:k, 1./Valp((ii-1)*k+1:ii*k));
        alpi0    = alp0((ii-1)*k+1:ii*k);                     % zero here, but the
        Kalpi    = iValpi + Wi'*Wi;                           %   term is kept so
        CKalpi   = chol(Kalpi, 'lower');                      %   this diffs against
        alpi_hat = (CKalpi')\(CKalpi\(iValpi*alpi0 + Wi'*yi));%   the legacy file
        A(:,ii)  = alpi_hat + CKalpi'\randn(k,1);
    end
    if isim == 1
        fprintf('  block 1 stacks %s rows for equations 1..%d (n*T = %d if B0 were dense)\n', ...
            mat2str(T*(n:-1:1)), n, n*T);
    end

    %% ---- BLOCK 2: the free elements of B0, equation by equation ----
    %
    % Given A, the reduced-form residuals E = Y - X*A are data. Row ii of B0
    % says eps_ii + sum_{j<ii} B0(ii,j)*eps_j has variance exp(h(:,ii)), i.e.
    % regress E(:,ii) on -E(:,1:ii-1) with weights exp(-h(:,ii)) and a normal
    % prior. bvt.samplers.alp_tri_cs is that loop, byte for byte (it was
    % extracted from the Cholesky-SV sampler of the OISV package, where the same
    % block appears under the name `alp` - the naming trap again).
    E    = Y - X*A;
    beta = bvt.samplers.alp_tri_cs(E, h, Vbeta)';
    B0(B0_id) = beta;

    %% ---- BLOCK 3: the n log-volatility paths ----
    %
    % Orthogonalize the residuals, then run one KSC auxiliary-mixture pass per
    % equation. bvt.sv.ksc_ar1_mean is ex02's sampler with a stationary AR(1)
    % state equation instead of a random walk; the call signature below matches
    % the legacy `sample_SV(ystar,h(:,ii),mu(ii),phi(ii),sig2(ii))` exactly.
    B0E = E*B0';
    for ii = 1:n
        ystar   = log(B0E(:,ii).^2 + .0001);
        h(:,ii) = bvt.sv.ksc_ar1_mean(ystar, h(:,ii), mu(ii), phi(ii), sig2(ii));
    end

    %% ---- BLOCK 4: the SV state-equation parameters (mu, phi, sig2) ----
    %
    % sig2 conjugate inverse gamma, phi an independence-chain MH step against
    % the stationary initial condition, mu Gaussian. See the caveat in the
    % header: this is the OISV copy of the block, not the ml_varsv one.
    [mu, phi, sig2, flag_phi] = bvt.sv.sv_params(h, mu, phi, Hyper);
    count_phi = count_phi + flag_phi;

    if isim > burnin
        isave = isim - burnin;
        store_A(isave,:)     = reshape(A, 1, k_alp);
        store_beta(isave,:)  = beta';
        store_h(isave,:,:)   = h;
        store_hpara(isave,:) = [mu', phi', sig2'];
    end
end
elapsed = toc;

A_hat     = reshape(mean(store_A)', k, n);
A_sd      = reshape(std(store_A)',  k, n);
beta_hat  = mean(store_beta)';
h_hat     = squeeze(mean(store_h));
hpara_hat = mean(store_hpara)';

%% ------------------------------------------------------------------
%  6. Read the output
%  ------------------------------------------------------------------
fprintf('done in %.1f seconds (%.1f ms per sweep)\n', elapsed, 1e3*elapsed/(nsim+burnin));

fprintf('\nequation 1, the whole column A(:,1):\n');
rownames = [{'intercept'}, arrayfun(@(j) sprintf('y%d(t-1)',j), 1:n, 'uni', 0), ...
    arrayfun(@(j) sprintf('y%d(t-2)',j), 1:n, 'uni', 0)];
fprintf('  %-12s %10s %12s %10s\n', 'coefficient', 'true', 'post mean', 'post sd');
for ii = 1:k
    fprintf('  %-12s %10.3f %12.3f %10.3f\n', rownames{ii}, ...
        A_true(ii,1), A_hat(ii,1), A_sd(ii,1));
end

    % Ordinary least squares on the same design is the natural comparator: it
    % keeps the reduced form but drops BOTH the prior and the SV weighting.
A_ols = X\Y;
iz  = A_true(:) == 0;
rms = @(M) [sqrt(mean((M(:)-A_true(:)).^2)), sqrt(mean(M(iz).^2)), ...
            sqrt(mean((M(~iz)-A_true(~iz)).^2))];
fprintf('\nRMSE against the truth over all %d coefficients of A:\n', k_alp);
fprintf('  %-40s %10s %10s %12s\n', '', 'all', 'true zero', 'true nonzero');
fprintf('  %-40s %10.4f %10.4f %12.4f\n', 'posterior mean (Minnesota + SV)', rms(A_hat));
fprintf('  %-40s %10.4f %10.4f %12.4f\n', 'ordinary least squares', rms(A_ols));
fprintf('  The gap sits mostly in the "true zero" column: that is the Minnesota\n');
fprintf('  prior pulling cross-lag coefficients to zero, where least squares\n');
fprintf('  scatters them. With %d observations and %d coefficients the two are\n', T, k_alp);
fprintf('  still close; at the paper''s n = 15, p = 4 (%d coefficients for ~240\n', 15^2*4+15);
fprintf('  quarters) least squares is simply unusable.\n');

fprintf('\nfree elements of the impact matrix B0:\n');
fprintf('  %-12s %10s %12s %12s\n', 'element', 'true', 'post mean', 'post sd');
blab = {'B0(2,1)','B0(3,1)','B0(3,2)'};
for ii = 1:k_beta
    fprintf('  %-12s %10.3f %12.3f %12.3f\n', blab{ii}, beta_true(ii), ...
        beta_hat(ii), std(store_beta(:,ii)));
end

fprintf('\nstochastic volatility, by equation:\n');
fprintf('  %-4s %14s %10s %20s %20s %20s\n', 'eq', 'corr(h,h_true)', 'RMSE(h)', ...
    'mu (true)', 'phi (true)', 'sig2 (true)');
for ii = 1:n
    fprintf('  %-4d %14.3f %10.3f %11.3f (%6.3f) %11.3f (%6.3f) %11.3f (%6.3f)\n', ii, ...
        corr(h_hat(:,ii), h_true(:,ii)), sqrt(mean((h_hat(:,ii)-h_true(:,ii)).^2)), ...
        hpara_hat(ii), mu_true(ii), hpara_hat(n+ii), phi_true(ii), ...
        hpara_hat(2*n+ii), sig2_true(ii));
end
h_ci = quantile(store_h, [.05 .95]);
fprintf('  share of periods with h_true inside the 90%% band: %s\n', ...
    mat2str(round(mean(squeeze(h_ci(1,:,:)) <= h_true & h_true <= squeeze(h_ci(2,:,:))), 2)));
fprintf('  phi MH acceptance rates: %s\n', ...
    mat2str(round(count_phi'/(nsim+burnin), 2)));
fprintf('  Two things to read here. mu is estimated far less sharply than phi or\n');
fprintf('  the path itself: it is the LEVEL of a near-unit-root series, so T = %d\n', T);
fprintf('  contributes few effective observations about it. And the path\n');
fprintf('  correlations sit around 0.8, not 0.99, because a single squared\n');
fprintf('  observation carries log-chi^2(1) measurement error of variance\n');
fprintf('  pi^2/2 = %.2f against a log-volatility whose own variance is about\n', pi^2/2);
fprintf('  %.2f - all the sharpness comes from smoothing across t. What the\n', ...
    mean(sig2_true./(1-phi_true.^2)));
fprintf('  coefficient draw uses is the SHAPE of h, and that is recovered well.\n');

%% ------------------------------------------------------------------
%  7. Figure (non-blocking)
%  ------------------------------------------------------------------
figure('Name','ex04 reduced-form BVAR-SV');
for ii = 1:n
    subplot(n,1,ii)
    hq = squeeze(quantile(store_h(:,:,ii), [.05 .95]))';
    plot(1:T, h_true(:,ii), 'k', 'LineWidth', 1.2); hold on
    plot(1:T, h_hat(:,ii), 'r', 'LineWidth', 1.2);
    plot(1:T, hq, 'r:'); hold off; box off
    title(sprintf('orthogonalized error %d, log-volatility: truth (black), posterior mean + 90%% band (red)', ii))
end
drawnow

fprintf('\nex04 done. The published version of this sampler - 15 variables, 4 lags,\n');
fprintf('10000 draws, drawn kappa, and the marginal likelihood that the paper is\n');
fprintf('actually about - is\n');
fprintf('  replications/chan2023_joe_mlvarsv/legacy/VAR_ARSV_redu.m  (via main_varsv.m)\n');
fprintf('Read it next: this loop is its lines 40-90 with the kappa block removed.\n');
fprintf('That package has NOT been functionized yet - there is no run_all.m under\n');
fprintf('replications/chan2023_joe_mlvarsv/, only the frozen legacy/ folder, so run\n');
fprintf('it the legacy way (cd into legacy/, then main_varsv).\n');
fprintf('Next: ex05_marginal_likelihood.\n');
