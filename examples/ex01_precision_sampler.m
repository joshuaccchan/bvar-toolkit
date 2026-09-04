%% ex01 - The precision sampler: drawing a whole state path in one block
%
% This is the computational foundation of everything else in the toolkit.
%
% THE IDEA (Chan and Jeliazkov, 2009). A linear Gaussian state space model can
% be written without any filtering recursion at all. Stack the states into one
% long vector and note that the state equation
%
%       tau_t = tau_{t-1} + u_t,     u_t ~ N(0, sig2_tau)
%
% is, for the whole path at once, just a linear map:  H*tau = alpha + u, with
%
%       H = |  1                |          (first difference matrix)
%           | -1   1            |
%           |     -1   1        |
%           |          .   .    |
%
% so tau ~ N(H\alpha, (H' S^-1 H)^-1). The prior precision K_tau = H'S^-1H is
% TRIDIAGONAL, and adding the (diagonal) measurement precision keeps it
% tridiagonal. A banded Cholesky factorization costs O(T) flops, not O(T^3),
% and the whole path is drawn in ONE multivariate normal draw. No Kalman
% filter, no forward-backward pass, no loop over t.
%
% Everything in core/+bvar/+sv/ (bvar.sv.ksc_rw_h0, bvar.sv.ksc_rw_diffuse,
% bvar.sv.ksc_ar1_mean, bvar.sv.csv_armh) is this same three-line construction
% with a different H and a different diagonal.
%
% WHAT TO LOOK AT when you run this:
%   1. the sparsity numbers: K has ~3T nonzeros out of T^2;
%   2. the timing table: doubling T roughly doubles the time (sparse) but
%      multiplies it by ~8 (dense) - that is O(T) versus O(T^3);
%   3. the last section: bvar.util.surform builds the same machinery for a
%      k-dimensional state (a time-varying-parameter regression).
%
% See: Chan, J.C.C. and I. Jeliazkov (2009). Efficient Simulation and
% Integrated Likelihood Estimation in State Space Models, International
% Journal of Mathematical Modelling and Numerical Optimisation, 1: 101-120.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

rng(20260903, 'twister')
fprintf('\n=== ex01: the precision sampler on a local-level model ===\n');

%% ------------------------------------------------------------------
%  1. Simulate a local-level (unobserved-components) model
%     y_t   = tau_t + eps_t,      eps_t ~ N(0, sig2_eps)     [measurement]
%     tau_t = tau_{t-1} + u_t,    u_t   ~ N(0, sig2_tau)     [state]
%     tau_1 ~ N(tau0, sig2_tau)                              [known start]
%  ------------------------------------------------------------------
T        = 500;
sig2_eps = 1.00;          % noise variance   (large: the trend is hard to see)
sig2_tau = 0.02;          % state variance   (small: the trend is smooth)
tau0     = 0;             % known initial level

tau_true = tau0 + cumsum(sqrt(sig2_tau)*randn(T,1));
y        = tau_true + sqrt(sig2_eps)*randn(T,1);

fprintf('\nsimulated T = %d observations; signal-to-noise sig2_tau/sig2_eps = %.3f\n', ...
    T, sig2_tau/sig2_eps);

%% ------------------------------------------------------------------
%  2. Build the banded precision matrix and draw tau in one block
%  ------------------------------------------------------------------
% H is the first difference matrix. Note it is built with sparse(), never
% with eye(T) - a dense T x T identity is exactly what we are avoiding.
H = speye(T) - sparse(2:T, 1:T-1, ones(T-1,1), T, T);

% alpha carries the initial condition into the stacked state equation:
% H*tau = [tau0; 0; ...; 0] + u, so the state prior mean is H\[tau0;0;...].
alph = H\[tau0; sparse(T-1,1)];

invS  = spdiags(1/sig2_tau*ones(T,1), 0, T, T);   % state precision
K_tau = H'*invS*H;                                % PRIOR precision: tridiagonal
invOm = spdiags(1/sig2_eps*ones(T,1), 0, T, T);   % measurement precision
K     = K_tau + invOm;                            % POSTERIOR precision: still tridiagonal

% This is the whole sampler. Three lines.
C      = chol(K, 'lower');                        % banded: L has 2 nonzero bands
tau_hat = K\(K_tau*alph + invOm*y);               % posterior mean of the path
tau_draw = tau_hat + C'\randn(T,1);               % one draw of the WHOLE path

fprintf('\nsparsity of the T x T posterior precision K (T = %d):\n', T);
fprintf('  nonzeros in K            : %6d   (%.4f%% of T^2 = %d)\n', ...
    nnz(K), 100*nnz(K)/T^2, T^2);
fprintf('  nonzeros in chol(K)      : %6d   (a dense factor would have %d)\n', ...
    nnz(C), T*(T+1)/2);
fprintf('  top-left 4x4 block of K  :\n');
disp(full(K(1:4,1:4)));
fprintf('  -> every row has at most 3 nonzeros. That is what makes it O(T).\n');

%% ------------------------------------------------------------------
%  3. Put it inside a Gibbs sampler: also learn the two variances
%     Priors: sig2_eps ~ IG(nu_eps0, S_eps0),  sig2_tau ~ IG(nu_tau0, S_tau0)
%  ------------------------------------------------------------------
nu_eps0 = 3;  S_eps0 = 1*(nu_eps0-1);      % prior mean 1
nu_tau0 = 3;  S_tau0 = 0.1*(nu_tau0-1);    % prior mean 0.1

nsim = 400; burnin = 100;
store_tau   = zeros(nsim, T);
store_sig2  = zeros(nsim, 2);
sig2_eps_d  = 1;  sig2_tau_d = 0.1;        % chain initial values

for isim = 1:nsim+burnin
        % (a) draw the whole path tau | y, sig2_eps, sig2_tau   <- the block above
    invS  = spdiags(1/sig2_tau_d*ones(T,1), 0, T, T);
    K_tau = H'*invS*H;
    invOm = spdiags(1/sig2_eps_d*ones(T,1), 0, T, T);
    K     = K_tau + invOm;
    C     = chol(K, 'lower');
    tau_h = K\(K_tau*alph + invOm*y);
    tau   = tau_h + C'\randn(T,1);

        % (b) draw sig2_eps | tau, y      (conjugate inverse gamma)
    e_eps      = y - tau;
    sig2_eps_d = 1/gamrnd(nu_eps0 + T/2, 1/(S_eps0 + sum(e_eps.^2)/2));

        % (c) draw sig2_tau | tau         (conjugate inverse gamma)
    e_tau      = tau - [tau0; tau(1:T-1)];
    sig2_tau_d = 1/gamrnd(nu_tau0 + T/2, 1/(S_tau0 + sum(e_tau.^2)/2));

    if isim > burnin
        store_tau(isim-burnin,:)  = tau';
        store_sig2(isim-burnin,:) = [sig2_eps_d sig2_tau_d];
    end
end

tau_mean = mean(store_tau)';
tau_ci   = quantile(store_tau, [.05 .95])';
cover    = mean(tau_true >= tau_ci(:,1) & tau_true <= tau_ci(:,2));

fprintf('\nGibbs output (%d draws after %d burn-in):\n', nsim, burnin);
fprintf('  sig2_eps  true %.3f   posterior mean %.3f   90%% CI [%.3f, %.3f]\n', ...
    sig2_eps, mean(store_sig2(:,1)), quantile(store_sig2(:,1),.05), quantile(store_sig2(:,1),.95));
fprintf('  sig2_tau  true %.3f   posterior mean %.3f   90%% CI [%.3f, %.3f]\n', ...
    sig2_tau, mean(store_sig2(:,2)), quantile(store_sig2(:,2),.05), quantile(store_sig2(:,2),.95));
fprintf('  RMSE of the smoothed trend vs the truth : %.4f\n', ...
    sqrt(mean((tau_mean-tau_true).^2)));
fprintf('  RMSE of the raw data      vs the truth : %.4f  (the filter is doing work)\n', ...
    sqrt(mean((y-tau_true).^2)));
fprintf('  share of periods where the truth lies in the 90%% band: %.2f\n', cover);

%% ------------------------------------------------------------------
%  4. WHY IT IS O(T): time one draw as T grows, sparse versus dense
%  ------------------------------------------------------------------
fprintf('\ncost of ONE state-path draw as T grows (median of 20 draws):\n');
fprintf('  %8s %12s %10s %12s %10s\n', 'T', 'sparse (ms)', 'ratio', 'dense (ms)', 'ratio');
Tgrid   = [500 1000 2000 4000 8000 16000];
t_sp    = nan(size(Tgrid));
t_de    = nan(size(Tgrid));
for ii = 1:numel(Tgrid)
    Ti  = Tgrid(ii);
    Hi  = speye(Ti) - sparse(2:Ti, 1:Ti-1, ones(Ti-1,1), Ti, Ti);
    Ki  = Hi'*spdiags(1/sig2_tau*ones(Ti,1),0,Ti,Ti)*Hi ...
        + spdiags(1/sig2_eps*ones(Ti,1),0,Ti,Ti);
    ri  = randn(Ti,1);
    reps = zeros(20,1);
    for r = 1:20
        tic; Ci = chol(Ki,'lower'); x = Ci'\ri; reps(r) = toc;
    end
    t_sp(ii) = 1e3*median(reps);

        % the same draw with a DENSE factorization, for contrast only
    if Ti <= 2000
        Kd = full(Ki);
        reps = zeros(5,1);
        for r = 1:5
            tic; Cd = chol(Kd,'lower'); x = Cd'\ri; reps(r) = toc;
        end
        t_de(ii) = 1e3*median(reps);
    end
    rsp = ''; rde = '';
    if ii > 1
        rsp = sprintf('%.1fx', t_sp(ii)/t_sp(ii-1));
        if ~isnan(t_de(ii)) && ~isnan(t_de(ii-1)), rde = sprintf('%.1fx', t_de(ii)/t_de(ii-1)); end
    end
    if isnan(t_de(ii))
        fprintf('  %8d %12.3f %10s %12s %10s\n', Ti, t_sp(ii), rsp, '(skipped)', '');
    else
        fprintf('  %8d %12.3f %10s %12.3f %10s\n', Ti, t_sp(ii), rsp, t_de(ii), rde);
    end
end
fprintf('  read the ratio columns: doubling T roughly DOUBLES the sparse cost\n');
fprintf('  (linear), while the dense cost grows several times over - it is\n');
fprintf('  headed for the 8x per doubling of an O(T^3) factorization, and only\n');
fprintf('  LAPACK''s fixed overheads keep it below that at these small sizes.\n');
fprintf('  At T = 16000 the dense route needs a %.1f GB matrix, which is why it\n', ...
    8*16000^2/2^30);
fprintf('  is skipped above; the sparse draw there still takes under a\n');
fprintf('  millisecond. That is the difference between a feasible and an\n');
fprintf('  infeasible sampler.\n');

%% ------------------------------------------------------------------
%  5. The same trick with a k-DIMENSIONAL state: bvar.util.surform
%     A time-varying-parameter regression
%       y_t = b_{1t} + b_{2t}*z_t + e_t,   b_t = b_{t-1} + u_t
%     bvar.util.surform(X) turns the T x k matrix of regressors into the
%     T x kT block-diagonal matrix that multiplies the STACKED state
%     vector [b_1; b_2; ...; b_T]. Precision is then banded with bandwidth
%     2k instead of 2 - still O(T).
%  ------------------------------------------------------------------
T2 = 300; kk = 2;
z  = randn(T2,1);
Xt = [ones(T2,1) z];
b0 = [0; 1];                       % known initial coefficient vector
sig2_b = [0.005; 0.02];            % state variances, one per coefficient
sig2_e = 0.25;

b_true = zeros(T2, kk);
b_true(1,:) = (b0 + sqrt(sig2_b).*randn(kk,1))';
for t = 2:T2
    b_true(t,:) = b_true(t-1,:) + (sqrt(sig2_b).*randn(kk,1))';
end
y2 = sum(Xt.*b_true, 2) + sqrt(sig2_e)*randn(T2,1);

Xbig = bvar.util.surform(Xt);                     % T2 x (kk*T2), sparse
Hb   = speye(kk*T2) - sparse(kk+1:kk*T2, 1:kk*T2-kk, ones(kk*T2-kk,1), kk*T2, kk*T2);
invSb = spdiags(repmat(1./sig2_b, T2, 1), 0, kk*T2, kk*T2);
Kb    = Hb'*invSb*Hb;
alphb = Hb\[b0; sparse(kk*T2-kk,1)];
invOm2 = spdiags(1/sig2_e*ones(T2,1), 0, T2, T2);
Kpost  = Kb + Xbig'*invOm2*Xbig;
bhat   = Kpost\(Kb*alphb + Xbig'*invOm2*y2);
Cb     = chol(Kpost,'lower');

nrep = 200; store_b = zeros(nrep, kk*T2);
for r = 1:nrep
    store_b(r,:) = (bhat + Cb'\randn(kk*T2,1))';
end
b_mean = reshape(mean(store_b)', kk, T2)';       % T2 x kk posterior mean paths

fprintf('\nTVP regression with a %d-dimensional state, T = %d:\n', kk, T2);
fprintf('  bvar.util.surform(Xt) is %d x %d with %d nonzeros (%.3f%% dense)\n', ...
    size(Xbig,1), size(Xbig,2), nnz(Xbig), 100*nnz(Xbig)/numel(Xbig));
fprintf('  posterior precision Kpost: %d x %d, %d nonzeros, bandwidth %d\n', ...
    size(Kpost,1), size(Kpost,2), nnz(Kpost), bandwidth(Kpost,'lower'));
fprintf('  RMSE of the smoothed coefficient paths vs the truth:\n');
fprintf('    intercept b_1t : %.4f\n', sqrt(mean((b_mean(:,1)-b_true(:,1)).^2)));
fprintf('    slope     b_2t : %.4f\n', sqrt(mean((b_mean(:,2)-b_true(:,2)).^2)));

%% ------------------------------------------------------------------
%  6. Figures (non-blocking)
%  ------------------------------------------------------------------
figure('Name','ex01 precision sampler');
subplot(2,2,1)
plot(1:T, y, 'Color', [.75 .75 .75]); hold on
plot(1:T, tau_true, 'k', 'LineWidth', 1.2);
plot(1:T, tau_mean, 'r', 'LineWidth', 1.2);
plot(1:T, tau_ci, 'r:');
hold off; box off
title('local level: data, truth, posterior mean and 90% band')
legend({'y_t','\tau_t true','posterior mean'}, 'Location','best'); legend boxoff

subplot(2,2,2)
spy(K(1:60,1:60)); title('sparsity of K (first 60 x 60 block)')

subplot(2,2,3)
loglog(Tgrid, t_sp, 'o-'); hold on
loglog(Tgrid(~isnan(t_de)), t_de(~isnan(t_de)), 's-'); hold off; box off
xlabel('T'); ylabel('ms per draw'); title('cost of one path draw')
legend({'sparse','dense'}, 'Location','best'); legend boxoff

subplot(2,2,4)
plot(1:T2, b_true, 'k'); hold on
plot(1:T2, b_mean, 'r'); hold off; box off
title('TVP coefficients: truth (black) vs smoothed (red)')

drawnow

fprintf('\nex01 done. Next: ex02_sv_ksc (stochastic volatility on top of this).\n');
