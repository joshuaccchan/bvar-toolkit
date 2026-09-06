% bvar.samplers.eq_var_oi - equation-by-equation Gaussian draw of the VAR
% coefficient matrix A (k x n, intercept first) in the order-invariant VAR-SV
% model
%     y_t = A' x_t + u_t,   B0 u_t = eps_t,   eps_t ~ N(0, diag(exp(h_t))),
% i.e. the same conditional as bvar.samplers.eq_svar_oi, computed in
% O(T k^2 + k^3) per equation instead of O(T n k^2 + k^3).
%
%   A = bvar.samplers.eq_var_oi(Y, X, B0, h, A, tmpdV)
%
%   Y     : T x n observations
%   X     : T x k regressors [1, y_{t-1}, ..., y_{t-p}] (bvar.util.build_lags)
%   B0    : n x n INVERSE impact matrix - it maps reduced-form innovations to
%           structural ones, eps_t = B0 u_t, so the impact matrix itself is
%           B0^{-1} and the reduced-form covariance is B0^{-1} diag(exp(h_t))
%           B0^{-T}. Pass B0, not B0^{-1}.
%   h     : T x n log-variances of the structural innovations, column j = eps_jt
%   A     : current k x n coefficients (columns updated in place, in order)
%   tmpdV : k*n stacked prior variances, column-major (zero prior mean; to use
%           a prior mean A0, call on Y - X*A0 with A - A0 and add A0 back -
%           an exact reparameterization, tested in tests/unit/test_eq_var_oi.m)
%
% rng consumption: randn(k,1) per equation, equations in order ii = 1:n -
% identical to eq_svar_oi, so under a common seed the two functions return the
% same draw to floating-point precision (tests/unit/test_eq_var_oi.m).
%
% New in the consolidated toolkit (2026-09-06); it canonicalizes no legacy file
% and replaces nothing. eq_svar_oi remains the verbatim CKY24 block
% (SVARSV_MH.m lines 76-87) and the bitwise anchor of the
% chan_koop_yu2024_jbes_oisv replication; this function is for new code. The
% name says what is drawn - the reduced-form coefficients A - whereas "svar" in
% the legacy name refers to the structural parameterization of the error
% covariance, not to the coefficients.
%
% Why it is faster. Column ii of A enters structural equation j with the
% coefficient B0(j,ii), so its conditional precision is
%     sum_j B0(j,ii)^2 X' diag(exp(-h_j)) X  =  X' diag(w) X,
%     w_t = sum_j B0(j,ii)^2 exp(-h_jt),
% a single weighted cross-product with T rows, and its right-hand side is
% X' c with c_t = sum_j B0(j,ii) exp(-h_jt) etil_jt, where etil_t = B0 (y_t -
% A_{-ii}' x_t). eq_svar_oi writes the same system as a stacked SUR,
% Wi = kron(B0(:,ii), X)./Lambda, a dense (T n) x k matrix, and forms Wi'*Wi
% with T n rows: the factor n is the redundancy of stacking the same X once
% per structural equation. Both then pay the same O(k^3) Cholesky, which is
% what keeps the observed ratio well below n and can dominate the sweep once
% k is large relative to T. Measured ratios, and the shapes they were measured
% at, are tabulated in tests/variant_map.md; the unit test prints the figures
% for the machine it runs on. At macroeconomic dimensions the saving is about
% an order of magnitude.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function A = eq_var_oi(Y, X, B0, h, A, tmpdV)
[T, n] = size(Y);
k = size(X, 2);
if ~isequal(size(h), [T n])
    error('bvar:samplers:eq_var_oi:badH', ...
        'h must be %d x %d (T x n), not %d x %d', T, n, size(h,1), size(h,2));
end
if ~isequal(size(B0), [n n])
    error('bvar:samplers:eq_var_oi:badB0', ...
        'B0 must be %d x %d, not %d x %d', n, n, size(B0,1), size(B0,2));
end
if numel(tmpdV) ~= k*n
    error('bvar:samplers:eq_var_oi:badV', ...
        'tmpdV must have k*n = %d elements, not %d', k*n, numel(tmpdV));
end
eh_inv = exp(-h);                          % T x n, column j = exp(-h_jt)
for ii = 1:n
    A(:, ii) = 0;
    Etil = (Y - X*A)*B0';                  % row t = (B0 (y_t - A_{-ii}' x_t))'
    b = B0(:, ii);                         % loadings of A(:,ii) in the n structural equations
    w = eh_inv*(b.^2);                     % T x 1 precision weights
    c = (eh_inv.*Etil)*b;                  % T x 1
    iVi = 1./tmpdV((ii-1)*k+1:ii*k);
    Kai = X'*(w.*X);
    Kai(1:k+1:end) = Kai(1:k+1:end) + iVi(:)';
    CKai = chol(Kai, 'lower');
    ai_hat = CKai'\(CKai\(X'*c));
    A(:, ii) = ai_hat + CKai'\randn(k, 1);
end
end
