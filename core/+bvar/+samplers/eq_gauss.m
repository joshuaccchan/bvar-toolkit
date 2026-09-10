% bvar.samplers.eq_gauss - one full equation-by-equation sweep of the Gaussian
% coefficient draw in the structural-form BVAR of Chan (2021, IJF): for each
% equation ii = 1..n, draw thetai = (beta_ii; alp_ii) jointly from its Gaussian
% conditional through the sparse precision
%   Kthetai = iVthetai + Xi' * diag(exp(-h(:,ii))) * Xi,   Xi = [Z -Y(:,1:ii-1)],
% via chol(Kthetai,'lower'), and accumulate the structural residuals U.
% rng consumption: exactly one randn(ki,1) per equation, ki = n*p+ii - nothing else.
%
%   [beta, alp, U] = bvar.samplers.eq_gauss(Y, Z, h, Valp, Vbeta)
%
% Inputs:  Y     - T x n observations (equation ii regresses Y(:,ii) on
%                  [Z -Y(:,1:ii-1)], the structural triangular form)
%          Z     - T x (n*p+1) design [1, y_{t-1}, ..., y_{t-p}]
%          h     - T x n log-volatilities
%          Valp  - n*(n-1)/2 x 1 prior variances of the impact-matrix elements
%          Vbeta - n*(n*p+1) x 1 prior variances of the VAR coefficients
% Outputs: beta  - n*(n*p+1) x 1 stacked coefficient draw
%          alp   - n*(n-1)/2 x 1 impact-matrix draw
%          U     - T x n structural residuals yi - Xi*thetai
%
% Valp and Vbeta are used exactly as passed. Any scaling of the prior variances
% - some callers double both before this step - stays with the caller.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_mahp_equivalence.m.
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3), 1212-1226.

function [beta, alp, U] = eq_gauss(Y, Z, h, Valp, Vbeta)
[T, n] = size(Y);
np = size(Z, 2) - 1;   % = n*p
beta = zeros(n*(np+1), 1);
alp = zeros(n*(n-1)/2, 1);
count_alp = 0;
U = zeros(T, n);
for ii = 1:n
    yi = Y(:, ii);
    ki = np + ii;
    Xi = [Z -Y(:, 1:ii-1)];

    iVthetai = sparse(1:ki, 1:ki, 1./[Vbeta((ii-1)*(np+1)+1:ii*(np+1)); ...
        Valp(count_alp+1:count_alp+ii-1)]);
    XiiSighi = Xi' * sparse(1:T, 1:T, exp(-h(:, ii)));
    Kthetai = iVthetai + XiiSighi * Xi;
    CKthetai = chol(Kthetai, 'lower');
    thetai_hat = (CKthetai') \ (CKthetai \ (XiiSighi * yi));
    thetai = thetai_hat + CKthetai' \ randn(ki, 1);
    U(:, ii) = yi - Xi * thetai;

    beta((ii-1)*(np+1)+1:ii*(np+1)) = thetai(1:np+1);
    alp(count_alp+1:count_alp+ii-1) = thetai(np+2:end);
    count_alp = count_alp + ii - 1;
end
end
