% bvar.samplers.eq_gauss - one full equation-by-equation sweep of the Gaussian
% coefficient draw in the structural-form BVAR of Chan (2021, IJF): for each
% equation ii = 1..n, draw thetai = (beta_ii; alp_ii) jointly from its Gaussian
% conditional through the sparse precision
%   Kthetai = iVthetai + Xi' * diag(exp(-h(:,ii))) * Xi,   Xi = [Z -Y(:,1:ii-1)],
% via chol(Kthetai,'lower'), and accumulate the structural residuals U.
% rng consumption: exactly one randn(ki,1) per equation, ki = n*p+ii - nothing else.
%
% Body from chan2021_ijf_mahp/legacy/BVAR_MNG.m lines 40-59 (the inline "sample
% alp and beta" block), wrapped as a function (Y,Z,h,Valp,Vbeta in; beta,alp,U
% out) with np = size(Z,2)-1 replacing the literal n*p. The same block appears
% verbatim in BVAR_NG.m, BVAR_Minn.m and the three forecast_BVAR_*.m scripts
% (Y/Z/T -> Yt/Zt/Tt only). Prior-variance scaling stays with the CALLER:
% BVAR_MNG and all three forecast scripts double Valp/Vbeta before this block,
% estimation BVAR_NG/BVAR_Minn do not - pass Valp/Vbeta exactly as the legacy
% script has them at entry to the equation loop.
% Equivalence: tests/unit/test_mahp_equivalence.m. Record: tests/variant_map.md.
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
