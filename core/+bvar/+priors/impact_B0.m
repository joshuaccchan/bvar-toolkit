% bvar.priors.impact_B0 - data-based prior variances for the free elements of
% the impact matrix B0 (AR(4) residual-variance scaling, Minnesota-style).
%
%   [beta0,Vbeta] = bvar.priors.impact_B0(Y0, Y, kappa)
%
%   Y0    : presample rows; the last 4 are prepended to Y for the univariate
%           AR(4) fits that produce the scaling variances sig2
%   Y     : T x n estimation sample
%   kappa : scalar shrinkage on the impact-matrix block
%   beta0 : prior mean of the n(n-1)/2 free elements, zeros
%   Vbeta : their prior variances, kappa*sig2_i/sig2_j, stacked equation by
%           equation (row 2 of B0 first, then row 3, ...)
%
% Do not merge with the role-equivalent impact-matrix priors elsewhere in the
% library - bvar.priors.vtheta (Valp), the kappa(3)/sig2 term inside
% bvar.priors.acp_stru, and the inline unit prior variances of the OISV package,
% which carry no data-based scaling: they are numerically different priors on
% the same object.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_impact_B0.m.
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function [beta0,Vbeta] = impact_B0(Y0,Y,kappa)
[T,n] = size(Y);
k_beta = n*(n-1)/2;
beta0 = zeros(k_beta,1);
Vbeta = zeros(k_beta,1);

sig2 = zeros(n,1);
tmpY = [Y0(end-4+1:end,:); Y];
U_hat = zeros(T,n);
for i=1:n
    Z = [ones(T,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
        tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    U_hat(:,i) = tmpY(5:end,i)-Z*tmpb;
    sig2(i) = mean(U_hat(:,i).^2);
end
count_beta = 0;
for i = 2:n
    Vbeta(count_beta+1:count_beta+i-1) = kappa*repmat(sig2(i),i-1,1)./sig2(1:i-1);
    count_beta = count_beta + i-1;
end
end
