% bvar.priors.impact_B0 - data-based prior variances for the free elements of
% the impact matrix B0 (AR(4) residual-variance scaling, Minnesota-style).
%
% Body from chan2023_joe_mlvarsv/legacy/utility/prior_B0.m,
% renamed. Do not merge with the role-equivalent impact-matrix priors elsewhere -
% oisv (inline unit prior variances, no data-based scaling), mahp/hybtvp (Valp via
% getVtheta) and the acp packages (the kappa(3)/sig2(idx) term in prior_ACPi):
% they are numerically different priors on the same object.
% Equivalence: tests/unit/test_impact_B0.m. Record: tests/variant_map.md.
%
% This function constructs the prior for the impact matrix B0
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
