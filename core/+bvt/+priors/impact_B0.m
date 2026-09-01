% bvt.priors.impact_B0 - data-based prior variances for the free elements of
% the impact matrix B0 (AR(4) residual-variance scaling, Minnesota-style).
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical source:
% replications/chan2023_joe_mlvarsv/legacy/utility/prior_B0.m (single copy,
% glob-verified). Body verbatim; only the function was renamed
% prior_B0 -> impact_B0. No parameterization was added.
% Role-equivalent constructions elsewhere that are NOT canonicalized here
% (numerically different priors on the same object):
% - chan_koop_yu2024_jbes_oisv: inline Hyper.B0 = eye(n), Hyper.VB0 = ones(n)
%   (full n x n impact matrix, unit prior variances, no data-based scaling);
% - chan2021_ijf_mahp: Valp via getVtheta.m, kappa(3)*sig2(ii)/sig2(1:ii-1)
%   with sig2 supplied by the caller (and the HYB getVtheta copy hard-codes
%   kappa3 = .2 - see the never-merge entry for getVtheta.m);
% - chan2022_qe_acp / chan2019wp_acp: the kappa(3)/sig2(idx) term inside
%   prior_ACPi (no sig2(i) numerator, sig2 from the caller's AR fits).
%
% This function constructs the prior for the impact matrix B0
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

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
