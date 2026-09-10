% bvar.priors.vtheta - conditional Minnesota-type prior variances of the VAR
% coefficients (Vbeta) and the free elements of the impact matrix (Valp), given
% the shrinkage hyperparameters kappa and the AR(4) residual variances sig2.
%
%   [Valp,Vbeta] = bvar.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C, sig2)
%
%   idx_kappa1, idx_kappa2, C : the own-lag indices, other-lag indices and
%            second moments returned by bvar.priors.minnesota_C
%   kappa  : ALWAYS a 4-vector - kappa(1) own lags, kappa(2) other lags,
%            kappa(3) impact matrix, kappa(4) intercepts. Callers that fix the
%            last two (the hybrid TVP-VAR setting kappa3 = .2, kappa4 = 1) pass
%            [kappa1, kappa2, .2, 1].
%   sig2   : n-vector of AR(4) residual variances
%   Vbeta  : prior variances of the VAR coefficients, stacked equation by
%            equation, intercept first in each block
%   Valp   : prior variances of the n(n-1)/2 free elements of the impact matrix,
%            kappa(3)*sig2_i/sig2_j. A caller wanting only Vbeta may pass
%            kappa(3) = NaN and take [~,Vbeta]; Valp is then NaN throughout.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_vtheta.m, tests/unit/test_oisv_equivalence.m.
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3): 1212-1226

function [Valp,Vbeta] = vtheta(idx_kappa1,idx_kappa2,kappa,C,sig2)
np = length(idx_kappa1);
n = length(idx_kappa2)/np + 1;
k_beta = length(C);
k_alp = n*(n-1)/2;
Vbeta = zeros(k_beta,1);
Valp = zeros(k_alp,1);

Vbeta(1:np+1:end) = kappa(4)*sig2;          % intercepts
Vbeta(idx_kappa1) = kappa(1)*C(idx_kappa1); % own lags
Vbeta(idx_kappa2) = kappa(2)*C(idx_kappa2); % other lags

count_alp = 0;
for ii = 1:n
    Valpi = kappa(3)*repmat(sig2(ii),ii-1,1)./sig2(1:ii-1);
    Valp(count_alp+1:count_alp+ii-1) = Valpi;
    count_alp = count_alp + ii - 1;
end
end
