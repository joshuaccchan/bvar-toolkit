% bvt.priors.vtheta - conditional Minnesota-type prior variances of the VAR
% coefficients (Vbeta) and the free elements of the impact matrix (Valp), given
% the shrinkage hyperparameters kappa and the AR(4) residual variances sig2.
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical source (body verbatim,
% renamed getVtheta -> vtheta, MAHP semantics):
%   chan2021_ijf_mahp/legacy/getVtheta.m.
% Also canonicalizes chan2023_jbes_hybtvp/legacy/utility/getVtheta.m, whose only
% differences (comment-stripped diff) are that it hard-codes kappa_3 = .2 (impact
% matrix) and kappa_4 = 1 (intercepts) inside the body and reads only kappa(1:2).
% Parameterization: kappa here is ALWAYS a 4-vector,
%   kappa(1) own lags, kappa(2) other lags, kappa(3) impact matrix, kappa(4) intercepts.
% Documented settings reproducing each legacy copy exactly:
%   MAHP: pass its kappa 4-vector unchanged;
%   HYB : pass [kappa(1), kappa(2), .2, 1].
% (kappa(3)=.2 and kappa(4)=1 give bit-identical products to the HYB hard-coded
% constants; verified by unit test under both settings.)
%
% Added 2026-09-02 (step 7, OISV family pass): the Vbeta output also
% canonicalizes chan_koop_yu2024_jbes_oisv/legacy/utility/getVbeta.m - that
% function is exactly the three Vbeta assignment lines below (same inputs,
% same order, no Valp); OISV callers use `[~,Vbeta] = bvt.priors.vtheta(...)`
% and discard Valp (NaN under the OI kappa(3) = NaN, never read). Verified
% draw-for-draw through tests/unit/test_oisv_equivalence.m.
%
% This function constructs the conditional prior of the VAR coefficients
% and the impact matrix
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
