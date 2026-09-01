% bvt.priors.resid_var_ar4 - residual variances of univariate AR(4) models,
% used to set the Minnesota-prior scalings.
% Extracted 2026-09-01 (step 4, SV/prior core): verified identical (modulo header
% comments and a closing `end`; comment-stripped diff + unit test) across the four
% legacy copies of get_resid_var.m:
%   chan2021_ijf_mahp/legacy/get_resid_var.m (canonical),
%   chan2019wp_acp/legacy/get_resid_var.m,
%   chan2022_qe_acp/legacy/utility/get_resid_var.m,
%   chan_koop_yu2024_jbes_oisv/legacy/utility/get_resid_var.m.
% Function renamed get_resid_var -> resid_var_ar4; body verbatim, nothing
% parameterized. NEVER merge with bvt.priors.resid_var_allvars_ridge (legacy
% get_resid_var_v2, HYB): that one regresses each variable on 4 lags of ALL
% variables with a 1e-4 ridge - numerically different sig2.
%
% This function computes the residuals of univariate AR(4) models
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3): 1212-1226

function sig2 = resid_var_ar4(Y0,Y)
[T,n] = size(Y);
sig2 = zeros(n,1);
tmpY = [Y0(end-4+1:end,:); Y];
for i=1:n
    Z = [ones(T,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i) tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
end
end
