% bvt.sv.ksc_ar1_mean - KSC auxiliary-mixture sampler for the log-volatility
% path, stationary AR(1)-with-mean state equation:
%   ystar_t = h_t + eps_t,  eps_t approximated by the Kim-Shephard-Chib (1998)
%             7-component normal mixture,
%   h_t = mu + rho*(h_{t-1} - mu) + v_t,  v_t ~ N(0,sig2),
%   h_1 ~ N(mu, sig2/(1-rho^2))  (stationary initialization).
% Returns the new path h AND the mixture indicators S.
% Consumes rand(T,1) then randn(T,1) - one of each per call.
%
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical body:
% chan2023_joe_mlvarsv/legacy/utility/sample_SV.m, verbatim. Also canonicalizes
% chan_koop_yu2024_jbes_oisv/legacy/utility/sample_SV.m (header-comment-only
% diff; bodies byte-identical, draw-for-draw bitwise equality verified in
% R2025b and by tests/unit/test_ksc_ar1_mean.m).
% Function renamed sample_SV -> ksc_ar1_mean; nothing parameterized.
% NEVER merge with the random-walk variants (bvt.sv.ksc_rw_h0 / ksc_rw_diffuse):
% different state equations. The companion parameter samplers sample_SVpara /
% sample_SV0para differ across packages and are NOT consolidated here (see
% tests/variant_map.md never-merge list).
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [h,S] = ksc_ar1_mean(ystar,h,mu,rho,sig2)
T = length(h);
    % 7-component normal mixture
p_N = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
m_N = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  % means already adjusted!!
sig2_N = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];

    % sample S from a 7-point distrete distribution
tmprand = rand(T,1);
q = repmat(p_N,T,1).*normpdf(repmat(ystar,1,7),repmat(h,1,7)+repmat(m_N,T,1),...
    repmat(sqrt(sig2_N),T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(tmprand,1,7)<cumsum(q,2),2)+1;

    % sample h using the precision sampler
% y^* = h + d + \epsilon, \epsilon \sim N(0,\Omega)
% h  ~ N(mu 1_T, sig2(Hrho' S^{-1} Hrho)^{-1}),
% where S^{-1} = diag(1-rho^2,1,1,...,1)

Hrho = speye(T) - sparse(2:T,1:(T-1),rho*ones(1,T-1),T,T);
HiSH = Hrho'*sparse(1:T,1:T,[1-rho^2, ones(1,T-1)])*Hrho;
d = m_N(S)'; iOmega = sparse(1:T,1:T,1./sig2_N(S));
Kh = HiSH/sig2 + iOmega;
h_hat = Kh\(mu/sig2*HiSH*ones(T,1) + iOmega*(ystar-d));
h = h_hat + chol(Kh,'lower')'\randn(T,1);
end
