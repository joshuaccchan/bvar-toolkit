% bvar.sv.ksc_rw_h0 - KSC auxiliary-mixture sampler for the log-volatility path,
% random-walk state equation with KNOWN initial value h0 (h_1 ~ N(h0,sig)):
%   ystar_t = h_t + eps_t,  eps_t = log of a chi^2_1 draw, approximated by the
%             Kim-Shephard-Chib (1998) 7-component normal mixture,
%   h_t = h_{t-1} + v_t,    v_t ~ N(0,sig),   with sig a VARIANCE.
% Consumes rand(T,1) then randn(T,1) - one of each per call.
%
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical body:
% chan2020_springer_largebvar/legacy/SVRW.m, verbatim. Also canonicalizes
%   - chan2020_jbes_kronecker/legacy/realtime_forecasts/SVRW.m (comment-only diff),
%   - chan2021_ijf_mahp/legacy/SVRW.m (comment-only diff),
%   - chan2023_jbes_hybtvp/legacy/utility/sample_SVRW.m (algebraically identical
%     code path: its sparse()/spdiags() constructors build bitwise-identical
%     matrices and its chol(Kh,'lower')' equals this chol(Ph) bitwise for sparse
%     inputs (CHOLMOD); verified draw-for-draw bitwise in R2025b over 20 seeded
%     trials and by tests/unit/test_ksc_rw_h0.m).
% Function renamed SVRW -> ksc_rw_h0; nothing parameterized.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

function h = ksc_rw_h0(Ystar,h,sig,h0)
T = length(h);
%% normal mixture
pi = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mi = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sqrtsigi = sqrt(sigi);

%% sample S from a 7-point distrete distribution
temprand = rand(T,1);
q = repmat(pi,T,1).*normpdf(repmat(Ystar,1,7),repmat(h,1,7)+repmat(mi,T,1), repmat(sqrtsigi,T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2)+1;

%% sample h
% y^* = h + d + \epsilon, \epsilon \sim N(0,\Omega),
% Hh = \alpha + \nu, \nu \ sim N(0,S),
% where d_t = Ez_t, \Omega = diag(\omega_1,\ldots,\omega_n),
% \omega_t = var z_t, S = diag(sig, \ldots, sig)
Hh =  speye(T) - spdiags(ones(T-1,1),-1,T,T);
invSh = spdiags(1/sig*ones(T,1),0,T,T);
dconst = mi(S)'; invOmega = spdiags(1./sigi(S)',0,T,T);
alph = Hh\[h0;sparse(T-1,1)];
Kh = Hh'*invSh*Hh;
Ph = Kh + invOmega;
Ch = chol(Ph);
hhat = Ph\(Kh*alph + invOmega*(Ystar-dconst));
h = hhat + Ch\randn(T,1);
end
