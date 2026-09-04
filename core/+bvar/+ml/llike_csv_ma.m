% bvar.ml.llike_csv_ma - log likelihood of the BVAR-CSV-MA(1) model at MA
% coefficient psi, given untransformed residuals U, error covariance Sig and
% the log-volatility path h (common stochastic volatility exp(h_t) scaling):
% includes the -n/2*sum(h) volatility normalizing term and the first
% observation's (1+psi^2) initialization variance.
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/llike_CSV_MA.m  (the package-ROOT copy;
%   body verbatim).
%
% NEVER-MERGE / PATH HAZARD (recorded in tests/variant_map.md): the same
% package ships a SECOND llike_CSV_MA.m in legacy/realtime_forecasts/ (and
% chan2020_springer_largebvar carries the same reduced form) that OMITS the
% -n/2*sum(h) term. Inside a psi-MH at fixed h the term is constant in psi
% and cancels from acceptance ratios and normalized grid densities
% (mathematically interchangeable, bitwise different through fminunc/fminbnd),
% but as a likelihood ORDINATE the two differ by n/2*sum(h) - the marginal-
% likelihood path requires THIS root version. The legacy main_forecasting.m
% line 20 runs addpath('./realtime_forecasts'), after which the resolution of
% an unqualified llike_CSV_MA call is path-order/cwd-dependent (the root copy
% wins only while the package root is the current folder). This core function
% pins the root semantics; all bvar.ml consumers call it fully qualified.
% Verified by tests/unit/test_kron_ml_densities.m (bitwise against the root
% copy, and asserted to differ from the realtime copy by exactly n/2*sum(h))
% and end-to-end through the ML equivalence tests.
%
% The t-model reduced runs reuse this same function with h := log(lam)
% (legacy ml_BVAR_t_MA.m line 82) and with U pre-scaled by sqrt(lam)
% (legacy ml_BVAR_CSV_t_MA.m lines 107-108), exactly as the legacy scripts do.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function llike = llike_csv_ma(psi,U,Sig,h)
[T,n] = size(U);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Utld = Hpsi\U;
CSig = chol(Sig,'lower');
c = -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) - n/2*log(1+psi^2) -n/2*sum(h);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
llike = c -.5*(s2(1)/((1+psi^2)*exp(h(1))) + sum(s2(2:end)./exp(h(2:end))));
end
