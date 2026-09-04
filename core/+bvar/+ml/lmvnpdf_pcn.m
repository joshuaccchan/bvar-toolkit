% bvar.ml.lmvnpdf_pcn - log density of N(mu, inv(K)) at x, parameterized by the
% PRECISION K (not the covariance).
%
% Extracted 2026-09-03 (step 10, ml_varsv marginal-likelihood pass). Canonical
% source: chan2023_joe_mlvarsv/legacy/utility/lmvnpdf_pcn.m (body verbatim; the
% only copy). Used for the log-volatility importance-sampling ordinate of
% bvar.ml.mlvarsv_csv, whose density comes out of bvar.ml.isden_arss in precision
% form. It re-factors K on every call, as the legacy does.
% Verified bitwise by tests/unit/test_mlvarsv_ml_densities.m.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function lden = lmvnpdf_pcn(x,mu,K)
    n = length(mu);
    CK = chol(K,'lower');
    e = CK'*(x-mu);
    lden = -n/2*log(2*pi) + sum(log(diag(CK))) - .5*(e'*e);
end
