% bvar.ml.lgampdf - log density of the gamma distribution at x (shape a, RATE
% b; elementwise over array inputs).
%
% Extracted 2026-09-03 (step 10, ml_varsv marginal-likelihood pass). Canonical
% source: chan2023_joe_mlvarsv/legacy/utility/lgampdf.m (body verbatim; the
% only copy). Used by the kappa prior and importance-sampling ordinates of
% bvar.ml.mlvarsv_csv. Rate, not scale: gamfit returns a scale, and the legacy
% callers invert it before calling (ckappa_hat = [a; 1/scale]).
% Verified bitwise by tests/unit/test_mlvarsv_ml_densities.m.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function lden = lgampdf(x,a,b)
    lden = a.*log(b) -gammaln(a) +(a-1).*log(x) -b.*x;
end
