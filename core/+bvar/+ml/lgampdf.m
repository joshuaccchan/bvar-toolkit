% bvar.ml.lgampdf - log density of the gamma distribution at x (shape a, RATE
% b; elementwise over array inputs).
%
% Body from chan2023_joe_mlvarsv/legacy/utility/lgampdf.m. Rate, not
% scale: gamfit returns a scale, and the legacy callers invert it before calling
% (ckappa_hat = [a; 1/scale]).
% Equivalence: tests/unit/test_mlvarsv_ml_densities.m. Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function lden = lgampdf(x,a,b)
    lden = a.*log(b) -gammaln(a) +(a-1).*log(x) -b.*x;
end
