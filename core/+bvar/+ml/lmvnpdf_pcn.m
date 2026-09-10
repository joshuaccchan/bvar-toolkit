% bvar.ml.lmvnpdf_pcn - log density of N(mu, inv(K)) at x, parameterized by the
% PRECISION K (not the covariance).
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_mlvarsv_ml_densities.m.
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
