% bvar.ml.lmvnpdf_pcn - log density of N(mu, inv(K)) at x, parameterized by the
% PRECISION K (not the covariance).
%
% Body from chan2023_joe_mlvarsv/legacy/utility/lmvnpdf_pcn.m.
% Equivalence: tests/unit/test_mlvarsv_ml_densities.m. Record: tests/variant_map.md.
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
