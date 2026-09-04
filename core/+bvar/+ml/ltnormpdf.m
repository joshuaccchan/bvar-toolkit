% bvar.ml.ltnormpdf - log density of the normal N(mu,sig2) truncated to (lb,ub).
%
% Extracted 2026-09-03 (step 10, ml_varsv marginal-likelihood pass). Canonical
% source: chan2023_joe_mlvarsv/legacy/utility/ltnormpdf.m (body verbatim; the
% only copy). Used by the phi prior and importance-sampling ordinates of
% bvar.ml.mlvarsv_csv. The density counterpart of bvar.util.tnormrnd.
% Verified bitwise by tests/unit/test_mlvarsv_ml_densities.m.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function lden = ltnormpdf(x, mu, sig2, lb, ub)
c = -.5*log(2*pi*sig2) - log(normcdf((ub-mu)/sqrt(sig2))-normcdf((lb-mu)/sqrt(sig2)));
lden = c -.5/sig2*(x-mu).^2 ;
end
