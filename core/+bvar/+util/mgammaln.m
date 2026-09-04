% Extracted 2026-09-01 (step 3, zero-risk core): verified identical (modulo comments/whitespace). Canonical source: chan2023_joe_mlvarsv/legacy/utility/mgammaln.m (also code-identical in cjz2019_ad_opthyper/legacy).
% This function evaluates the multivariate gamma function (in log)
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function k = mgammaln(n,x)    
    k = n*(n-1)/4*log(pi) + sum(gammaln((x+(0:-.5:(1-n)/2))));
end