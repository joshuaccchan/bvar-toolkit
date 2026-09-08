% Body from chan2023_joe_mlvarsv/legacy/utility/mgammaln.m (identical modulo comments).
% Equivalence: tests/unit/test_util_small.m. Record: tests/variant_map.md.
% This function evaluates the multivariate gamma function (in log)
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function k = mgammaln(n,x)    
    k = n*(n-1)/4*log(pi) + sum(gammaln((x+(0:-.5:(1-n)/2))));
end