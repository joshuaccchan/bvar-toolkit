% bvar.util.mgammaln - log of the multivariate gamma function Gamma_n(x).
%
%   k = bvar.util.mgammaln(n, x)
%
%   n : dimension, a positive integer
%   x : scalar argument; needs x > (n-1)/2, else gammaln returns Inf/NaN
%   k : log Gamma_n(x)
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function k = mgammaln(n,x)    
    k = n*(n-1)/4*log(pi) + sum(gammaln((x+(0:-.5:(1-n)/2))));
end