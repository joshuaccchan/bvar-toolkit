% This function evaluates the truncated normal pdf (in log)
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function lden = ltnormpdf(x, mu, sig2, lb, ub)
c = -.5*log(2*pi*sig2) - log(normcdf((ub-mu)/sqrt(sig2))-normcdf((lb-mu)/sqrt(sig2)));
lden = c -.5/sig2*(x-mu).^2 ;
end