% This function evaluates the inverse-gamma pdf (in log)
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function lden = ligampdf(x,a,b)
    lden = a.*log(b) - gammaln(a) - (a+1).*log(x) - b./x;
end