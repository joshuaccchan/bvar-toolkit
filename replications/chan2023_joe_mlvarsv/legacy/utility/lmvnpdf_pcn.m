% This function evaluates the normal pdf (in log) 
% the input takes the precision matrix instead of the covariance matrix
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function lden = lmvnpdf_pcn(x,mu,K)
    n = length(mu);
    CK = chol(K,'lower');
    e = CK'*(x-mu); 
    lden = -n/2*log(2*pi) + sum(log(diag(CK))) - .5*(e'*e);
end