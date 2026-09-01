% This function evaluates the log determinant
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function k = ldet(Omega)
    k = 2*sum(log(diag(chol(Omega,'lower'))));
end