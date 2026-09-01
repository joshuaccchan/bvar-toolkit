% This function computes the gradient of the log multivariate gamma
% function with respect to kappas
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function k=mgammaln_AD(n,x)
 k.v = n*(n-1)/4*log(pi) + sum(gammaln((x.v+(0:-.5:(1-n)/2))));
 k.d=sum(psi( x.v+(0:-.5:(1-n)/2)))*x.d;
 