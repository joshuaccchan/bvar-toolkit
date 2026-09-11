% bvar.util.ldet - log determinant of a symmetric positive definite matrix, taken
% from its Cholesky factor.
%
%   k = bvar.util.ldet(Omega)
%
%   Omega : symmetric positive definite matrix; chol errors if it is not
%   k     : log(det(Omega))
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function k = ldet(Omega)
    k = 2*sum(log(diag(chol(Omega,'lower'))));
end