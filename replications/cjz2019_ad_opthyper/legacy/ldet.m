% This function evaluates the log determinant
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function k = ldet(Omega)
    k = 2*sum(log(diag(chol(Omega,'lower'))));
end