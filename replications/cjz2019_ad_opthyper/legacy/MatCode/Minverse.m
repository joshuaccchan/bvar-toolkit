% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function [inv_A] = Minverse(A, dim)
inv_A.v=A.v\speye(dim);
f1 = sparse(kron(inv_A.v', inv_A.v));
 inv_A.d = -f1 * A.d;
