% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function [C] = Mtimes(A, B)

  I_a = sparse(eye(size(A.v,1)));
  I_b =sparse( eye(size(B.v,2)));
  f1 =sparse( kron(I_b, A.v));
  f2 = sparse(kron(B.v', I_a));
  C.d = f1 * B.d + f2 * A.d;
  C.v=A.v*B.v;
