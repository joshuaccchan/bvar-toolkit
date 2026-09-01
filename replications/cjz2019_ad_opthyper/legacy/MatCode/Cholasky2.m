% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019

function [L] = Cholasky2(A,K_nn,E_n,I_n,I_nn)
  L.v=chol(A.v,'lower');   
  D_n = E_n';
  f1 = D_n * ((E_n * (I_nn + K_nn) * kron(L.v, I_n) * D_n) \ E_n);
  L.d= f1 * A.d;