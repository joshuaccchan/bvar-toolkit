% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of June 2018
% This function computes lower cholasky decomposition as well as its derivatives.


function [L] = Cholasky(A,K_nn,E_n,I_n,I_nn)
n=size(A.v,1);
if isempty(K_nn)
 I_n = speye(n);
  I_nn = speye(n^2);
  
  K_nn = commutation_matrix(n, n);
  E_n = elimination_matrix(n);
end
L.v=chol(A.v,'lower');   
  D_n = E_n';
  f1 = D_n * ((E_n * (I_nn + K_nn) * kron(L.v, I_n) * D_n) \ E_n);
  L.d= f1 * A.d;