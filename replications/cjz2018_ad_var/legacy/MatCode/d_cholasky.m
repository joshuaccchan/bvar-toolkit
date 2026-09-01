function [dL] = d_cholasky(L, dA)
  %  Note that LL' = A or L = chol(A)'
  n = size(L,1);
  I_n = eye(n);
  I_nn = eye(n^2);
  K_nn = commutation_matrix(n, n);
  E_n = elimination_matrix(n);
  D_n = E_n';
 f1 = D_n * inv(E_n * (I_nn + K_nn) * kron(L, I_n) * D_n) * E_n;
 dL= f1 * dA;

