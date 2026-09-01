function [dinv_A] = d_minverse(inv_A, dA)
  f1 = sparse(kron(inv_A', inv_A));
 dinv_A = -f1 * dA;
