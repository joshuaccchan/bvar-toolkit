function [dC] = d_prod(A, dA, B, dB)
  I_a = sparse(eye(size(A,1)));
  I_b =sparse( eye(size(B,2)));
  f1 =sparse( kron(I_b, A));
  f2 = sparse(kron(B', I_a));
  dC = f1 * dB + f2 * dA;
