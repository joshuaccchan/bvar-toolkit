function [dC] = d_kron(A, dA, B, dB)
  [m,n] = size(A);
  [p,q] = size(B);
  I_n = eye(n);
  K_qm = commutation_matrix(q, m);
  I_p = eye(p);
  f1 = sparse(kron(kron(I_n, K_qm), I_p));
      dC = f1 * (sparse(kron(A(:), dB)) + sparse(kron(dA, B(:))));

