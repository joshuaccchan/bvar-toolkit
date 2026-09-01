function [dXT] = d_trans(X, dX)

[r,c]=size(X);
K_nq=commutation_matrix(r,c);
      dXT = K_nq * dX;

