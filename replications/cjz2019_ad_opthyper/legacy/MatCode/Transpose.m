% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function Xt=Transpose(X,r,c)
K_nq=commutation_matrix(r,c);
     Xt.d = K_nq * X.d;
     Xt.v=X.v';