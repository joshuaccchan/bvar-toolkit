% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function S=dmatrix(S,n,k)
S.inv=S.v\speye(n);
S.invd=d_minverse(S.inv,S.d);
S.det=det(S.v);
S.ddet=d_det(S.det,S.v,S.d,n,k);