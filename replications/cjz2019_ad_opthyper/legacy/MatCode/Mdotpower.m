% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function A=Mdotpower(C,n)

A.v=C.v.^n;
D=n*C.v.^(n-1);
A.d=diag(D(:))*C.d;

