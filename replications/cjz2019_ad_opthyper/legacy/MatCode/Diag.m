% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function C=Diag(X,k)
s=k*(0:k-1)+(1:k);
In=sparse(s,1:k,ones(k,1));
C.d =sparse(In * X.d);
C.v=sparse(1:k,1:k,X.v);