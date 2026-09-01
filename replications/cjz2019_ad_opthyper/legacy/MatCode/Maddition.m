% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function C=Maddition(A,B)

C.v=A.v+B.v;
C.d=A.d+B.d;