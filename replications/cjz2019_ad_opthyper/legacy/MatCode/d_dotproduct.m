% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function dC=d_dotproduct(A,dA, B,dB)

Da=diag(A);
Db=diag(B);

dC =Da*dB +Db*dA;
