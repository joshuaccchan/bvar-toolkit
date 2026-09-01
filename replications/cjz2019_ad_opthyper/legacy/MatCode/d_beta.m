% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function [beta,bg, error]=d_beta(Vg,bg,Z,d)
global X Y BigX K n
 bg.d=d_prod(Vg.v,Vg.d,bg.v,bg.d);
   bg.v=Vg.v*bg.v;

   C=chol(Vg.v,'lower');
   beta.v=bg.v+C*Z;
   beta.d=bg.d+kron(Z',eye(d))*d_cholasky(C,Vg.d);
    b=reshape(beta.v,K,n);
   error.v=(Y-X*b); 
   error.d=-BigX*beta.d;
   