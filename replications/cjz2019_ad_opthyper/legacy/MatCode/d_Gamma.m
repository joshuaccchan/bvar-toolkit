% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function [ G] = d_Gamma(alpha,delta,dim, dim1) 
  g=gamrnd(alpha,1);
  fun = @(t) log(t) .*t.^(alpha-1).*exp(-t);
  num_1 = integral(fun, 0, g)/gamma(alpha);
  num_2 = psi(alpha) *gamcdf(g,alpha,1);
  G.d =[sparse(1,dim-1), - (num_1 - num_2) ./ gampdf(g, alpha, 1),sparse(dim1-dim)];  
  G.v=g/delta.v;
  G.d=1/delta.v*G.d-g/delta.v^2*delta.d;
 
  
  

  