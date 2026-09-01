% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function Sigma=d_Sigma(delta, n, nu_1) 
   Sigma.inv=wishrnd(delta.v,nu_1);
   L.v=chol(delta.v)';
   L.d=d_cholasky(L.v,delta.d);
   A=(L.v\speye(n))*Sigma.inv;
   A2=Sigma.inv*(L.v'\speye(n));
   
   Sigma.invd=kron(A',eye(n))*L.d+kron(eye(n),A2)*d_trans(L.v,L.d);
   
   Sigma.v=Sigma.inv\speye(n);  
   Sigma.d=d_minverse(Sigma.v,Sigma.invd);