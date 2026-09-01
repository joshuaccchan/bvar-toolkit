% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function Sigma=d_Sigma2(S,L,T,n,k)
   Sigma.inv=wishrnd(S.inv,T);
    A=L.inv*Sigma.inv;
    A2=Sigma.inv*L.tinv;
    Sigma.invd=kron(A',eye(n))*L.d+kron(eye(n),A2)*L.dt;
    Sigma.v=Sigma.inv\speye(n);
    Sigma.d=d_minverse(Sigma.v,Sigma.invd);
    Sigma.det=det(Sigma.v);Sigma.ddet=d_det(Sigma.det,Sigma.v,Sigma.d,n,k);   
end

