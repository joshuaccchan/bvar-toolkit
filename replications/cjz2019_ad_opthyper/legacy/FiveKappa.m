% This function evaluates the log margina likelihood and its gradient with
% respect to kappa1 - kappa5
%
% See:
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019.
function [ML, dML]=FiveKappa(kappa,Y,Y0,T,n,k,p)
global K_kn
%#of coefficients in each equation
%prior
[A0,VA,nu0,S0] = prior_NCP_AD(p,kappa,Y0,Y,false);
S0.inv=S0.v\speye(n);
iVA=Minverse(VA,k);VA.inv=iVA.v;
% compute the parameters for the posterior density
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for i=1:p
    Z(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
Z = [ones(T,1) Z];
KA.v = iVA.v + Z'*Z;
KA.d=iVA.d;
Ahat=Mtimes(iVA,A0);
Ahat.v=Ahat.v+Z'*Y;
iKA=Minverse(KA,k);
KA.inv=iKA.v;
Ahat=Mtimes(iKA,Ahat);
S_hat=Maddition(S0,Msubtraction(Mtimes(Transpose2(A0,K_kn),Mtimes(iVA,A0)),Mtimes(Transpose2(Ahat,K_kn),Mtimes(KA,Ahat))));
S_hat.v = S_hat.v + Y'*Y ;
S_hat=Maddition(S_hat,Transpose(S_hat,n,n));
S_hat.v=S_hat.v*0.5;S_hat.d=S_hat.d*0.5;
S_hat.inv=S_hat.v\speye(n);
result=ml_VAR_NCP_AD(VA,S0,nu0,KA,S_hat,T,n,false);
ML=-result.v;

dML=-result.d';
