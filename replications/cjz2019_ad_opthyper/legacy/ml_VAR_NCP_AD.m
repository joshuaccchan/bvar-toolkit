% This function computes the gradient of the marginal likelihood under the 
% natural conjugate prior with respect to kappas
%
% See:
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019.

function lml = ml_VAR_NCP_AD(VA,S0,nu0,KA,S_hat,T,n,three)
    k=size(VA.v,1);
    lml.v = -n*T/2*log(pi);
    if three==false
    lml.d=zeros(1,5);
    else
        lml.d=zeros(1,3);
    end   
    V=Maddition(Ldet_AD(VA,k),Ldet_AD(KA,k));
    V.v=V.v*(-n/2);V.d=V.d*(-n/2);
    lml=Maddition(lml,V);
   if three==false
       nu1.v=nu0/2; nu1.d=0.5*[zeros(1,3),1,0];
       nu2.v=(nu0+T)/2; nu2.d=0.5*[zeros(1,3),1,0];
      lml=Maddition(lml,Mtimes(nu1,Ldet_AD(S0,n)));
      lml=Msubtraction(lml,Mtimes(nu2,Ldet_AD(S_hat,n)));
      lml=Msubtraction(Maddition(lml,mgammaln_AD(n,nu2)),mgammaln_AD(n,nu1));

   else
       lml.v=lml.v+ nu0/2*ldet(S0);
       V=Ldet_AD(S_hat,n);
       lml.v=lml.v-(nu0+T)/2*V.v;
       lml.d=lml.d-(nu0+T)/2*V.d;
       lml.v=lml.v + mgammaln(n,(nu0+T)/2) - mgammaln(n,nu0/2);
   end    