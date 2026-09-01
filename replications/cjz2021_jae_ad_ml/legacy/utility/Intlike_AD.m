function llike = Intlike_AD(Y,a,invSig,invOmega,Ind,Ind_a,K_rn,K_nn,E_n)
[T,n] = size(Y); 
r = size(invOmega.v,1);na=size(Ind_a,2);
A=Mstack(a,r,n,Ind,1:na);
A.v(1:r,1:r)=A.v(1:r,1:r)+eye(r);

In1=sparse(n*(0:n-1)+(1:n),1:n,ones(n,1));
In2=sparse(r*(0:r-1)+(1:r),1:r,ones(r,1));
 invSig.d =sparse(In1 * invSig.d);
invSig.v=sparse(1:n,1:n,invSig.v);
  
invOmega.d=sparse(In2*invOmega.d);
invOmega.v=sparse(1:r,1:r,invOmega.v);

% invSig=d_diag(invSig,true,n);
% invOmega=d_diag(invOmega,true,r);
AiSig = Mtimes(A, invSig);At=Mtrans(A,K_rn);
AiSigt=Mtrans(AiSig,K_rn);
    % iB = (A Omega A' + Sig)^{-1} obtained by the Woodbury formula
iB = Msubtraction(invSig, Mtimes(AiSigt,...
M_division(AiSig,Maddition(invOmega, Mtimes(AiSig,At)),true)));
CiB = Cholasky(iB,K_nn,E_n,speye(n),speye(n^2));
tmp=Mtimes(Y,CiB);
tmp.v=tmp.v(:);
tmp=dot_product(tmp,tmp);
CiB.v=diag(CiB.v);
CiB.d=In1'*CiB.d;
CiB=mlog(CiB);

llike.v=-0.5*tmp.v- T*n/2*log(2*pi)+ T*sum(CiB.v);
llike.d= T*sum(CiB.d) -.5*tmp.d;
end