function [A,a,Ft,K_a,mu_a, L]=Sample_A(F,invSig,Y,idx, Va,a0,na,Ind,Ind_a,K_rt,K_rn,...
    Ka,Ea,Ia,Iaa)
[~,n]=size(Y);
r=size(F.v,1);
 Ft=Mtrans(F,K_rt);
  BigF=d_kron2(invSig,Mtimes(F,Ft),K_rn);
  K_a=Mtimes(Ind_a',Mtimes(BigF,Ind_a));
  K_a.v=K_a.v+Va^(-1)*speye(na);
  K_a.d(idx,2)= K_a.d(idx,2)-Va^(-2);
  mu=Mtimes(Mtimes(F,Y),invSig);
  mu.v=mu.v(:);mu_a=Mtimes(Ind_a',mu);
  mu_a.v=mu_a.v+a0/Va*ones(na,1);
  mu_a.d(:,1)=mu_a.d(:,1)+1/Va;
  mu_a.d(:,2)=mu_a.d(:,2)-a0/Va^2;
  mu_a=M_division(mu_a,K_a,true);
  L=Cholasky(K_a,Ka,Ea,Ia,Iaa);
  a=Maddition(mu_a,M_division(randn(na,1),Mtrans(L,Ka),true));
  A=Mstack(a,r,n,Ind,1:na);
  A.v(1:r,1:r)=A.v(1:r,1:r)+speye(r);