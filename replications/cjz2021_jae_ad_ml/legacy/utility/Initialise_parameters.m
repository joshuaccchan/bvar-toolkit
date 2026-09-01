function [Sig,Omega,A,Ind_a,Ind,K_rt,K_rn,K_rr,I_r,I_rr,E_r,Ka,I_a,Iaa,Ea]=Initialise_parameters(Y,na,n,r,dim)
varY = var(Y)';
Sig.v = varY/2;Sig.d=zeros(n,dim);
Omega.v = mean(varY)/2*ones(r,1);Omega.d=zeros(r,dim);
a = rand(na,1);
A.v=[tril(ones(r),-1);ones(n-r,r)]';
Ind=find(A.v==1);
Ind_a=sparse(Ind,1:na,ones(na,1));
A.v(A.v==1)=a;
A.v(1:r,1:r)=A.v(1:r,1:r)+eye(r);
A.d=zeros(n*r,dim);

na = n*r - r*(r+1)/2;
Ka=commutation_matrix(na,na);
I_a=speye(na);
Iaa=speye(na^2);
Ea=elimination_matrix(na);

K_rn=commutation_matrix(r,n);
K_rr=commutation_matrix(r,r);
I_r=speye(r);
I_rr=speye(r^2);
E_r=elimination_matrix(r);
T=size(Y,1);
K_rt=commutation_matrix(r,T);
