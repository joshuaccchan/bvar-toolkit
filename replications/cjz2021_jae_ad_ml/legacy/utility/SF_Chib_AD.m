function ml = SF_Chib_AD(lpa,store_a,store_Sig,...
    store_Omega,Y,prior, Get_F,Get_Omega,Get_Sig,nusig2,nuomega2)
    nsim_re = 10000; % size of reduced runs
    [~,na] = size(store_a.v);
    [T,n] = size(Y);
    r = size(store_Omega.v,2);    
    a_s.v = mean(store_a.v)';
    a_s.d=mean(store_a.d,3);% posterior ordinate to evaluate the ML identity   
    Sig_s.v = mean(store_Sig.v); Sig_s.v=Sig_s.v';
    Sig_s.d=mean(store_Sig.d,3);
    Omega_s.v= mean(store_Omega.v);Omega_s.v=Omega_s.v';
    Omega_s.d=mean(store_Omega.d,3);
       a_s.v=a_s.v';
   prior=prior(a_s,Sig_s,Omega_s);
    
    
    invSig=matrix_dot_division(ones(n,1),Sig_s);
    invOmega=matrix_dot_division(ones(r,1),Omega_s);
  A.v=[tril(ones(r),-1);ones(n-r,r)]';
Ind=find(A.v==1);
Ind_a=sparse(Ind,1:na,ones(na,1));
 A=Mstack(a_s,r,n,Ind,1:na);
     A.v(1:r,1:r)=A.v(1:r,1:r)+eye(r);
        % evaluate the log integrated likelihood at a_s,Sig_s,Omega_s
    llike = Intlike_AD(Y,a_s,invSig,invOmega,Ind,Ind_a,commutation_matrix(r,n),...
        commutation_matrix(n,n),elimination_matrix(n));      
        % obtain draws of (Sig,Omega) given a = a_s and Y using a reduced run
        % initialize the chain at (Sig_s,Omega_s)
    invS=invSig;
    invO=invOmega; 
    store_lpSig.v=zeros(nsim_re,1);store_lpSig.d=zeros(nsim_re,4);
       store_lpOmega.v=zeros(nsim_re,1);store_lpOmega.d=zeros(nsim_re,4);
   K_rt=commutation_matrix(r,T);
   ind_s=[0,0,1,0];ind_o=[0,0,0,1];
   In1=sparse(n*(0:n-1)+(1:n),1:n,ones(n,1));

In2=sparse(r*(0:r-1)+(1:r),1:r,ones(r,1));
for isim = 1:nsim_re
     
        invS.d =sparse(In1 * invS.d);
invS.v=sparse(1:n,1:n,invS.v);
  
invO.d=sparse(In2*invO.d);
invO.v=sparse(1:r,1:r,invO.v);
%         invS=d_diag(invS,true,n);
%         invO=d_diag(invO,true,r);
        
        
    F=Get_F(invS,invO,A);
    Ft=Mtrans(F,K_rt);
               % sample Sig
   [~, invS, S]=Get_Sig(A,Ft);   
   l=Ligamma(nusig2+T/2,S,Sig_s,n, ind_s);
   store_lpSig.v(isim) = l.v;store_lpSig.d(isim,:)=l.d;
    %Sample Omega
    [~,invO, S2]=Get_Omega(F);

   l=Ligamma(nuomega2+T/2,S2,Omega_s,r,ind_o);
      
     store_lpOmega.v(isim) = l.v;store_lpOmega.d(isim,:)=l.d;             
end
    
   lpSig=logmean(store_lpSig,nsim_re);
   lpOmega=logmean(store_lpOmega,nsim_re);
    ml.v = llike.v + prior.v - (lpa.v + lpSig.v + lpOmega.v);
    ml.d=llike.d+prior.d-(lpa.d+lpSig.d+lpOmega.d);
end