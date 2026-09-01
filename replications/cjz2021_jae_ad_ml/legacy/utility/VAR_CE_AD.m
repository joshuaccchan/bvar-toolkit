function MLCE=VAR_CE_AD(beta,S)
global K n   B_0 S_0 nu_0 T  X Y BigX b_0
  [q,k,M]=size(beta.d); 
   [mbeta,cbeta]=PosteriorNormal(beta,k);
    Z=randn(q,M);
    S.v=mean(S.v,3)*T; 
    S.d=mean(S.d,3)*T;
    S=dmatrix(S,n,k);
  L=GettheL(S,n);
    for i=1:M
   Sigma=d_Sigma2(S,L,T,n,k);    
    B.v=mbeta.v+cbeta.L*Z(:,i); 
    B.d=mbeta.d+kron(Z(:,i)',eye(q))*cbeta.Ld;
    error.v=(Y-X*reshape(B.v,K,n)); 
    error.d=-BigX*B.d;
    W.d=d_prod(error.v',d_trans(error.v,error.d),error.v,error.d);
    W.v=error.v'*error.v;
    
    l.v(i,1)= -T*n/2*log(2*pi) -T/2*log(Sigma.det) - 1/2*trace(Sigma.inv*W.v); 
    l.d(i,:)=-T/(2*Sigma.det)*Sigma.ddet-0.5*d_trace(d_prod(Sigma.inv,Sigma.invd,W.v,W.d),n);
   
    [Pw.v,Pw.d]=linvwishpdfAD(Sigma,nu_0,S_0,n); 
    [Pn.v,Pn.d]=lmvnpdfAD(B,b_0,B_0,q);
    p.v(i,1)=Pw.v+Pn.v; p.d(i,:)=Pw.d+Pn.d;
    [Pw.v,Pw.d]=linvwishpdfAD(Sigma,T,S,n); 
    [Pn.v,Pn.d]=lmvnpdfAD(B,mbeta,cbeta,q);
    CE.v(i,1)=Pw.v+Pn.v;
    CE.d(i,:)=Pw.d+Pn.d;
    end
    MLCE=d_CE(p,l,CE);
end    
    