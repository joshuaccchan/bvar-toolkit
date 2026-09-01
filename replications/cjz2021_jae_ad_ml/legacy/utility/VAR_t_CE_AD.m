function MLCE=VAR_t_CE_AD(beta,Sig,T, nu, K)
global BigX X Y n v_0 S_0 B_0 b_0
   [q,k,M]=size(beta.d);
   [mbeta,cbeta]=PosteriorNormal(beta,k);
    Z=randn(q,M);
    S.v=mean(Sig.v,3)*T; 
    S.d=mean(Sig.d,3)*T;
    S=dmatrix(S,n,k);  L=GettheL(S,n);  
for i=1:M
    Beta.v=mbeta.v+cbeta.L*Z(:,i);
    Beta.d=mbeta.d+kron(Z(:,i)',eye(q))*cbeta.Ld;
    Sigma=d_Sigma2(S,L,T,n,k);
    B=reshape(Beta.v,n,K);
    err.v=Y-X*B';
    err.d=-BigX*d_trans(B,Beta.d);
    
    [l.v(i,1),l.d(i,:)]=logmvtpdfAD(err,Sigma,nu);
    
    [Pw.v,Pw.d]=linvwishpdfAD(Sigma,v_0,S_0,n); 
    [Pn.v,Pn.d]=lmvnpdfAD(Beta,b_0,B_0,q);
    p.v(i,1)=Pw.v+Pn.v; p.d(i,:)=Pw.d+Pn.d;
    
    [Pw.v,Pw.d]=linvwishpdfAD(Sigma,T,S,n); 
    [Pn.v,Pn.d]=lmvnpdfAD(Beta,mbeta,cbeta,q);
    CE.v(i,1)=Pw.v+Pn.v;
    CE.d(i,:)=Pw.d+Pn.d;
end
   [mlog, ind]=max(l.v+p.v-CE.v);
    md=l.d(ind,:)+p.d(ind,:)-CE.d(ind,:);
    E.v=l.v+p.v-CE.v-mlog;
    E.d=l.d+p.d-CE.d-repmat(md,M,1);
    E.v=exp(E.v);
    E.d=repmat(E.v,1,k).*E.d;
    E.v=mean(E.v);E.d=mean(E.d);
    MLCE.v=log(E.v)+mlog; 
    MLCE.d=1/E.v*E.d+md;
    