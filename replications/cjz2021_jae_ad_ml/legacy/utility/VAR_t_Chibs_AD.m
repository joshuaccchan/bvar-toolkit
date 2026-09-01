function ML=VAR_t_Chibs_AD(beta,Sig,b,Vb,T,nu, K)
global BigX X Y n v_0 S_0 B_0 b_0
newnu=(nu+n)/2;
[q,k,M]=size(beta.d);
Bm.v=mean(beta.v)';Bm.d=mean(beta.d,3);
Sm.v=mean(Sig.v,3);  Sm.d=mean(Sig.d,3); Sm=dmatrix(Sm,n,k);
B=reshape(Bm.v,n,K);
error.v=Y-X*B';error.d=-BigX*d_trans(B,Bm.d);
[l.v,l.d]=logmvtpdfAD(error,Sm,nu);
[Pw.v,Pw.d]=linvwishpdfAD(Sm,v_0,S_0,n); 
[Pn.v,Pn.d]=lmvnpdfAD(Bm,b_0,B_0,q);
prior.v=Pw.v+Pn.v; prior.d=Pw.d+Pn.d;

G=gamrnd(newnu,1,T,M+500);

Sigma=Sm;
Burn=100;
for i=1:M+Burn
    %sample lambda
    L=Getlambda(G(:,i),error,Sigma,nu,k,T);
    %sample Sigma
    C.v=L.v*error.v;
    C.d=d_prod(L.v,L.d,error.v,error.d);
    delta.v=S_0.v+error.v'*C.v;
    delta.d=S_0.d+d_prod(error.v',d_trans(error.v,error.d),C.v,C.d);
    delta.v=delta.v\speye(n);
    delta.d=d_minverse(delta.v,delta.d);
    Sigma=d_Sigma(delta,n,v_0+T);
    if i>Burn
        ind=i-Burn;
        mu.v=b.v(ind,:)';mu.d=b.d(:,:,ind);
        Cov.v=Vb.v(:,:,ind);Cov.d=Vb.d(:,:,ind);Cov=dmatrix(Cov,q,k);
        [p_b.v(ind,1),p_b.d(ind,:)]=lmvnpdfAD(Bm,mu,Cov,q);
        delta.v=delta.v\speye(n);
        delta.d=d_minverse(delta.v,delta.d);
        delta=dmatrix(delta,n,k);
        [p_s.v(ind,1),p_s.d(ind,:)]=linvwishpdfAD(Sm,v_0+T,delta,n);
    end
end  
p_b=mdensity(p_b,k);

p_s=mdensity(p_s,k);

ML.v=l.v+prior.v-p_s.v-p_b.v;
ML.d=l.d+prior.d-p_s.d-p_b.d;
end

