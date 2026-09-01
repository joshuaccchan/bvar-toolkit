function ML=VAR_Chib_AD(beta,sigma, Varg,Bg)
global  nu_1  n K B_0 S_0 nu_0 T  Y X BigX b_0
   [M,q]=size(beta.v);
    b.v=reshape(mean(beta.v),K,n);
    b.d=mean(beta.d,3);
    k=size(b.d,2);
    sigm.v=mean(sigma.v,3);
    sigm.d=mean(sigma.d,3);
    sigm=dmatrix(sigm,n,k);
    error.v=(Y-X*b.v); 
    error.d=-BigX*b.d;
    W.v=error.v'*error.v; W.d=d_prod(error.v',d_trans(error.v,error.d),error.v,error.d);W=dmatrix(W,n,k);
    S.v=S_0.v+W.v;S.d=S_0.d+ W.d;S=dmatrix(S,n,k);
    b.v=reshape(b.v,q,1);
    
    log_B_Sigma.v=0; 
    log_B_Sigma.d=zeros(1,k);
    for i=1:M
        mu.v=Bg.v(:,i);
        mu.d=Bg.d(:,:,i);
        V.v=Varg.v(:,:,i);V.d=Varg.d(:,:,i);V=dmatrix(V,q,k);
       [Bn.v,Bn.d]=lmvnpdfAD(b,mu, V,q);  
       Bn.v=exp(Bn.v); Bn.d=repmat(Bn.v,1,k).*Bn.d;
       log_B_Sigma.v=log_B_Sigma.v+Bn.v;
       log_B_Sigma.d=log_B_Sigma.d+Bn.d;
    end
    log_B_Sigma.v=log_B_Sigma.v/M;
   log_B_Sigma.d=log_B_Sigma.d./(M*log_B_Sigma.v);
   log_B_Sigma.v=log(log_B_Sigma.v);
    l.v= -T*n/2*log(2*pi) -T/2*log(sigm.det) - 1/2*trace(sigm.v\W.v);
    l.d=T/(2*sigm.det)*d_det(sigm.det,sigm.v,sigm.d,n,k)-0.5*d_trace(d_prod(sigm.inv,sigm.invd,W.v,W.d),n);

    [Pw.v,Pw.d]=linvwishpdfAD(sigm,nu_0,S_0,n); 
    [Pn.v,Pn.d]=lmvnpdfAD(b,b_0,B_0,q);
    prior.v=Pw.v+Pn.v; prior.d=Pw.d+Pn.d;
    [Pow.v,Pow.d]=linvwishpdfAD(sigm,nu_1,S,n);
    post.v=Pow.v+log_B_Sigma.v; 
    post.d=Pow.d+log_B_Sigma.d; 
    ML.v=l.v+prior.v-post.v;
     ML.d=l.d+prior.d-post.d;
    