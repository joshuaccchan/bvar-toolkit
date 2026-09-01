function L=Getlambda(G,error,Sigma,nu,k,T)
C.v=error.v*Sigma.inv;C.d=d_prod(error.v,error.d,Sigma.inv,Sigma.invd);
    gam.v=0.5*nu+0.5*sum(C.v.*error.v,2);
    gam.d=0.5*d_diagg(d_prod(C.v,C.d,error.v',d_trans(error.v,error.d)),T);
    invlambda.v=G./gam.v;invlambda.d=-repmat(invlambda.v./gam.v,1,k).*gam.d;
    %lambda.v=1./invlambda.v;lambda.d=repmat(-lambda.v.^2,1,k).*invlambda.d;
    L.v=sparse(1:T,1:T,invlambda.v);
    L.d=d_diag(invlambda.d,T);