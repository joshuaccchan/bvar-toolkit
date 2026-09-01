function [l,dl]=logmvtpdfAD(X,C, d)
[T,n]=size(X.v);
k=size(C.d,2);
nu=0.5*d;
newnu=0.5*(d+n);

l=T*log(gamma(n/2)/beta(nu,n/2))-T*n/2*log(d*pi)-T/2*log(C.det);
  
dl=-T/(2*C.det)*C.ddet;

%V= -newnu*sum(log(1+1/d*sum(X.v*C.inv.*X.v,2)));

V.v=X.v*C.inv;V.d=d_prod(X.v,X.d,C.inv,C.invd);
    gam.v=sum(V.v.*X.v,2);
    gam.d=d_diagg(d_prod(V.v,V.d,X.v',d_trans(X.v,X.d)),T);
    
gam.d=gam.d./d;
gam.v=gam.v/d+1;
gam.d=repmat(1./gam.v,1,k).*gam.d;
gam.v=log(gam.v);
dl=dl-newnu*sum(gam.d);
l=l-newnu*sum(gam.v);