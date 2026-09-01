function L= lmvnpdfAD2(X,mu,K, L, N)
err=Msubtraction(X,mu);
Kernel=dot_product(Mtimes(K,err),err);

Kernel.v=-0.5*Kernel.v;
Kernel.d=-0.5*Kernel.d;
n=size(L.v,1);
In1=sparse(n*(0:n-1)+(1:n),1:n,ones(n,1));
L.d =sparse(In1' *L.d);
L.v=diag(L.v);
L.d=sparse(1:n,1:n,1./L.v)*L.d;L.v=log(L.v);
%L=mlog(L);
L.v=sum(L.v);
L.d=sum(L.d);
L.v = (-N/2)*log(2*pi)+L.v+Kernel.v;
L.d=L.d+Kernel.d;
end
