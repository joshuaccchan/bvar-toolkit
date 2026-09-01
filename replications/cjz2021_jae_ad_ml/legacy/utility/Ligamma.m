
function l=Ligamma(nu,S,X,n, ind)

l.v=nu*log(S.v) - gammaln(nu) - (nu+1)*log(X.v) - S.v./X.v;
l.d=(log(S.v)-psi(nu)-log(X.v))*ind;
l.d=l.d+sparse(1:n,1:n,nu./S.v-1./X.v)*S.d;
l.d=l.d+sparse(1:n,1:n,-(nu+1)./X.v+S.v./X.v.^2)*X.d;
l.v=sum(l.v);
l.d=sum(l.d);

