function m=mlog(S)
[n,s]=size(S.v);

m.d=sparse(1:n*s,1:n*s,1./S.v(:))*S.d;
m.v=log(S.v);