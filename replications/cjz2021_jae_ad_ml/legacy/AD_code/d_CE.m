function MLCE=d_CE(p,l,CE)
M=size(p.v,1);k=size(p.d,2);
[mlog, ind]=max(l.v+p.v-CE.v);
    md=l.d(ind,:)+p.d(ind,:)-CE.d(ind,:);
    E.v=l.v+p.v-CE.v-mlog;
    E.d=l.d+p.d-CE.d-repmat(md,M,1);
    E.v=exp(E.v);
    E.d=repmat(E.v,1,k).*E.d;
    E.v=mean(E.v);E.d=mean(E.d);
    MLCE.v=log(E.v)+mlog; 
    MLCE.d=1/E.v*E.d+md;