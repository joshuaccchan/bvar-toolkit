function prior=prior_AD_ordinates(ax,s,o,a0,Va,nusig2,Ssig2,nuomega2,Somega2)
[R,na]=size(ax.v);
d=size(ax.d,2);
n=size(s.v,1);
r=size(o.v,1);
err.v=ax.v-a0;
err.d=ax.d;
err.d(:,1)=err.d(:,1)-1;
err=Sum_matrix(matrix_dot_times(err,err),2);
err.d=0.5/Va*err.d-0.5/Va^2*err.v*[0,1,sparse(1,d-2)];
err.v=0.5/Va*err.v;
prior.v=-na/2*log(2*pi*Va)-err.v;
prior.d=-na/(2*Va)*ones(R,1)*[0,1,sparse(1,d-2)]-err.d;

L=mlog(s);
L2=matrix_dot_division(ones(n,R),s);
l_sig.v=nusig2*log(nusig2-1)-gammaln(nusig2)...
    -(nusig2+1)*L.v-(nusig2-1)*L2.v;
v=log(nusig2-1)+nusig2/(nusig2-1)-psi(nusig2)-L.v(:)-L2.v(:);
l_sig.d=v*[0,0,1,sparse(1,d-3)]-(nusig2+1)*L.d-(nusig2-1)*L2.d;

l_sig=Sum_matrix(l_sig,1);l_sig.v=l_sig.v';

L=mlog(o);
L2=matrix_dot_division(ones(r,R),o);
l_omega.v=nuomega2*log(nuomega2-1)-gammaln(nuomega2)...
    -(nuomega2+1)*L.v-(nuomega2-1)*L2.v;
v=log(nuomega2-1)+nuomega2/(nuomega2-1)-psi(nuomega2)-L.v(:)-L2.v(:);
l_omega.d=v*[0,0,0,1,sparse(1,d-4)]-(nuomega2+1)*L.d-(nuomega2-1)*L2.d;
if size(l_omega.v,1)>1
l_omega=Sum_matrix(l_omega,1);
end
l_omega.v=l_omega.v';

prior.v=prior.v+l_sig.v+l_omega.v;
prior.d=prior.d+l_sig.d+l_omega.d;
