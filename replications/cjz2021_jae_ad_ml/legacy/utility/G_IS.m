function g_IS =G_IS(ax,s,o,a_bar,CDa_bar,nusig2_bar,Ssig2_bar,nuomega2_bar,Somega2_bar)

[R,na]=size(ax.v);n=size(s.v,1);
r=size(o.v,1);
err.v=ax.v-a_bar.v';
err.d=ax.d-kron(speye(na),ones(R,1))*a_bar.d;
err=M_division(err,CDa_bar,false);
C=Sum_matrix(matrix_dot_times(err,err),2);
L=d_diag2(CDa_bar,false,na);
L=mlog(L);L.v=sum(L.v);L.d=sum(L.d);
g_IS.v=-0.5*C.v-  na/2*log(2*pi)-L.v;
g_IS.d=-0.5*C.d-L.d;

L=mlog(s);
L2=matrix_dot_division(ones(n,R),s);
l_sig= matrix_dot_times(nusig2_bar,mlog(Ssig2_bar));
l_sig.v=l_sig.v-gammaln(nusig2_bar.v);
l_sig.d=l_sig.d-sparse(1:n,1:n,psi(nusig2_bar.v))*nusig2_bar.d;
l_sig.v=l_sig.v*ones(1,R);l_sig.d=repmat(l_sig.d,R,1);
l_sig=Msubtraction(l_sig,Mtimes(d_diag2(Ssig2_bar,true,n),L2));
nusig2_bar.v=nusig2_bar.v+1;
l_sig=Msubtraction(l_sig,Mtimes(d_diag2(nusig2_bar,true,n),L));
l_sig=Sum_matrix(l_sig,1);

L=mlog(o);
L2=matrix_dot_division(ones(r,R),o);
l_omega= matrix_dot_times(nuomega2_bar,mlog(Somega2_bar));
l_omega.v=l_omega.v-gammaln(nuomega2_bar.v);
l_omega.d=l_omega.d-sparse(1:r,1:r,psi(nuomega2_bar.v))*nuomega2_bar.d;
l_omega.v=l_omega.v*ones(1,R);l_omega.d=repmat(l_omega.d,R,1);
l_omega=Msubtraction(l_omega,Mtimes(d_diag2(Somega2_bar,true,r),L2));
nuomega2_bar.v=nuomega2_bar.v+1;
l_omega=Msubtraction(l_omega,...
    Mtimes(d_diag2(nuomega2_bar,true,r),L));
if size(l_omega.v,1)>1
l_omega=Sum_matrix(l_omega,1); 
end

g_IS.v=g_IS.v+l_sig.v'+l_omega.v';
g_IS.d=g_IS.d+l_sig.d+l_omega.d;