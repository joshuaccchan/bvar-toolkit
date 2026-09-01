
function [invSig_IS,invOmega_IS,mA]=Sample_IS(a_bar,CDa_bar,...
    nusig2_bar,Ssig2_bar,...
    nuomega2_bar,Somega2_bar,R)
n=size(nusig2_bar.v,1);r=size(nuomega2_bar.v,1);
G.v=zeros(R,n);
G.d=zeros(R*n,1);
for i=1:n
G.v(:,i)=gamrnd(nusig2_bar.v(i),1,R,1);
dif=(gamcdf(G.v(:,i),nusig2_bar.v(i)+1e-8,1)-gamcdf(G.v(:,i),nusig2_bar.v(i)-1e-8,1))./2e-8;
pdf=gampdf(G.v(:,i),nusig2_bar.v(i));
G.d((i-1)*R+1:i*R)=-dif./pdf;
end
invSig_IS=Mtrans(G,[]);
invSig_IS.d=sparse(1:R*n,1:R*n,invSig_IS.d)*repmat(nusig2_bar.d,R,1);
invSig_IS=Mtimes(d_diag2(Ssig2_bar,true,n),invSig_IS);

G.v=zeros(R,r);
G.d=zeros(R*r,1);
for i=1:r
G.v(:,i)=gamrnd(nuomega2_bar.v(i),1,R,1);
dif=(gamcdf(G.v(:,i),nuomega2_bar.v(i)+1e-8,1)-gamcdf(G.v(:,i),nuomega2_bar.v(i)-1e-8,1))./2e-8;
pdf=gampdf(G.v(:,i),nuomega2_bar.v(i));
G.d((i-1)*R+1:i*R)=-dif./pdf;
end
invOmega_IS=Mtrans(G,[]);
invOmega_IS.d=sparse(1:R*r,1:R*r,invOmega_IS.d)*repmat(nuomega2_bar.d,R,1);
invOmega_IS=Mtimes(d_diag2(Somega2_bar,true,r),invOmega_IS);

% invOmega_IS=Gamma(nuomega2_bar.v,R);
% invOmega_IS.d=sparse(1:R*r,1:R*r,invOmega_IS.d)*repmat(nuomega2_bar.d,R,1);
% invOmega_IS=Mtimes(d_diag(Somega2_bar,true,r),invOmega_IS);

na=size(a_bar.v,1);
Z=randn(R,na);
mA.v=repmat(a_bar.v',R,1)+Z*CDa_bar.v;
mA.d=kron(speye(na),ones(R,1))*a_bar.d+kron(speye(na),Z)*CDa_bar.d;
