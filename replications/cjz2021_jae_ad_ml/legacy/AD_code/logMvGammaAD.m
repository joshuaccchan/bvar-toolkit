function [v,d]=logMvGammaAD(S,dS,a,da,b,db,k)
dS=dS';
l=log(S);
%v=(a-1)*l-S./b-a*log(b)-log(gamma(a));
v=log(gampdf(S,a,b));
d=repmat((a-1)./S-1./b,1,k).*dS+(l-log(b)-psi(a))*da-(a./b-S./b.^2)*db;
d=d';
end





