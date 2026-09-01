function [v,d]=logGammaAD(S,dS,a,da,b,db,ind)
l=log(S);
if ind==1
v=-S'*b.^(-1)+sum(a.*(l-log(b))-l)-length(S)*logMvGamma(a,1);
else
  v=sum(log(gampdf(S,a,b)));
end
d=((a-1)./S-1./b)'*dS+sum(l-log(b)-psi(a))*da-(a./b-S./b.^2)'*db;
%d=d';