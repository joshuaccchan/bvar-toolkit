function [L,D]= linvwishpdfAd2(W, nu,S,n, ind,v)

sum1 = sum(gammaln((nu+1-(1:n))/2));
const=n/2*log(2);
dv=const+sum(psi((nu+1-(1:n))/2))*0.5;
const = nu*const + n*(n-1)/4*log(pi) + sum1;

l1=log(S.det);
l2=log(W.det);
L= -const + nu/2 * l1 - (nu + n + 1)/2 * l2 - trace(S.v/W.v)/2;
dv=-dv+0.5*l1-0.5*l2;

D=0.5*nu/S.det*S.ddet-(nu + n + 1)/(2*W.det)*W.ddet-0.5*d_trace(d_prod(S.v,S.d,W.inv,W.invd),n);
D(1,ind)=dv*v;
end