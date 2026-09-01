function [L,D]= linvwishpdfAD(W, nu,S,n)

sum1 = sum(gammaln((nu+1-(1:n))/2));
const = nu*n/2*log(2) + n*(n-1)/4*log(pi) + sum1;
L= -const + nu/2 * log(S.det) - (nu + n + 1)/2 * log(W.det) - trace(S.v/W.v)/2;

D=0.5*nu/S.det*S.ddet-(nu + n + 1)/(2*W.det)*W.ddet-0.5*d_trace(d_prod(S.v,S.d,W.inv,W.invd),n);
end