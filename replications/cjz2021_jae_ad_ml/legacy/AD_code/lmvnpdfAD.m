function [L,D]= lmvnpdfAD(X,mu,sigma2,N)
error.v = X.v-mu.v;
error.d=X.d-mu.d;

kernel.v=-0.5*error.v'*sigma2.inv*error.v;
E=error.v*error.v';
kernel.d=-error.v'*sigma2.inv*error.d-0.5*E(:)'*sigma2.invd;
L = (-N/2)*log(2*pi) - 1/2*log(sigma2.det)+kernel.v;
D=-0.5/sigma2.det*sigma2.ddet+kernel.d;
end
