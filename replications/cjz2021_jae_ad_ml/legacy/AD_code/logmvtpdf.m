function l=logmvtpdf(X,C, d)
[T,n]=size(X);

nu=0.5*d;
newnu=0.5*(d+n);
invS=C\speye(n);
l=T*log(gamma(n/2)/beta(nu,n/2))-T*n/2*log(d*pi)-T/2*log(det(C))-newnu*sum(log(1+1/d*sum(X*invS.*X,2)));