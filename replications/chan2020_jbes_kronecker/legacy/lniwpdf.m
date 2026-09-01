% This script returns the log density of the normal-inverse-Wishart
% distribuiton
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

function lden = lniwpdf(A,Sig,A0,iVA0,nu0,S0)
[k,n] = size(A);
CSig = chol(Sig,'lower');
cA = -n*k/2*log(2*pi) + n*sum(log(diag(chol(iVA0))));
cSig = -nu0*n/2*log(2) -n*(n-1)/4*log(pi) -sum(gammaln((nu0+1-(1:n))/2))...
    + nu0*sum(log(diag(chol(S0))));
tmp = A-A0;
lden = cA + cSig - (n+nu0+k+1)*sum(log(diag(CSig))) ...
    - .5*trace(Sig\(S0+tmp'*iVA0*tmp));
end