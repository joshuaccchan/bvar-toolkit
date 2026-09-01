% This script evaluates log-likelihood of the BVAR-CSV-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

function llike = llike_CSV_MA(psi,U,Sig,h)
[T,n] = size(U);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T); 
Utld = Hpsi\U;
CSig = chol(Sig,'lower');
c = -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) - n/2*log(1+psi^2) -n/2*sum(h);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
llike = c -.5*(s2(1)/((1+psi^2)*exp(h(1))) + sum(s2(2:end)./exp(h(2:end))));
end
