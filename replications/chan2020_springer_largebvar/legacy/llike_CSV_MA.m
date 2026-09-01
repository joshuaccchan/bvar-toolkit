% This script evaluates the conditional likelihood of BVAR-CSV-MA
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

function llike = llike_CSV_MA(psi,U,Sig,h)
[T,n] = size(U);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T); 
Utld = Hpsi\U;
CSig = chol(Sig,'lower');
c = -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) - n/2*log(1+psi^2);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
llike = c -.5*(s2(1)/((1+psi^2)*exp(h(1))) + sum(s2(2:end)./exp(h(2:end))));
end
