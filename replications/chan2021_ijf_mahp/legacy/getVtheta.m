% This function constructs the conditional prior of the VAR coefficients
% and the impact matrix
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for 
% Large Bayesian VARs, International Journal of Forecasting, forthcoming

function [Valp,Vbeta] = getVtheta(idx_kappa1,idx_kappa2,kappa,C,sig2)
np = length(idx_kappa1);
n = length(idx_kappa2)/np + 1;
k_beta = length(C);
k_alp = n*(n-1)/2;
Vbeta = zeros(k_beta,1);
Valp = zeros(k_alp,1);

Vbeta(1:np+1:end) = kappa(4)*sig2;          % intercepts
Vbeta(idx_kappa1) = kappa(1)*C(idx_kappa1); % own lags
Vbeta(idx_kappa2) = kappa(2)*C(idx_kappa2); % other lags

count_alp = 0;
for ii = 1:n
    Valpi = kappa(3)*repmat(sig2(ii),ii-1,1)./sig2(1:ii-1);    
    Valp(count_alp+1:count_alp+ii-1) = Valpi;    
    count_alp = count_alp + ii - 1;
end
end