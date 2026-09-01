% This function constructs the conditional prior of the VAR coefficients

function Vbeta = getVbeta(idx_kappa1,idx_kappa2,kappa,C,sig2)
np = length(idx_kappa1);
k_beta = length(C);
Vbeta = zeros(k_beta,1);
Vbeta(1:np+1:end) = kappa(4)*sig2;          % intercepts
Vbeta(idx_kappa1) = kappa(1)*C(idx_kappa1); % own lags
Vbeta(idx_kappa2) = kappa(2)*C(idx_kappa2); % other lags
end