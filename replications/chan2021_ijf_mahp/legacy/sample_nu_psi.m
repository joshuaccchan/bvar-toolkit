% This script samples the parameter nu_psi in the Minnesota-type 
% normal-gamma prior
% 
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for 
% Large Bayesian VARs, International Journal of Forecasting, forthcoming

function [nu_psi, flag, f_nu_psi] = sample_nu_psi(psi_kappa1,psi_kappa2,...
    nu_psi,lam0_nu_psi)
flag = 0;
n2p = size(psi_kappa1,1) + size(psi_kappa2,1);
sum1 = sum(log(psi_kappa1)) + sum(log(psi_kappa2));
sum2 = sum(psi_kappa1) + sum(psi_kappa2);
f_nu_psi = @(x) n2p*(x*log(x/2) - gammaln(x)) + (x-1)*sum1 - x/2*sum2 - lam0_nu_psi*x;
df_nu_psi = @(x) n2p*(log(x/2) + 1 - psi(x)) + sum1 - .5*sum2 - lam0_nu_psi;
d2f_nu_psi = @(x) n2p/x - n2p*psi(1,x);
S_nu_psi = 1;
nu_psit = nu_psi;
count = 0;
while abs(S_nu_psi) > 1e-5 && count < 100  % stopping criteria
    S_nu_psi = df_nu_psi(nu_psit);
    H_nu_psi = d2f_nu_psi(nu_psit); 
    nu_psit = nu_psit - H_nu_psi\S_nu_psi;
    if nu_psit < 0
        nu_psit = fminbnd(@(x)-f_nu_psi(x),.01,5);
        H_nu_psi = d2f_nu_psi(nu_psit);
        break;
    end
end
if count == 100
    nu_psit = fminbnd(@(x)-f_nu_psi(x),.01,5);
    H_nu_psi = d2f_nu_psi(nu_psit);
end
        
Dnu_psi = -1/H_nu_psi;
nu_psic = nu_psit + sqrt(Dnu_psi)*randn; 
if nu_psic > 0
    alp_MH = f_nu_psi(nu_psic) - f_nu_psi(nu_psi) ...
        - .5*(nu_psi-nu_psit)^2/Dnu_psi + .5*(nu_psic-nu_psit)^2/Dnu_psi;
    if exp(alp_MH) > rand
        nu_psi = nu_psic;
        flag = 1;
    end    
end
end