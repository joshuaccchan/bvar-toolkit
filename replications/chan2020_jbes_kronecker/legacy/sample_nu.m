% This script samples the degree of freedom parameter in the t models
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

function [nu,flag,fnu] = sample_nu(lam,nu,nuub)
flag = 0;
T = size(lam,1);
sum1 = sum(log(lam));
sum2 = sum(1./lam);
fnu = @(x) T*(x/2.*log(x/2)-gammaln(x/2)) - (x/2+1)*sum1 - x/2*sum2;
f1 = @(x) T/2*(log(x/2)+1-psi(x/2)) - .5*(sum1+sum2);
f2 = @(x) T/(2*x) - T/4*psi(1,x/2);
err_nu = 1;
nut = nu;
while abs(err_nu) > 10^(-5) % stopping criteria
    err_nu = f1(nut);
    Knu = -f2(nut);         % infomation matrix
    nut = nut + Knu\err_nu;
    if nut<2
        nut = 5;
        Knu = -f2(nut);
        break;
    end
end
sqrtDnu = sqrt(1/Knu);
nuc = nut + sqrtDnu*randn; 
if nuc>2 && nuc < nuub
    alp_MH = exp(fnu(nuc)-fnu(nu)) ... 
        * normpdf(nu,nut,sqrtDnu)/normpdf(nuc,nut,sqrtDnu);
    if alp_MH > rand
        nu = nuc;
        flag = 1;
    end    
end    
    
end