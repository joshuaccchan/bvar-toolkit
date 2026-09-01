% This script samples the df parameter
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

function [nu, flag, f_nu] = sample_nu(lam,nu,nu_ub)
flag = 0;
T = size(lam,1);
sum1 = sum(log(lam));
sum2 = sum(1./lam);
f_nu = @(x) T*(x/2.*log(x/2)-gammaln(x/2)) - (x/2+1)*sum1 - x/2*sum2;
T = size(lam,1);
sum1 = sum(log(lam));
sum2 = sum(1./lam);
df_nu = @(x) T/2*(log(x/2) + 1 - psi(x/2)) - .5*(sum1+sum2);
d2f_nu = @(x) T/(2*x) - T/4*psi(1,x/2);
S_nu = 1;
nut = nu;
while abs(S_nu) > 10^(-5)   % stopping criteria
    S_nu = df_nu(nut);
    H_nu = d2f_nu(nut); 
    nut = nut - H_nu\S_nu;
    if nut<2
        nut = 5;
        H_nu = d2f_nu(nut);
        break;
    end
end
Dnu = -1/H_nu;
nuc = nut + sqrt(Dnu)*randn; 
if nuc > 2 && nuc < nu_ub
    lalp_MH = f_nu(nuc) - f_nu(nu) - .5*(nu-nut)^2/Dnu + .5*(nuc-nut)^2/Dnu;
    if exp(lalp_MH) > rand
        nu = nuc;
        flag = 1;
    end    
end
end