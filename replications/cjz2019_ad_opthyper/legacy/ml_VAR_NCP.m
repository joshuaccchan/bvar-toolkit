% This function evaluates the marginal likelihood under the natural 
% conjugate prior
%
% See:
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019.

function lml = ml_VAR_NCP(VA,S0,nu0,KA,S_hat,T)
    n = size(S0,1);
    lml = -n*T/2*log(pi) - n/2*(ldet(VA) + ldet(KA)) + nu0/2*ldet(S0)...
        - (nu0+T)/2*ldet(S_hat) + mgammaln(n,(nu0+T)/2) - mgammaln(n,nu0/2);
end